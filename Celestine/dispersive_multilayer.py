"""
N-layer THz-TDS transfer-matrix model with low-order polynomial n(f), k(f) per layer, fitted with Adam.

Generalisation of the two-layer `DispersiveTwoLayer` of QuartzSI_dispersive_Adam.ipynb to any number of layers,
with per-layer freeze flags so that e.g. an air gap can be kept at n = 1, k = 0 while its thickness is fitted.

    n_j(f) = c0 + a x + b x^2 + c x^3        (4 coefficients)
    k_j(f) = c0 + a x + b x^2                (3 coefficients)
    x      = (clip(f, F_LO, F_HI) - F_REF) / F_SCALE

Every physical coefficient is written as  init + raw * scale * mask  with the trainable `raw` starting at 0.
`Simulate.py` is imported unchanged (only RTm_batched is used).
"""
import numpy as np
import torch

from Simulate import RTm_batched

c_light = 299792458.0

N_COEF_NAMES = ['c0', 'a', 'b', 'c']     # n: cubic
K_COEF_NAMES = ['c0', 'a', 'b']          # k: quadratic
N_ORDER, K_ORDER = len(N_COEF_NAMES), len(K_COEF_NAMES)


# ---------------------------------------------------------------------------------------------------------------
# frequency grid / normalised frequency
# ---------------------------------------------------------------------------------------------------------------
class Grid:
    """FFT grid of the forward model (same convention as Simulate.simulate_parallel: M = 2L, N = 4L)."""

    def __init__(self, L, deltat, f_ref=1.0, f_scale=1.0, f_lo=0.2, f_hi=3.0):
        self.L, self.deltat = L, deltat
        self.M, self.N = 2 * L, 4 * L
        self.deltaf = 1.0 / (self.N * deltat)
        self.dk = 2 * np.pi * self.deltaf / c_light
        self.k_vals = torch.arange(self.M + 1, dtype=torch.float32)
        self.f_THz = self.k_vals * self.deltaf / 1e12
        self.f_np = self.f_THz.numpy()
        self.F_REF, self.F_SCALE, self.F_LO, self.F_HI = f_ref, f_scale, f_lo, f_hi

    def xnorm(self, f):
        if isinstance(f, torch.Tensor):
            return (torch.clamp(f, self.F_LO, self.F_HI) - self.F_REF) / self.F_SCALE
        return (np.clip(f, self.F_LO, self.F_HI) - self.F_REF) / self.F_SCALE


def poly(coefs, x):
    """coefs[0] + coefs[1] x + coefs[2] x^2 + ...  (Horner)."""
    out = coefs[len(coefs) - 1]
    for i in range(len(coefs) - 2, -1, -1):
        out = out * x + coefs[i]
    return out


def simulate_dispersive(grid, x, layers):
    """
    layers: list of (n_complex [M+1] or scalar, d) ; d may carry a gradient.
    Returns (T [N] complex, y [N] real); y[:L] is the simulated sample pulse.
    """
    M, N, L = grid.M, grid.N, grid.L
    xt = torch.as_tensor(x, dtype=torch.float32)
    device = xt.device
    n0 = torch.tensor(1.0, dtype=torch.cfloat, device=device)

    batched = []
    for nj, d in layers:
        nj = torch.as_tensor(nj, dtype=torch.cfloat, device=device)
        if nj.ndim == 0:
            nj = nj.expand(M + 1)
        batched.append((nj, grid.dk * grid.k_vals * d))

    _, T_half = RTm_batched(len(batched), n0, batched)

    T = torch.zeros(N, dtype=torch.cfloat, device=device)
    T[:M + 1] = T_half
    T[M + 1:] = torch.conj(torch.flip(T[1:M], dims=[0]))

    z = torch.zeros(N, dtype=torch.float32, device=device)
    z[:L] = xt
    X = torch.fft.fft(z) / N
    y = N * torch.fft.ifft(T * X).real
    return T, y


# ---------------------------------------------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------------------------------------------
class DispersiveMultiLayer:
    """
    layer_specs: list of dicts, one per layer, in beam order:
        name   : label
        n, k, d: initial values (n, k at F_REF; d in m)
        fit_n, fit_k, fit_d : bool  (False -> frozen at its init value, e.g. air: fit_n=False, fit_k=False)
        n_scale, k_scale, d_scale : optional per-layer scales for c0 (default: the init value itself, i.e. a relative step)
    """

    def __init__(self, grid, layer_specs,
                 n_scale_disp=(0.1, 0.1, 0.1),       # scales of a, b, c of n
                 k_scale_disp=(0.005, 0.005)):       # scales of a, b of k
        self.grid = grid
        self.specs = [dict(s) for s in layer_specs]
        self.m = len(self.specs)
        self.names = [s['name'] for s in self.specs]

        self.n_init = torch.tensor([[s['n'], 0., 0., 0.] for s in self.specs])
        self.k_init = torch.tensor([[s['k'], 0., 0.] for s in self.specs])
        self.d_init = torch.tensor([s['d'] for s in self.specs], dtype=torch.float32)

        self.n_scale = torch.tensor([[s.get('n_scale') or s['n'], *n_scale_disp] for s in self.specs])
        self.k_scale = torch.tensor([[s.get('k_scale') or (s['k'] if s['k'] != 0 else 1e-3), *k_scale_disp] for s in self.specs])
        self.d_scale = torch.tensor([s.get('d_scale') or s['d'] for s in self.specs], dtype=torch.float32)

        self.n_raw = torch.zeros(self.m, N_ORDER, requires_grad=True)
        self.k_raw = torch.zeros(self.m, K_ORDER, requires_grad=True)
        self.d_raw = torch.zeros(self.m, requires_grad=True)

        # per-layer "may this quantity ever move" flags (hard freezes, e.g. air)
        self.fit_n = torch.tensor([float(s.get('fit_n', True)) for s in self.specs])
        self.fit_k = torch.tensor([float(s.get('fit_k', True)) for s in self.specs])
        self.fit_d = torch.tensor([float(s.get('fit_d', True)) for s in self.specs])

        # stage masks (1 = free, 0 = frozen) ; combined with the hard flags in *_mask_eff()
        self.n_mask = torch.ones(self.m, N_ORDER)
        self.k_mask = torch.ones(self.m, K_ORDER)
        self.d_mask = torch.ones(self.m)

        self.xf = grid.xnorm(grid.f_THz)
        self.band = (grid.f_THz >= grid.F_LO) & (grid.f_THz <= grid.F_HI)

    # ---- configuration -------------------------------------------------------------------------------------
    def parameters(self):
        return [self.n_raw, self.k_raw, self.d_raw]

    def set_dispersion(self, free: bool):
        v = 1.0 if free else 0.0
        self.n_mask[:, 1:] = v
        self.k_mask[:, 1:] = v

    def set_thickness(self, free: bool):
        self.d_mask[:] = 1.0 if free else 0.0

    def set_absorption(self, free: bool):
        self.k_mask[:] = 1.0 if free else 0.0

    def _n_mask_eff(self): return self.n_mask * self.fit_n[:, None]
    def _k_mask_eff(self): return self.k_mask * self.fit_k[:, None]
    def _d_mask_eff(self): return self.d_mask * self.fit_d

    # ---- physical quantities -------------------------------------------------------------------------------
    def n_coefs(self): return self.n_init + self.n_raw * self.n_scale * self._n_mask_eff()
    def k_coefs(self): return self.k_init + self.k_raw * self.k_scale * self._k_mask_eff()
    def d(self):       return self.d_init + self.d_raw * self.d_scale * self._d_mask_eff()

    def nk(self, f=None):
        x = self.xf if f is None else self.grid.xnorm(torch.as_tensor(f, dtype=torch.float32))
        nc, kc = self.n_coefs(), self.k_coefs()
        n = torch.stack([poly(nc[j], x) for j in range(self.m)])
        k = torch.stack([poly(kc[j], x) for j in range(self.m)])
        return n, k

    def layers(self):
        n, k = self.nk()
        d = self.d()
        return [(n[j] - 1j * k[j], d[j]) for j in range(self.m)]

    def forward(self, x):
        T, y = simulate_dispersive(self.grid, x, self.layers())
        return T, y[:self.grid.L]

    def k_negativity(self):
        """Mean of max(0,-k(f)) over the band, over the layers whose k is fitted."""
        _, k = self.nk()
        sel = self.fit_k > 0
        return torch.relu(-k[sel][:, self.band]).mean()

    def k_negative_fraction(self):
        with torch.no_grad():
            _, k = self.nk()
            return (k[:, self.band] < 0).float().mean(dim=1).numpy()

    # ---- bookkeeping ---------------------------------------------------------------------------------------
    def state(self):
        with torch.no_grad():
            return dict(n=self.n_coefs().clone(), k=self.k_coefs().clone(), d=self.d().clone())

    def load_state(self, st):
        """Make the current physical values the new init (raw -> 0). Used between stages."""
        with torch.no_grad():
            self.n_init = st['n'].clone(); self.k_init = st['k'].clone(); self.d_init = st['d'].clone()
            self.n_raw.zero_(); self.k_raw.zero_(); self.d_raw.zero_()
            self.d_scale = torch.where(self.d_init > 0, self.d_init, self.d_scale)

    def report(self):
        st = self.state()
        neg = self.k_negative_fraction()
        for j, nm in enumerate(self.names):
            n = ', '.join(f"{v:+.5f}" for v in st['n'][j].tolist())
            k = ', '.join(f"{v:+.5f}" for v in st['k'][j].tolist())
            tag = '' if (self.fit_n[j] > 0 or self.fit_k[j] > 0) else '   (n, k fixed)'
            print(f"{nm:>8s}:  n(x) coefs [c0,a,b,c] = [{n}]{tag}")
            print(f"           k(x) coefs [c0,a,b]   = [{k}]      d = {st['d'][j].item()*1e6:8.2f} um"
                  f"   k<0 in band: {neg[j]*100:.1f} %")


# ---------------------------------------------------------------------------------------------------------------
# loss / fit
# ---------------------------------------------------------------------------------------------------------------
def log_cosh_loss(y_pred, y_test):
    diff = y_pred - y_test
    return torch.mean(torch.log(torch.cosh(diff + 1e-12)))


def data_loss(y_pred, y_test):
    return (y_test - y_pred).abs().mean() + log_cosh_loss(y_pred, y_test)


def total_loss(y_pred, y_test, model=None, k_pos_weight=1.0):
    loss = data_loss(y_pred, y_test)
    if model is not None and k_pos_weight > 0:
        loss = loss + k_pos_weight * model.k_negativity()
    return loss


def fit(model, x_ref, y_meas, n_iters=20000, lr=1e-4, print_step=2000, label='', k_pos_weight=1.0, progress=True):
    """Adam on the raw parameters. Returns per-iteration histories n [it,m,4], k [it,m,3], d [it,m], loss, penalty."""
    try:
        from tqdm.auto import tqdm
    except ImportError:
        tqdm = lambda it, **kw: it
    y_meas = torch.as_tensor(y_meas, dtype=torch.float32)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    hist = {'n': [], 'k': [], 'd': [], 'loss': [], 'penalty': []}

    def record(loss):
        st = model.state()
        hist['n'].append(st['n'].numpy()); hist['k'].append(st['k'].numpy()); hist['d'].append(st['d'].numpy())
        hist['loss'].append(loss.item())
        with torch.no_grad():
            hist['penalty'].append(k_pos_weight * model.k_negativity().item())

    def line(i, loss):
        st = model.state()
        s = ' '.join(f"{nm}: n={st['n'][j,0]:.5f} k={st['k'][j,0]:.5f} d={st['d'][j]*1e6:.1f}um |"
                     for j, nm in enumerate(model.names))
        print(f"[{label}] Epoch={i:6d}  loss={loss.item():.4e} (k-pen {hist['penalty'][-1]:.1e})  {s}")

    _, y = model.forward(x_ref)
    loss = total_loss(y, y_meas, model, k_pos_weight)
    record(loss); line(1, loss)

    it = range(1, n_iters)
    if progress:
        it = tqdm(it, desc=f'fit {label}', leave=False)
    for i in it:
        opt.zero_grad()
        loss.backward()
        opt.step()
        _, y = model.forward(x_ref)
        loss = total_loss(y, y_meas, model, k_pos_weight)
        with torch.no_grad():
            record(loss)
            if i == n_iters - 1 or (i + 1) % print_step == 0:
                line(i + 1, loss)
    return {key: np.array(v) for key, v in hist.items()}
