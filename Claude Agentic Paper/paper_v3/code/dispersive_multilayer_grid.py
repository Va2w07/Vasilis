"""
N-layer THz-TDS transfer-matrix model with COARSE-GRID n(f), k(f): each curve is given by its values at N_KNOTS knot
frequencies inside the fit band and interpolated in between (natural cubic spline or linear). No polynomial.

Companion to dispersive_multilayer.py (polynomial version); same per-layer freeze flags, so an air gap can be held at
n = 1, k = 0 while its thickness is fitted. Same forward model (`Grid`, `simulate_dispersive`) and loss / `fit` loop.

    n_j(f) = sum_q A_q(f) n_{j,q},   k_j(f) = sum_q A_q(f) k_{j,q}

A_q(f) are fixed basis functions obtained by interpolating the unit vectors (so any interpolation that is linear in the
knot values works). Outside [F_LO, F_HI] the frequency is clipped -> the curve is held at its edge value.

Stage 1 ("flat"): knots of a curve tied together, only one common offset per curve moves.
Stage 2 ("dispersive"): every knot free.
"""
import numpy as np
import torch
from scipy.interpolate import CubicSpline

from dispersive_multilayer import Grid, simulate_dispersive, data_loss, total_loss   # noqa: F401 (re-exported)


def interp_basis(f, f_knots, kind='cubic'):
    """Basis matrix A [len(f), N_KNOTS]: curve(f) = A @ values_at_knots. f clipped to the knot range."""
    f = np.clip(np.asarray(f, dtype=np.float64), f_knots[0], f_knots[-1])
    A = np.zeros((len(f), len(f_knots)))
    for q in range(len(f_knots)):
        e = np.zeros(len(f_knots)); e[q] = 1.0
        if kind == 'linear':
            A[:, q] = np.interp(f, f_knots, e)
        elif kind == 'cubic':
            A[:, q] = CubicSpline(f_knots, e, bc_type='natural')(f)
        else:
            raise ValueError(kind)
    return A


class DispersiveMultiLayerGrid:
    """
    layer_specs: list of dicts (beam order): name, n, k, d (initial values; n, k are used for every knot),
                 fit_n / fit_k / fit_d (False -> frozen, e.g. air), optional d_scale.
    """

    def __init__(self, grid, layer_specs, n_knots=6, interp='cubic',
                 n_scale=0.05,      # absolute Adam step scale for n knots (moves ~ lr * scale per iteration)
                 k_scale=0.005):    # absolute step scale for k knots
        self.grid = grid
        self.specs = [dict(s) for s in layer_specs]
        self.m = len(self.specs)
        self.names = [s['name'] for s in self.specs]
        self.N_KNOTS, self.interp = n_knots, interp
        self.f_knots = np.linspace(grid.F_LO, grid.F_HI, n_knots)
        self.A_grid = torch.as_tensor(interp_basis(grid.f_np, self.f_knots, interp), dtype=torch.float32)   # [M+1, Q]

        self.n_init = torch.tensor([[s['n']] * n_knots for s in self.specs])        # [m, Q]
        self.k_init = torch.tensor([[s['k']] * n_knots for s in self.specs])
        self.d_init = torch.tensor([s['d'] for s in self.specs], dtype=torch.float32)

        self.n_scale = torch.full((self.m, n_knots), float(n_scale))
        self.k_scale = torch.full((self.m, n_knots), float(k_scale))
        self.d_scale = torch.tensor([s.get('d_scale') or s['d'] for s in self.specs], dtype=torch.float32)

        # trainable: one common offset per curve (stage 1) + individual knot offsets (stage 2)
        self.n_off = torch.zeros(self.m, requires_grad=True)
        self.k_off = torch.zeros(self.m, requires_grad=True)
        self.n_raw = torch.zeros(self.m, n_knots, requires_grad=True)
        self.k_raw = torch.zeros(self.m, n_knots, requires_grad=True)
        self.d_raw = torch.zeros(self.m, requires_grad=True)

        self.fit_n = torch.tensor([float(s.get('fit_n', True)) for s in self.specs])
        self.fit_k = torch.tensor([float(s.get('fit_k', True)) for s in self.specs])
        self.fit_d = torch.tensor([float(s.get('fit_d', True)) for s in self.specs])

        self.knots_free = 0.0       # 0 -> knots tied (flat curves), 1 -> every knot free
        self.k_free = 1.0
        self.d_mask = torch.ones(self.m)
        self.band = (grid.f_THz >= grid.F_LO) & (grid.f_THz <= grid.F_HI)

    # ---- configuration -------------------------------------------------------------------------------------
    def parameters(self):
        return [self.n_off, self.k_off, self.n_raw, self.k_raw, self.d_raw]

    def set_dispersion(self, free: bool):
        self.knots_free = 1.0 if free else 0.0

    def set_thickness(self, free: bool):
        self.d_mask[:] = 1.0 if free else 0.0

    def set_absorption(self, free: bool):
        self.k_free = 1.0 if free else 0.0

    # ---- physical quantities -------------------------------------------------------------------------------
    def n_knots(self):
        return self.n_init + self.fit_n[:, None] * self.n_scale * (self.n_off[:, None] + self.knots_free * self.n_raw)

    def k_knots(self):
        return self.k_init + self.fit_k[:, None] * self.k_free * self.k_scale * (self.k_off[:, None] + self.knots_free * self.k_raw)

    def d(self):
        return self.d_init + self.d_raw * self.d_scale * self.d_mask * self.fit_d

    def nk(self, f=None):
        A = self.A_grid if f is None else torch.as_tensor(interp_basis(np.asarray(f), self.f_knots, self.interp), dtype=torch.float32)
        return self.n_knots() @ A.T, self.k_knots() @ A.T

    def layers(self):
        n, k = self.nk()
        d = self.d()
        return [(n[j] - 1j * k[j], d[j]) for j in range(self.m)]

    def forward(self, x):
        T, y = simulate_dispersive(self.grid, x, self.layers())
        return T, y[:self.grid.L]

    def k_negativity(self):
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
            return dict(n=self.n_knots().clone(), k=self.k_knots().clone(), d=self.d().clone())

    def load_state(self, st):
        with torch.no_grad():
            self.n_init = st['n'].clone(); self.k_init = st['k'].clone(); self.d_init = st['d'].clone()
            self.n_off.zero_(); self.k_off.zero_(); self.n_raw.zero_(); self.k_raw.zero_(); self.d_raw.zero_()
            self.d_scale = torch.where(self.d_init > 0, self.d_init, self.d_scale)

    def report(self):
        st = self.state()
        neg = self.k_negative_fraction()
        print('knots (THz): ' + ', '.join(f'{f:.3f}' for f in self.f_knots))
        for j, nm in enumerate(self.names):
            n = ', '.join(f"{v:.5f}" for v in st['n'][j].tolist())
            k = ', '.join(f"{v:+.5f}" for v in st['k'][j].tolist())
            tag = '' if (self.fit_n[j] > 0 or self.fit_k[j] > 0) else '   (n, k fixed)'
            print(f"{nm:>8s}:  n knots = [{n}]{tag}")
            print(f"           k knots = [{k}]      d = {st['d'][j].item()*1e6:8.2f} um   k<0 in band: {neg[j]*100:.1f} %")


def fit(model, x_ref, y_meas, n_iters=20000, lr=1e-4, print_step=2000, label='', k_pos_weight=1.0, progress=True):
    """Adam on the raw parameters. Histories: n [it,m,Q], k [it,m,Q], d [it,m], loss, penalty."""
    try:
        from tqdm.auto import tqdm
    except ImportError:
        tqdm = lambda it, **kw: it
    y_meas = torch.as_tensor(y_meas, dtype=torch.float32)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    hist = {'n': [], 'k': [], 'd': [], 'loss': [], 'penalty': []}
    q1 = int(np.argmin(np.abs(model.f_knots - 1.0)))      # knot closest to 1 THz for the progress line

    def record(loss):
        st = model.state()
        hist['n'].append(st['n'].numpy()); hist['k'].append(st['k'].numpy()); hist['d'].append(st['d'].numpy())
        hist['loss'].append(loss.item())
        with torch.no_grad():
            hist['penalty'].append(k_pos_weight * model.k_negativity().item())

    def line(i, loss):
        st = model.state()
        s = ' '.join(f"{nm}: n={st['n'][j,q1]:.5f} k={st['k'][j,q1]:.5f} d={st['d'][j]*1e6:.1f}um |"
                     for j, nm in enumerate(model.names))
        print(f"[{label}] Epoch={i:6d}  loss={loss.item():.4e} (k-pen {hist['penalty'][-1]:.1e})  [knot {model.f_knots[q1]:.2f} THz] {s}")

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
