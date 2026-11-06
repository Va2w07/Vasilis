"""
RobustExtractor.py — drop-in replacement for the Bayesian -> Adam pipeline.

Same interface as LayeredExtractor / BayesianLayeredExtractor:
    RobustLayeredExtractor(reference_pulse, experimental_pulse, deltat, layers_init, optimize_mask=None)
    .fit(...) -> [(complex n_j, D_j), ...]      (same return format as LayeredExtractor.optimize)

What it does differently (and why — see the 5-layer diagnosis notes):
  * float64 forward model (identical formulae / sign convention to Simulate.py; k < 0 = loss),
    so the loss floor is ~1e-10 instead of the float32 ~1e-6.
  * scaled parameters: n = n0 + 0.1 x,  k = k0 + 0.005 x,  D = D0 exp(0.1 x)
    -> every coordinate moves a physically sensible amount per optimiser step
       (raw Adam with lr=1e-3 moves k by 0.1 per 100 iterations, i.e. 10-100x its physical size).
  * hinge penalty  lam * sum(relu(k))  forbids gain (k > 0), which the unconstrained fit used
    as free per-layer amplitude knobs (k3 = -0.26, k4 = +0.40 in the stuck run).
  * coarse-to-fine: L-BFGS on a Gaussian low-passed version of both pulses
    (sigma_f = 0.25, 0.5, 1 THz, then full band). The smoothed loss depends mainly on the
    total optical path, so the echo-assignment local minima (wrong per-layer split of n*d)
    are not there yet; each stage then starts in the right basin of the next.
  * optional multi-start around layers_init, keep the best. 12/12 random starts
    (n +-0.15, D +-35 %) reached the exact truth on the 5-layer test; plain full-band L-BFGS
    from the same starts: 3/12.

Typical use in the notebook (replaces the Bayesian + Adam cells):
    RE = RobustLayeredExtractor(reference_pulse, experimental_pulse, deltat, layers_init,
                                optimize_mask=optimization_mask)
    optim_params = RE.fit(n_starts=6)          # ~10-20 s per start on CPU
    optim_params = RE.fit(n_starts=1)          # single start from layers_init
"""
import numpy as np
import torch

c = 299792458.0
DT = torch.float64
CT = torch.complex128


# ----------------------------------------------------------------------------- forward model
def _rts(n0, nj, Dj):
    cj = torch.cos(nj * Dj)
    sj = torch.sin(nj * Dj)
    d = cj + 0.5j * (nj / n0 + n0 / nj) * sj
    r = 0.5j * sj * (n0 / nj - nj / n0) / d
    t = 1.0 / d
    return r, t * torch.exp(1j * n0 * Dj), t * t - r * r


def _transmission(n_c, kD):
    F, m = kD.shape
    n0 = torch.tensor(1.0, dtype=CT)
    U = torch.zeros(F, dtype=CT)
    V = torch.ones(F, dtype=CT)
    T = torch.ones(F, dtype=CT)
    for j in range(m):
        r, t, s = _rts(n0, n_c[j].expand(F), kD[:, j].to(CT))
        Vlast = V
        U, V = r * V + s * U, V - r * U
        T = T * t * Vlast / V
    return T


class Forward64:
    """y = forward(n, k, d) for real (m,) tensors; same zero-padding / spectrum construction as simulate_parallel."""

    def __init__(self, x, deltat):
        self.x = torch.as_tensor(x).detach().to(DT)
        self.L = len(self.x)
        self.M = 2 * self.L
        self.N = 4 * self.L
        deltaf = 1.0 / (self.N * deltat)
        self.dk = 2 * torch.pi * deltaf / c
        self.kvals = torch.arange(self.M + 1, dtype=DT)
        z = torch.zeros(self.N, dtype=DT)
        z[: self.L] = self.x
        self.X = torch.fft.fft(z) / self.N
        self.freqs = torch.fft.fftfreq(self.N, d=deltat)

    def T_full(self, n, k, d):
        n_c = n.to(CT) + 1j * k.to(CT)
        kD = self.dk * self.kvals[:, None] * d[None, :]
        T_half = _transmission(n_c, kD)
        return torch.cat([T_half, torch.conj(torch.flip(T_half[1: self.M], dims=[0]))])

    def __call__(self, n, k, d, H=None):
        Y = self.T_full(n, k, d) * self.X
        if H is not None:
            Y = Y * H
        return (self.N * torch.fft.ifft(Y).real)[: self.L]


# ----------------------------------------------------------------------------- extractor
class RobustLayeredExtractor:
    def __init__(self, reference_pulse, experimental_pulse, deltat, layers_init, optimize_mask=None,
                 n_scale=0.1, k_scale=0.005, d_scale=0.1, lam_hinge=1.0):
        self.fwd = Forward64(reference_pulse, deltat)
        self.y_exp = torch.as_tensor(experimental_pulse).detach().to(DT)[: self.fwd.L]
        self.deltat = deltat
        self.m = len(layers_init)
        self.n0 = torch.tensor([float(np.real(l[0])) for l in layers_init], dtype=DT)
        self.k0 = torch.tensor([float(np.imag(l[0])) for l in layers_init], dtype=DT)
        self.d0 = torch.tensor([float(l[1]) for l in layers_init], dtype=DT)
        self.mask = optimize_mask or [(True, True, True)] * self.m
        self.scales = (n_scale, k_scale, d_scale)
        self.lam = lam_hinge
        # boolean selector of the free coordinates in the full (3m,) x-vector
        free = [m_[0] for m_ in self.mask] + [m_[1] for m_ in self.mask] + [m_[2] for m_ in self.mask]
        self.free = torch.tensor(free, dtype=torch.bool)
        self.n_free = int(self.free.sum())
        self.loss_history = []

    # --- parametrisation ---------------------------------------------------
    def _unpack(self, x_free, n0=None, k0=None, d0=None):
        n0 = self.n0 if n0 is None else n0
        k0 = self.k0 if k0 is None else k0
        d0 = self.d0 if d0 is None else d0
        x = torch.zeros(3 * self.m, dtype=DT)
        x = x.masked_scatter(self.free, x_free)
        m = self.m
        sn, sk, sd = self.scales
        n = torch.clamp(n0 + sn * x[:m], min=1.0)
        k = k0 + sk * x[m: 2 * m]
        d = d0 * torch.exp(torch.clamp(sd * x[2 * m:], -3.0, 3.0))
        return n, k, d

    def _pack(self, n, k, d):
        sn, sk, sd = self.scales
        n = torch.as_tensor(n, dtype=DT); k = torch.as_tensor(k, dtype=DT); d = torch.as_tensor(d, dtype=DT)
        x = torch.cat([(n - self.n0) / sn, (k - self.k0) / sk, torch.log(d / self.d0) / sd])
        return x[self.free]

    # --- loss --------------------------------------------------------------
    def _set_lowpass(self, f_lp):
        if f_lp is None:
            self.H, self.y_target = None, self.y_exp
        else:
            self.H = torch.exp(-0.5 * (self.fwd.freqs / f_lp) ** 2)
            z = torch.zeros(self.fwd.N, dtype=DT); z[: self.fwd.L] = self.y_exp
            self.y_target = torch.fft.ifft(torch.fft.fft(z) * self.H).real[: self.fwd.L]

    def _loss(self, x_free):
        n, k, d = self._unpack(x_free)
        y = self.fwd(n, k, d, self.H)
        data = torch.sqrt(torch.mean((y - self.y_target) ** 2))
        return data + self.lam * torch.relu(k).sum(), data

    def loss_at(self, layers):
        """rMSE (full band, no penalty) of an arbitrary [(complex n, D), ...] list."""
        n = torch.tensor([float(np.real(l[0])) for l in layers], dtype=DT)
        k = torch.tensor([float(np.imag(l[0])) for l in layers], dtype=DT)
        d = torch.tensor([float(l[1]) for l in layers], dtype=DT)
        with torch.no_grad():
            return torch.sqrt(torch.mean((self.fwd(n, k, d) - self.y_exp) ** 2)).item()

    # --- optimisers ----------------------------------------------------------
    def _lbfgs(self, x0, max_iter):
        x = x0.clone().detach().requires_grad_(True)
        opt = torch.optim.LBFGS([x], lr=1.0, max_iter=max_iter, history_size=50, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-12, tolerance_change=1e-12)
        hist = []

        def closure():
            opt.zero_grad()
            loss, data = self._loss(x)
            loss.backward()
            hist.append(data.item())
            return loss
        opt.step(closure)
        with torch.no_grad():
            loss, _ = self._loss(x)
        if not torch.isfinite(loss):
            return x0.clone(), float('inf'), hist
        return x.detach(), loss.item(), hist

    def fit_from(self, x0, schedule=(0.25e12, 0.5e12, 1.0e12, None), max_iter=300, verbose=True):
        """coarse-to-fine L-BFGS from a free-coordinate start vector."""
        x = x0.clone()
        for f_lp in schedule:
            self._set_lowpass(f_lp)
            x, l, hist = self._lbfgs(x, max_iter)
            self.loss_history += hist
            if verbose:
                tag = 'full band' if f_lp is None else f'low-pass {f_lp/1e12:.2f} THz'
                print(f"    {tag:22s} -> loss {l:.3e}  ({len(hist)} evaluations)")
        self._set_lowpass(None)
        return x, l

    def fit(self, n_starts=1, n_jitter=0.15, d_jitter=0.35, k_jitter=0.01, seed=0,
            schedule=(0.25e12, 0.5e12, 1.0e12, None), max_iter=300, adam_polish=0, verbose=True):
        """
        n_starts = 1 : single coarse-to-fine fit from layers_init.
        n_starts > 1 : start 0 = layers_init, the rest = random perturbations
                       (n +- n_jitter, D * exp(+- d_jitter), k in [-k_jitter, 0]); best loss wins.
        adam_polish  : optional extra Adam iterations (lr 1e-3 in scaled coordinates) on the winner.
        """
        rng = np.random.default_rng(seed)
        starts = [self._pack(self.n0, self.k0, self.d0)]
        for _ in range(n_starts - 1):
            n = self.n0 + torch.tensor(rng.uniform(-n_jitter, n_jitter, self.m))
            d = self.d0 * torch.exp(torch.tensor(rng.uniform(-d_jitter, d_jitter, self.m)))
            k = torch.tensor(rng.uniform(-k_jitter, 0.0, self.m))
            starts.append(self._pack(n, k, d))
        print(f"Fitting {self.n_free} free parameters, {n_starts} start(s), coarse-to-fine L-BFGS.")
        best = (float('inf'), None)
        self.start_results = []
        for s, x0 in enumerate(starts):
            if verbose:
                print(f"  start {s}:")
            x, l = self.fit_from(x0, schedule, max_iter, verbose)
            self.start_results.append((l, self._to_layers(x)))
            if verbose:
                print(f"  start {s}: loss {l:.3e}")
            if l < best[0]:
                best = (l, x)
        x = best[1]
        if adam_polish:
            x = x.clone().requires_grad_(True)
            opt = torch.optim.Adam([x], lr=1e-3)
            for _ in range(adam_polish):
                opt.zero_grad()
                loss, data = self._loss(x)
                loss.backward()
                self.loss_history.append(data.item())
                opt.step()
            x = x.detach()
        self.best_loss = self._loss(x)[1].item()
        self.best_x = x
        print(f"best loss (rMSE, full band): {self.best_loss:.3e}")
        return self._to_layers(x)

    def _to_layers(self, x):
        n, k, d = [t.detach() for t in self._unpack(x)]
        return [(complex(n[i].item(), k[i].item()), d[i].item()) for i in range(self.m)]

    def simulate(self, layers=None):
        """time-domain pulse (numpy) for the fitted (or given) layers — for plotting."""
        if layers is None:
            n, k, d = [t.detach() for t in self._unpack(self.best_x)]
        else:
            n = torch.tensor([float(np.real(l[0])) for l in layers], dtype=DT)
            k = torch.tensor([float(np.imag(l[0])) for l in layers], dtype=DT)
            d = torch.tensor([float(l[1]) for l in layers], dtype=DT)
        with torch.no_grad():
            return self.fwd(n, k, d).numpy()
