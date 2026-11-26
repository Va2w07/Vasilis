"""
loop1_tools.py — helpers for the loop-1 checks C1–C14 (domain review, 2026-11-24).

    noise_stats(t_ps, y)              C14: pre-pulse noise table of one trace
    TauAmpExtractor                   C1:  RobustLayeredExtractor + free reference delay tau and amplitude a
    jacobian_tau_amp(fwd, n, k, d)    C6:  Jacobian of y wrt (n, kappa, d, tau, a) in physical units
    cr_from_J(J, sigma)               Cramér–Rao covariance / bounds
    toeplitz_from_acf(...)            C5:  noise covariance matrices
    lbfgs_minimize(fun, u0)           generic float64 L-BFGS (strong Wolfe) used by C4/C5
    moving_block_bootstrap(r, block)  C5
    pulse_fwhm(t, x)                  C11a

Conventions as everywhere in the paper: ñ = n − iκ (κ > 0 = loss); Forward64 takes k = −κ.
"""
import numpy as np
import torch
from scipy.linalg import toeplitz, cho_factor, cho_solve
from RobustExtractor import RobustLayeredExtractor, Forward64

DT = torch.float64


# ----------------------------------------------------------------------------------------------- C14 noise
def noise_stats(t_ps, y, seg_ps=2.0, tmax_ps=10.0):
    """Pre-pulse noise of one trace.  Std over 0–5, 0–10, 5–10 ps; per 2-ps segment: std, mean and the
    first-difference std / sqrt(2) (the white component: a slowly varying offset cancels in the difference);
    white level = minimum and median of the segment values; autocorrelation length of the demeaned 0–10 ps
    window: first lag where the ACF drops below 1/e, and the integral time up to the first zero crossing."""
    t_ps = np.asarray(t_ps) - t_ps[0]; y = np.asarray(y, float)
    dt = t_ps[1] - t_ps[0]
    w = lambda a, b: y[(t_ps >= a) & (t_ps < b)]
    segs = []
    for a in np.arange(0, tmax_ps - 1e-9, seg_ps):
        s = w(a, a + seg_ps)
        segs.append(dict(t0_ps=float(a), t1_ps=float(a + seg_ps), std=float(np.std(s)), mean=float(np.mean(s)),
                         diff_std_over_sqrt2=float(np.std(np.diff(s)) / np.sqrt(2))))
    pre = w(0, tmax_ps); pre = pre - pre.mean()
    n = len(pre)
    acf = np.correlate(pre, pre, 'full')[n - 1:] / n; acf = acf / acf[0]
    below = np.where(acf < np.exp(-1))[0]; zero = np.where(acf <= 0)[0]
    l_e = float(below[0] * dt) if len(below) else float('nan')
    iz = zero[0] if len(zero) else n
    tau_int = float(dt * (1 + 2 * np.sum(acf[1:iz])))       # integrated correlation time (two-sided)
    white = np.array([s['diff_std_over_sqrt2'] for s in segs])
    return dict(std_0_5ps=float(np.std(w(0, 5))), std_0_10ps=float(np.std(w(0, 10))), std_5_10ps=float(np.std(w(5, 10))),
                segments=segs, white_min=float(white.min()), white_median=float(np.median(white)),
                max_abs_segment_mean=float(np.max(np.abs([s['mean'] for s in segs]))),
                acf_1_over_e_ps=l_e, acf_first_zero_ps=float(iz * dt), integrated_corr_time_ps=tau_int, dt_ps=float(dt),
                ratio_std10_over_white_min=float(np.std(w(0, 10)) / white.min()))


def prepulse_acf(y, t_ps, tmax_ps=10.0):
    """Biased autocovariance (lags 0 … n−1) of the demeaned pre-pulse; positive semi-definite as a sequence."""
    t_ps = np.asarray(t_ps) - t_ps[0]
    pre = np.asarray(y, float)[t_ps < tmax_ps]; pre = pre - pre.mean(); n = len(pre)
    return np.correlate(pre, pre, 'full')[n - 1:] / n


def toeplitz_from_acf(acf, L, taper_lags=None, nugget=0.0):
    """L x L Toeplitz covariance from an autocovariance sequence (zero beyond its length or beyond taper_lags,
    Bartlett taper if taper_lags is given), plus nugget * I."""
    r = np.zeros(L)
    m = min(len(acf), L)
    r[:m] = acf[:m]
    if taper_lags is not None:
        k = np.arange(L); r = r * np.clip(1 - k / taper_lags, 0, None)
    C = toeplitz(r)
    C[np.diag_indices(L)] += nugget
    return C


def gls_cov(J, C):
    """(J^T C^-1 J)^-1 : Cramér–Rao covariance for correlated Gaussian noise with covariance C."""
    cf = cho_factor(C, lower=True)
    return np.linalg.inv(J.T @ cho_solve(cf, J))


def sandwich_cov(J, Sigma):
    """(J^T J)^-1 J^T Sigma J (J^T J)^-1 : least-squares estimator covariance for noise covariance Sigma."""
    A = np.linalg.inv(J.T @ J)
    return A @ (J.T @ Sigma @ J) @ A


def cr_from_J(J, sigma):
    cov = sigma ** 2 * np.linalg.inv(J.T @ J)
    return cov, np.sqrt(np.diag(cov))


def moving_block_bootstrap(r, block, rng):
    L = len(r); nb = int(np.ceil(L / block))
    starts = rng.integers(0, L - block + 1, nb)
    return np.concatenate([r[s:s + block] for s in starts])[:L]


def pulse_fwhm(t_ps, x):
    """FWHM of the main field peak (linear interpolation of the half-maximum crossings) and of the
    intensity envelope |hilbert(x)|^2."""
    from scipy.signal import hilbert
    t_ps = np.asarray(t_ps); x = np.asarray(x, float)

    def fwhm(s):
        i = int(np.argmax(s)); h = s[i] / 2
        j = i
        while j > 0 and s[j] > h: j -= 1
        tl = np.interp(h, [s[j], s[j + 1]], [t_ps[j], t_ps[j + 1]])
        j = i
        while j < len(s) - 1 and s[j] > h: j += 1
        tr = np.interp(h, [s[j], s[j - 1]], [t_ps[j], t_ps[j - 1]])
        return float(tr - tl)
    sgn = 1 if x.max() >= -x.min() else -1
    return dict(field_peak_fwhm_ps=fwhm(sgn * x), intensity_envelope_fwhm_ps=fwhm(np.abs(hilbert(x)) ** 2))


# ----------------------------------------------------------------------------------------------- C1 tau / amplitude
def tau_amp_filter(freqs, tau, a):
    """a * exp(-i 2 pi f tau) on the signed FFT grid of Forward64 (positive tau = sample delayed)."""
    return a * torch.exp(-2j * torch.pi * freqs.to(torch.complex128) * tau)


class TauAmpExtractor(RobustLayeredExtractor):
    """RobustLayeredExtractor with two extra free coordinates appended to the free vector:
    tau = tau_scale * x_tau (reference/sample time-zero offset) and a = exp(a_scale * x_a) (amplitude factor).
    y = a * F^-1[ T X exp(-i w tau) ].  Low-pass stages, hinge and rMSE exactly as in the parent class."""

    def __init__(self, *args, free_tau_amp=True, tau_scale=10e-15, a_scale=0.01, **kw):
        super().__init__(*args, **kw)
        self.n_lay = self.n_free
        self.free_ta = free_tau_amp
        self.tau_scale, self.a_scale = tau_scale, a_scale

    def _unpack(self, x_free, *a, **kw):
        return super()._unpack(x_free[: self.n_lay], *a, **kw)

    def _pack(self, n, k, d, tau=0.0, amp=1.0):
        x = super()._pack(n, k, d)
        if self.free_ta:
            x = torch.cat([x, torch.tensor([tau / self.tau_scale, np.log(amp) / self.a_scale], dtype=DT)])
        return x

    def tau_amp(self, x):
        if not self.free_ta:
            return 0.0, 1.0
        return float(x[self.n_lay] * self.tau_scale), float(torch.exp(x[self.n_lay + 1] * self.a_scale))

    def _loss(self, x_free):
        n, k, d = self._unpack(x_free)
        if self.free_ta:
            H = tau_amp_filter(self.fwd.freqs, x_free[self.n_lay] * self.tau_scale, torch.exp(x_free[self.n_lay + 1] * self.a_scale))
            if self.H is not None:
                H = H * self.H
        else:
            H = self.H
        y = self.fwd(n, k, d, H)
        data = torch.sqrt(torch.mean((y - self.y_target) ** 2))
        return data + self.lam * torch.relu(k).sum(), data


def jacobian_tau_amp(fwd, n, kappa, d):
    """J [L, 3m+2] of y wrt (n_1..n_m, kappa_1..kappa_m, d_1..d_m, tau, a) at tau = 0, a = 1 (physical units)."""
    m = len(n)
    p0 = torch.tensor(np.concatenate([n, kappa, d, [0.0, 1.0]]), dtype=DT)

    def f(p):
        H = tau_amp_filter(fwd.freqs, p[3 * m], p[3 * m + 1])
        return fwd(p[:m], -p[m:2 * m], p[2 * m:3 * m], H)
    return torch.autograd.functional.jacobian(f, p0).numpy()


# ----------------------------------------------------------------------------------------------- generic L-BFGS
def lbfgs_minimize(fun, u0, max_iter=600, history=20, tol=1e-14, lr=1.0):
    """Minimise fun(u) (torch scalar, float64) with torch L-BFGS (strong Wolfe) in one .step() call, then restart
    once from the end point (fresh history) to make sure the stop is not a stale-curvature artefact.
    Returns (u_end numpy, f_end, n_evals)."""
    u = torch.tensor(np.asarray(u0, float), dtype=DT, requires_grad=True)
    nev = [0]
    for _ in range(2):
        opt = torch.optim.LBFGS([u], lr=lr, max_iter=max_iter, max_eval=int(1.25 * max_iter), history_size=history,
                                line_search_fn='strong_wolfe', tolerance_grad=1e-13, tolerance_change=tol)

        def closure():
            opt.zero_grad(); l = fun(u); l.backward(); nev[0] += 1; return l
        opt.step(closure)
    with torch.no_grad():
        fe = fun(u).item()
    return u.detach().numpy().copy(), fe, nev[0]


def grad_norm(fun, u):
    u = torch.tensor(np.asarray(u, float), dtype=DT, requires_grad=True)
    l = fun(u); l.backward()
    return float(torch.linalg.norm(u.grad)), u.grad.numpy().copy()
