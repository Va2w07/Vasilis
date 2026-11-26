"""
quartzsi_model.py — the quartz/Si two-layer problem shared by nb3 (Adam + landscape) and nb4 (L-BFGS).

Reproduces the protocol of the reference notebook `QuartzSI_loss_landscape_Adam_visualization.ipynb`:
    * measured reference pulse propagated through the two-layer TMM (same formulae as Simulate.py),
      here in float64 (`RobustExtractor.Forward64`);
    * parameter vector p = (n1, κ1, d1, n2, κ2, d2), complex index ñ = n − iκ (κ > 0 = loss);
    * loss = MAE + log-cosh of the time-domain residual;
    * the optimiser works on multiplicative factors u = p / p_start (all start at 1).
"""
import numpy as np
import torch
from RobustExtractor import Forward64

DT = torch.float64
PARAM_NAMES = ['n1', 'k1', 'd1', 'n2', 'k2', 'd2']
PARAM_TEX = {'n1': r'$n_1$', 'k1': r'$\kappa_1$', 'd1': r'$d_1$', 'n2': r'$n_2$', 'k2': r'$\kappa_2$', 'd2': r'$d_2$'}

# the three starts of the reference notebook: quartz n1 = 1.982 / 1.900 / 2.010, d1 = 1000 / 900 / 900 µm,
# κ1 = 0.0071; silicon always n2 = 3.40, κ2 = 0.0054, d2 = 400 µm
STARTS = np.array([[1.9820, 0.0071, 1.0e-3, 3.4000, 0.0054, 0.4e-3],
                   [1.9000, 0.0071, 0.9e-3, 3.4000, 0.0054, 0.4e-3],
                   [2.0100, 0.0071, 0.9e-3, 3.4000, 0.0054, 0.4e-3]])


def log_cosh_loss(y_pred, y_test):
    diff = y_pred - y_test
    return torch.mean(torch.log(torch.cosh(diff + 1e-12)))


class QuartzSiProblem:
    def __init__(self, ref_file='../data/Reference_15min_nitrogen.txt', sample_file='../data/Quartz_Si_15min_nitrogen.txt',
                 prepulse_ps=10.0):
        data_ref = np.loadtxt(ref_file); data_sam = np.loadtxt(sample_file)
        t = data_ref[:, 0]
        self.t_r = t - t[0]
        self.deltat = (t[1] - t[0]) * 1e-12
        self.x_exp = data_ref[:, 1].astype(np.float64)       # reference pulse
        self.true_x_exp = data_sam[:, 1].astype(np.float64)  # sample pulse
        self.L = len(self.x_exp)
        self.fwd = Forward64(self.x_exp, self.deltat)
        self.y_true = torch.tensor(self.true_x_exp, dtype=DT)
        self.prepulse_ps = prepulse_ps
        self.sigma = float(np.std(self.true_x_exp[self.t_r < prepulse_ps]))   # pre-pulse std of the sample trace
        self.f_THz = np.fft.rfftfreq(self.L, self.deltat) / 1e12

    # --- forward model -------------------------------------------------------
    def ysim(self, p):
        """p: tensor [6] = (n1, κ1, d1, n2, κ2, d2) -> simulated sample pulse [L] (float64, differentiable)."""
        p = torch.as_tensor(p, dtype=DT)
        return self.fwd(p[[0, 3]], -p[[1, 4]], p[[2, 5]])      # Forward64 uses n + i k with k = −κ

    def ysim_np(self, p):
        with torch.no_grad():
            return self.ysim(torch.as_tensor(np.asarray(p, float), dtype=DT)).numpy()

    # --- losses --------------------------------------------------------------
    def loss_from_y(self, y, kind='mae_logcosh'):
        r = self.y_true - y
        if kind == 'mae_logcosh':
            return r.abs().mean() + log_cosh_loss(y, self.y_true)
        if kind == 'mae':
            return r.abs().mean()
        if kind == 'rmse':
            return torch.sqrt(torch.mean(r ** 2))
        if kind == 'logcosh':
            return log_cosh_loss(y, self.y_true)
        raise ValueError(kind)

    def loss_torch(self, p, kind='mae_logcosh'):
        return self.loss_from_y(self.ysim(p), kind)

    def loss_np(self, p, kind='mae_logcosh'):
        with torch.no_grad():
            return self.loss_torch(torch.as_tensor(np.asarray(p, float), dtype=DT), kind).item()

    def mae_np(self, p):
        return float(np.mean(np.abs(self.true_x_exp - self.ysim_np(p))))

    # --- sensitivity / Cramér–Rao --------------------------------------------
    def jacobian(self, ref, scale):
        """J = d y_sim / d x  [L, 6] in scaled coordinates x = (p − ref)/scale, float64 autograd."""
        ref_t = torch.tensor(np.asarray(ref, float), dtype=DT); sc_t = torch.tensor(np.asarray(scale, float), dtype=DT)
        f = lambda x: self.ysim(ref_t + x * sc_t)
        return torch.autograd.functional.jacobian(f, torch.zeros(6, dtype=DT)).detach().numpy()
