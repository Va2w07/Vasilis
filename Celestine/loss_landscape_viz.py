"""
loss_landscape_viz.py
=====================
Visualising the loss landscape and the optimiser trajectory of a THz multilayer
(TMM) parameter extraction.

All functions work in *scaled* coordinates

        x_i = (p_i - p_ref_i) / scale_i

i.e. the deviation of each parameter from a reference point (normally the best
minimum found) in units of a nominal value.  This is exactly the coordinate
system Adam sees in the notebook (it optimises multiplicative factors), so
distances, PCA directions and Hessian eigenvalues are all meaningful at once,
even though n ~ 2, k ~ 0.005 and d ~ 1e-3 m differ by orders of magnitude.

Trajectory format: dict with keys
    'params' : ndarray [N, P]   parameter values per iteration
    'loss'   : ndarray [N]      loss per iteration
    'label'  : str

Typical use (see notebook):

    trajs  = load_trajectories('trajectories.npz')
    ref, scale = reference_and_scale(trajs)
    plane  = pca_plane(trajs, ref, scale)
    grid   = evaluate_plane(loss_np, plane, res=61)
    plot_pca_landscape(grid, trajs)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
from matplotlib import animation
from matplotlib.gridspec import GridSpec

# ----------------------------------------------------------------------------
# style
# ----------------------------------------------------------------------------
PARAM_NAMES = ['n1', 'k1', 'd1', 'n2', 'k2', 'd2']
PARAM_TEX = {'n1': r'$n_1$', 'k1': r'$\kappa_1$', 'd1': r'$d_1$',
             'n2': r'$n_2$', 'k2': r'$\kappa_2$', 'd2': r'$d_2$'}
# fixed categorical order (blue, orange, aqua) - colour-blind safe as a triple
TRAJ_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#4a3aa7']
LANDSCAPE_CMAP = 'Greys'       # sequential, single hue -> coloured paths stay legible
INK = '#0b0b0b'; INK2 = '#52514e'; GRID = '#d9d8d3'

plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'axes.edgecolor': INK2, 'axes.linewidth': 0.6,
    'xtick.color': INK2, 'ytick.color': INK2, 'axes.labelcolor': INK,
    'legend.frameon': False, 'figure.dpi': 110, 'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

MILESTONES = [0, 10, 100, 1000, 10000, 100000]   # iterations marked along paths


def _despine(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


# ----------------------------------------------------------------------------
# data handling
# ----------------------------------------------------------------------------
def load_trajectories(path):
    z = np.load(path, allow_pickle=True)
    names = [k[:-7] for k in z.files if k.endswith('_params')]
    return [{'params': z[f'{n}_params'], 'loss': z[f'{n}_loss'], 'label': n.replace('traj', 'Start ')}
            for n in sorted(names)]


def from_notebook_trajectory(traj_array, label):
    """Convert the notebook's np.array([n1s,k1s,d1s,n2s,k2s,d2s,losses]) layout."""
    a = np.asarray(traj_array, dtype=float)
    return {'params': a[:-1].T, 'loss': a[-1], 'label': label}


def reference_and_scale(trajs, scale=None):
    """Reference point = lowest-loss point of all trajectories.  scale = nominal magnitude per parameter."""
    best = min(trajs, key=lambda t: t['loss'].min())
    i = int(np.argmin(best['loss']))
    ref = best['params'][i].copy()
    if scale is None:
        scale = np.abs(np.vstack([t['params'][0] for t in trajs])).mean(axis=0)  # nominal = mean initial guess
    return ref, np.asarray(scale, float)


def to_scaled(p, ref, scale):   return (np.asarray(p) - ref) / scale
def from_scaled(x, ref, scale): return ref + np.asarray(x) * scale


# ----------------------------------------------------------------------------
# 1.  PCA plane of the trajectory  (Li et al., "Visualizing the loss landscape", 2018)
# ----------------------------------------------------------------------------
def split_phases(trajs, frac=0.05):
    """Iteration at which every trajectory has completed (1-frac) of its total loss decrease.
    Phase I (fast, n·d-driven) = before it; Phase II (slow, κ-driven) = after."""
    b = 0
    for t in trajs:
        L = t['loss']; target = L.min() + frac * (L[0] - L.min())
        b = max(b, int(np.argmax(L <= target)))
    return b


def curvature_metric(H):
    """Per-parameter weights w_i = sqrt(|H_ii|) / geometric-mean, so that a unit step along any
    axis raises the loss by the same amount.  Use as `metric=` in pca_plane: PCA is then done in
    the loss's own units instead of 'fraction of nominal', and directions the loss ignores
    (here κ) stop dominating the plane."""
    w = np.sqrt(np.abs(np.diag(H)))
    return w / np.exp(np.mean(np.log(w)))


def pca_plane(trajs, ref, scale, which='all', weight_late=True, iter_range=None, metric=None):
    """Two orthonormal directions that capture most of the path variance.

    Coordinates: x = (p - ref)/scale, optionally multiplied by `metric` (see curvature_metric).
    weight_late: Adam takes tiny steps near the minimum; weighting points equally
    would let the (few) early large-step points dominate.  We log-sample the
    iterations so each decade of the descent has similar weight.
    """
    w = np.ones(len(ref)) if metric is None else np.asarray(metric, float)
    pts = []
    for t in (trajs if which == 'all' else [trajs[which]]):
        X = to_scaled(t['params'], ref, scale) * w
        if iter_range is not None:
            X = X[iter_range[0]:iter_range[1]]
        if weight_late and len(X) > 400:
            idx = np.unique(np.logspace(0, np.log10(len(X) - 1), 400).astype(int))
            X = X[idx]
        pts.append(X)
    X = np.vstack(pts)
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S**2 / np.sum(S**2)
    return {'u': Vt[0], 'v': Vt[1], 'explained': var[:2], 'ref': ref, 'scale': scale,
            'all_var': var, 'components': Vt, 'iter_range': iter_range, 'w': w}


def plane_point(plane, a, b):
    """Parameter vector at plane coordinates (a, b)."""
    return from_scaled((a * plane['u'] + b * plane['v']) / plane['w'], plane['ref'], plane['scale'])


def project(plane, params):
    X = to_scaled(params, plane['ref'], plane['scale']) * plane['w']
    a = X @ plane['u']; b = X @ plane['v']
    resid = np.linalg.norm(X - np.outer(a, plane['u']) - np.outer(b, plane['v']), axis=1)
    return a, b, resid


def _sl(plane):
    r = plane.get('iter_range'); return slice(None) if r is None else slice(r[0], r[1])


def plane_extent(plane, trajs, margin=0.25):
    A, B = [], []
    for t in trajs:
        a, b, _ = project(plane, t['params'][_sl(plane)]); A.append(a); B.append(b)
    A = np.concatenate(A + [np.zeros(1)]); B = np.concatenate(B + [np.zeros(1)])
    da = np.ptp(A) * margin + 1e-9; db = np.ptp(B) * margin + 1e-9
    return (A.min() - da, A.max() + da, B.min() - db, B.max() + db)


def evaluate_plane(loss_fn, plane, extent=None, trajs=None, res=61, verbose=True):
    """Evaluate loss_fn(params) on the plane ref + a*u + b*v.  Returns dict for plotting."""
    if extent is None:
        extent = plane_extent(plane, trajs)
    a = np.linspace(extent[0], extent[1], res); b = np.linspace(extent[2], extent[3], res)
    Z = np.empty((res, res))
    for j, bj in enumerate(b):
        for i, ai in enumerate(a):
            Z[j, i] = loss_fn(plane_point(plane, ai, bj))
        if verbose and (j + 1) % 10 == 0:
            print(f'  plane row {j+1}/{res}', end='\r')
    if verbose: print()
    return {'a': a, 'b': b, 'Z': Z, 'plane': plane, 'kind': 'pca'}


def evaluate_pair_slice(loss_fn, ref, scale, i, j, half_width=(0.25, 0.25), res=61):
    """Axis-aligned 2-D slice through `ref` varying parameters i and j (scaled units).
    half_width may be a scalar or a (w_i, w_j) pair."""
    wi, wj = (half_width, half_width) if np.isscalar(half_width) else half_width
    a = np.linspace(-wi, wi, res); b = np.linspace(-wj, wj, res)
    Z = np.empty((res, res))
    for jj, bj in enumerate(b):
        for ii, ai in enumerate(a):
            x = np.zeros(len(ref)); x[i] = ai; x[j] = bj
            Z[jj, ii] = loss_fn(from_scaled(x, ref, scale))
    return {'a': a, 'b': b, 'Z': Z, 'ref': ref, 'scale': scale, 'i': i, 'j': j, 'kind': 'pair'}


# ----------------------------------------------------------------------------
# drawing primitives
# ----------------------------------------------------------------------------
def _landscape(ax, a, b, Z, levels=14, cbar=True, cbar_label='loss', vmin=None, vmax=None):
    Zp = np.clip(Z, np.nanmin(Z[Z > 0]) if np.any(Z > 0) else 1e-12, None)
    vmin = vmin or Zp.min(); vmax = vmax or Zp.max()
    lv = np.logspace(np.log10(vmin), np.log10(vmax), levels)
    cf = ax.contourf(a, b, Zp, levels=lv, norm=LogNorm(vmin, vmax), cmap=LANDSCAPE_CMAP, alpha=0.9, extend='both')
    ax.contour(a, b, Zp, levels=lv, colors='white', linewidths=0.4, alpha=0.7)
    cf.levels_used = lv
    if cbar:
        _colorbar(cf, ax, cbar_label)
    return cf


def _colorbar(cf, ax, label='loss', **kw):
    lv = cf.levels_used; step = max(1, len(lv) // 6)
    cb = plt.colorbar(cf, ax=ax, pad=kw.pop('pad', 0.02), fraction=kw.pop('fraction', 0.046), ticks=lv[::step], **kw)
    cb.ax.set_yticklabels([f'{v:.3g}' for v in lv[::step]])
    cb.set_label(label); cb.outline.set_visible(False)
    return cb


def _path(ax, a, b, color, label=None, lw=1.6, milestones=True, fade=True, zorder=5, offset=0):
    """Draw a trajectory as a line whose opacity increases with iteration + milestone dots."""
    if fade:
        pts = np.column_stack([a, b]).reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
        alphas = np.linspace(0.35, 1.0, len(seg))
        cols = np.array([plt.matplotlib.colors.to_rgba(color, al) for al in alphas])
        lc = LineCollection(seg, colors=cols, linewidths=lw, zorder=zorder, capstyle='round')
        ax.add_collection(lc)
    else:
        ax.plot(a, b, color=color, lw=lw, zorder=zorder)
    ax.plot(a[0], b[0], 'o', color=color, ms=6, mec='white', mew=0.8, zorder=zorder + 1, label=label)
    ax.plot(a[-1], b[-1], '*', color=color, ms=11, mec='white', mew=0.6, zorder=zorder + 2)
    if milestones:
        for m in MILESTONES:
            m -= offset
            if 0 < m < len(a) - 1:
                ax.plot(a[m], b[m], 'o', color='white', ms=4, mec=color, mew=1.0, zorder=zorder + 1)
                ax.annotate(f'{m + offset:,}', (a[m], b[m]), xytext=(4, 4), textcoords='offset points',
                            fontsize=6.5, color=INK2, zorder=zorder + 3)


def _fmt_dir(vec, names=PARAM_NAMES, k=3):
    vec = np.asarray(vec) / np.linalg.norm(vec)
    idx = np.argsort(-np.abs(vec))[:k]
    return ' '.join(f'{vec[i]:+.2f}{PARAM_TEX[names[i]]}' for i in idx)


# ----------------------------------------------------------------------------
# 2.  main figures
# ----------------------------------------------------------------------------
def plot_pca_landscape(grid, trajs, ax=None, title=None, show_dirs=True, cbar=True):
    """Log-colour contour map on the PCA plane with all trajectories projected onto it."""
    plane = grid['plane']
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
    _landscape(ax, grid['a'], grid['b'], grid['Z'], cbar=cbar)
    r = plane.get('iter_range')
    for t, c in zip(trajs, TRAJ_COLORS):
        a, b, _ = project(plane, t['params'])
        if r is not None:                      # remainder of the path, faint
            ax.plot(a, b, color=c, lw=0.8, alpha=0.3, zorder=4)
            a, b = a[r[0]:r[1]], b[r[0]:r[1]]
        _path(ax, a, b, c, label=t['label'], offset=(r[0] if r else 0))
    ax.plot(0, 0, '+', color=INK, ms=10, mew=1.2, zorder=20)
    ev = plane['explained'] * 100
    unit = ' [curvature-normalised units]' if not np.allclose(plane['w'], 1) else ' [fraction of nominal]'
    ax.set_xlabel(f'PC1  ({ev[0]:.1f} % of path variance){unit}' + (f'\n{_fmt_dir(plane["u"])}' if show_dirs else ''))
    ax.set_ylabel(f'PC2  ({ev[1]:.1f} %)' + (f'\n{_fmt_dir(plane["v"])}' if show_dirs else ''))
    ax.set_title(title or 'Loss landscape on the plane of the optimiser path')
    ax.legend(loc='upper right', fontsize=8, title='● start  ★ end  ○ iteration', title_fontsize=7)
    ax.set_xlim(grid['a'][0], grid['a'][-1]); ax.set_ylim(grid['b'][0], grid['b'][-1])
    _despine(ax); ax.grid(False)
    return ax


def plot_two_phase_landscape(loss_fn, trajs, ref, scale, frac=0.05, res=61, figsize=(13, 9),
                             metric=None, verbose=True):
    """Two-phase picture of the descent (see split_phases).

    Row 1 - Phase I, one panel per start: each path in *its own* PCA plane (curvature-normalised
            if `metric` is given), because the early descents from different starts do not share a
            plane; all three panels share one log colour scale.
    Row 2 - Phase II in one shared plane (plain scaled units, so the κ valley is visible), plus the
            convergence curves with the phase boundary marked.
    Returns fig, dict(grids1=[...], grid2=..., boundary=b).
    """
    b = split_phases(trajs, frac)
    N = max(len(t['loss']) for t in trajs)
    if verbose: print(f'phase boundary at iteration {b}')
    grids1 = []
    for k, t in enumerate(trajs):
        p = pca_plane(trajs, ref, scale, which=k, iter_range=(0, b + 1), weight_late=False, metric=metric)
        grids1.append(evaluate_plane(loss_fn, p, trajs=[t], res=res, verbose=verbose))
    p2 = pca_plane(trajs, ref, scale, iter_range=(b, N))
    g2 = evaluate_plane(loss_fn, p2, trajs=trajs, res=res, verbose=verbose)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1, 1.1], hspace=0.42, wspace=0.35)
    vmin = min(g['Z'].min() for g in grids1); vmax = max(g['Z'].max() for g in grids1)
    for k, (g, t) in enumerate(zip(grids1, trajs)):
        ax = fig.add_subplot(gs[0, k])
        cf = _landscape(ax, g['a'], g['b'], g['Z'], cbar=False, vmin=vmin, vmax=vmax)
        a, bb, _ = project(g['plane'], t['params'])
        ax.plot(a, bb, color=TRAJ_COLORS[k], lw=0.8, alpha=0.3, zorder=4)         # phase II remainder, faint
        _path(ax, a[:b + 1], bb[:b + 1], TRAJ_COLORS[k], label=t['label'])
        ax.plot(0, 0, '+', color=INK, ms=10, mew=1.2, zorder=20)
        ev = g['plane']['explained'] * 100
        ax.set_xlabel(f'PC1 ({ev[0]:.0f} %)  {_fmt_dir(g["plane"]["u"])}')
        ax.set_ylabel(f'PC2 ({ev[1]:.0f} %)  {_fmt_dir(g["plane"]["v"])}')
        ax.set_title(f'({"abc"[k]}) Phase I, {t["label"]}: iterations 1–{b:,}')
        ax.set_xlim(g['a'][0], g['a'][-1]); ax.set_ylim(g['b'][0], g['b'][-1])
        ax.legend(loc='best', fontsize=7.5, title='● start  ★ end of phase  ○ iteration', title_fontsize=6.5)
        _despine(ax); ax.grid(False)
    _colorbar(cf, [fig.axes[i] for i in range(3)], 'loss', fraction=0.015, pad=0.01)
    unit = 'curvature-normalised units' if metric is not None else 'fraction of nominal'
    fig.text(0.5, 0.505, f'Phase I planes: per-start PCA of the path, coordinates in {unit}; + marks the global optimum',
             ha='center', fontsize=8, color=INK2)

    ax = fig.add_subplot(gs[1, :2])
    plot_pca_landscape(g2, trajs, ax=ax, title=f'(d) Phase II, shared plane: iterations {b:,}–{N:,}  (slow drift along the κ valley)')
    ax.set_xlabel(ax.get_xlabel().replace(' [fraction of nominal]', ''))
    ax = fig.add_subplot(gs[1, 2])
    for t, c in zip(trajs, TRAJ_COLORS):
        ax.plot(np.arange(1, len(t['loss']) + 1), t['loss'], color=c, lw=1.4, label=t['label'])
    ax.axvline(b, color=INK2, lw=0.8, ls='--'); ax.text(b * 1.1, ax.get_ylim()[1] * 0.9, 'Phase I | II', fontsize=7.5, color=INK2)
    ax.set(xscale='log', yscale='log', xlabel='iteration', title='(e) Convergence: loss vs iteration'); ax.legend(fontsize=8); _despine(ax)
    return fig, {'grids1': grids1, 'grid2': g2, 'boundary': b}


def plot_pca_diagnostics(grid, trajs, loss_fn=None, figsize=(11, 3.4)):
    """Three panels: (a) loss vs iteration, (b) loss along path vs in-plane projection
    (how faithful the plane is), (c) out-of-plane distance vs iteration."""
    plane = grid['plane']
    fig, axs = plt.subplots(1, 3, figsize=figsize)
    for t, c in zip(trajs, TRAJ_COLORS):
        it = np.arange(1, len(t['loss']) + 1)
        axs[0].plot(it, t['loss'], color=c, lw=1.4, label=t['label'])
        a, b, resid = project(plane, t['params'])
        if loss_fn is not None:
            idx = np.unique(np.logspace(0, np.log10(len(a) - 1), 60).astype(int))
            Lp = [loss_fn(plane_point(plane, a[i], b[i])) for i in idx]
            axs[1].plot(it[idx], t['loss'][idx], color=c, lw=1.4)
            axs[1].plot(it[idx], Lp, color=c, lw=1.0, ls='--')
        axs[2].plot(it, resid, color=c, lw=1.4)
    axs[0].set(xscale='log', yscale='log', xlabel='iteration', ylabel='loss', title='(a) Convergence')
    axs[0].legend(fontsize=8)
    axs[1].set(xscale='log', yscale='log', xlabel='iteration', ylabel='loss',
               title='(b) True loss (—) vs loss of in-plane projection (- -)')
    axs[2].set(xscale='log', yscale='log', xlabel='iteration', ylabel='out-of-plane distance (scaled units)',
               title='(c) Distance of path from the PCA plane')
    for ax in axs: _despine(ax)
    fig.tight_layout()
    return fig


def plot_pair_slices(loss_fn, ref, scale, trajs, pairs=None, half_widths=None, res=41, figsize=(11.5, 6.6),
                     show_optical_thickness=True, verbose=True):
    """Classic 2x3 grid of axis-aligned slices *through the optimum*, improved:
    shared log colour scale, per-parameter axis ranges (κ needs a far wider range than n, d),
    trajectories as lines, and the n·d = const (optical-thickness) hyperbola in the n-d panels.
    Note: trajectory points are *projected* onto each slice - their true loss is not the slice value."""
    if pairs is None:
        pairs = [(0, 2), (0, 1), (2, 1), (3, 5), (3, 4), (5, 4)]
    if half_widths is None:   # scaled units (fraction of nominal value)
        half_widths = {'n': 0.15, 'k': 1.5, 'd': 0.15}
    hw = lambda idx: half_widths[PARAM_NAMES[idx][0]]
    grids = []
    for (i, j) in pairs:
        if verbose: print(f'slice {PARAM_NAMES[i]}-{PARAM_NAMES[j]}', end='  ', flush=True)
        grids.append(evaluate_pair_slice(loss_fn, ref, scale, i, j, (hw(i), hw(j)), res))
    if verbose: print()
    vmin = min(g['Z'].min() for g in grids); vmax = max(g['Z'].max() for g in grids)
    fig, axs = plt.subplots(2, 3, figsize=figsize, layout='constrained')
    for ax, g in zip(axs.ravel(), grids):
        i, j = g['i'], g['j']
        cf = _landscape(ax, 100 * g['a'], 100 * g['b'], g['Z'], cbar=False, vmin=vmin, vmax=vmax)
        for t, c in zip(trajs, TRAJ_COLORS):
            X = to_scaled(t['params'], ref, scale) * 100
            _path(ax, X[:, i], X[:, j], c, label=t['label'], milestones=False, lw=1.3)
        ax.plot(0, 0, '+', color=INK, ms=9, mew=1.2, zorder=20)
        if show_optical_thickness and {PARAM_NAMES[i][0], PARAM_NAMES[j][0]} == {'n', 'd'}:
            # constant optical thickness n*d through the optimum: (1+xi)(1+xj)=1
            xi = np.linspace(-hw(i), hw(i), 200)
            ax.plot(100 * xi, 100 * (1 / (1 + xi) - 1), ls=':', color='#e34948', lw=1.2, zorder=6)
            ax.text(0.03, 0.05, r'$n\,d$ = const', transform=ax.transAxes, color='#e34948', fontsize=7.5)
        ax.set_xlabel(f'{PARAM_TEX[PARAM_NAMES[i]]} deviation from optimum (% of nominal)')
        ax.set_ylabel(f'{PARAM_TEX[PARAM_NAMES[j]]} deviation (% of nominal)')
        ax.set_xlim(-100 * hw(i), 100 * hw(i)); ax.set_ylim(-100 * hw(j), 100 * hw(j))
        _despine(ax); ax.grid(False)
    axs[0, 1].legend(loc='upper left', fontsize=7.5)
    _colorbar(cf, axs, 'loss', fraction=0.025, pad=0.02)
    fig.suptitle('Axis-aligned slices of the loss through the optimum (shared log colour; paths are projections)')
    return fig


# ----------------------------------------------------------------------------
# 3.  Local geometry at the minimum: sensitivity (Gauss-Newton) matrix
# ----------------------------------------------------------------------------
# Why not the Hessian?  The loss here is MAE + log-cosh, i.e. V-shaped ("cone") around the
# optimum, not a bowl.  Its point-wise Hessian is ill-defined and a finite-difference Hessian
# of a cone even reports zero coupling between parameters that are in fact degenerate.
# The natural geometry is the sensitivity matrix  G = J^T J / L  with  J = d y_sim / d x
# (x = scaled parameters).  Its eigenvectors are the directions in which the simulated pulse
# changes least/most; small eigenvalues = poorly constrained combinations (e.g. n·d).
def sensitivity_matrix(ysim_torch_fn, ref, scale):
    """G = J^T J / L in scaled coordinates at `ref`.  ysim_torch_fn(p_tensor[P]) -> tensor [L]."""
    import torch
    ref_t = torch.tensor(ref, dtype=torch.float32); sc_t = torch.tensor(scale, dtype=torch.float32)
    f = lambda x: ysim_torch_fn(ref_t + x * sc_t)
    J = torch.autograd.functional.jacobian(f, torch.zeros(len(ref))).detach().numpy().astype(float)  # [L, P]
    return J.T @ J / J.shape[0], J


def curvature_metric(G):
    """Per-parameter weights w_i = sqrt(G_ii) / geometric-mean: a unit step along any axis then
    changes the simulated pulse by the same RMS amount.  Use as `metric=` in pca_plane so that
    directions the data cannot see (κ) stop dominating the plane."""
    w = np.sqrt(np.abs(np.diag(G)))
    return w / np.exp(np.mean(np.log(w)))


def hessian_fd(loss_fn, ref, scale, eps=0.02):
    """Finite-difference 'Hessian' at probe scale eps (fraction of nominal).  Kept for reference;
    for a V-shaped loss prefer sensitivity_matrix (see note above)."""
    P = len(ref)
    f0 = loss_fn(ref); H = np.zeros((P, P))
    E = np.eye(P) * eps
    fp = np.array([loss_fn(from_scaled(E[i], ref, scale)) for i in range(P)])
    fm = np.array([loss_fn(from_scaled(-E[i], ref, scale)) for i in range(P)])
    g = (fp - fm) / (2 * eps)
    for i in range(P):
        H[i, i] = (fp[i] - 2 * f0 + fm[i]) / eps**2
        for j in range(i + 1, P):
            fpp = loss_fn(from_scaled(E[i] + E[j], ref, scale)); fmm = loss_fn(from_scaled(-E[i] - E[j], ref, scale))
            H[i, j] = H[j, i] = (fpp - fp[i] - fp[j] + 2 * f0 - fm[i] - fm[j] + fmm) / (2 * eps**2)
    return H, g


def plot_geometry(G, ref, scale, loss_fn=None, names=PARAM_NAMES, figsize=(12.5, 3.8), probe_width=0.15):
    """(a) singular-value spectrum of the sensitivity matrix, (b) eigenvector composition,
    (c) true loss along each eigen-direction.  Shows anisotropy (which combinations the data
    constrain) and the V shape of the loss."""
    w, V = np.linalg.eigh(G)
    order = np.argsort(-w); w, V = w[order], V[:, order]
    sv = np.sqrt(np.clip(w, 0, None))
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 3, width_ratios=[1, 1.15, 1.35], figure=fig, wspace=0.38)
    # (a)
    ax = fig.add_subplot(gs[0])
    ax.bar(range(len(sv)), sv, color='#2a78d6', width=0.6)
    ax.set_yscale('log'); ax.set_xticks(range(len(sv))); ax.set_xticklabels([f'σ{i+1}' for i in range(len(sv))])
    ax.set_ylabel('RMS pulse change per unit step  (σ = √eig JᵀJ/L)')
    ax.set_title(f'(a) Sensitivity spectrum   σ₁/σ₆ = {sv[0]/max(sv[-1],1e-300):.0f}')
    _despine(ax)
    # (b)
    ax = fig.add_subplot(gs[1])
    im = ax.imshow(V, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(len(sv))); ax.set_xticklabels([f'σ{i+1}' for i in range(len(sv))])
    ax.set_yticks(range(len(names))); ax.set_yticklabels([PARAM_TEX[n] for n in names])
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            if abs(V[i, j]) > 0.25:
                ax.text(j, i, f'{V[i, j]:+.2f}', ha='center', va='center', fontsize=6.5,
                        color='white' if abs(V[i, j]) > 0.6 else INK)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.outline.set_visible(False); cb.set_label('component')
    ax.set_title('(b) Eigen-directions (well → poorly constrained)')
    # (c)
    ax = fig.add_subplot(gs[2])
    if loss_fn is not None:
        s_ = np.linspace(-probe_width, probe_width, 41)
        for k in range(len(sv)):
            c = plt.cm.viridis(k / max(len(sv) - 1, 1))
            Ltrue = [loss_fn(from_scaled(si * V[:, k], ref, scale)) for si in s_]
            ax.plot(100 * s_, Ltrue, color=c, lw=1.4, label=f'σ{k+1}')
        ax.set_yscale('log'); ax.set_xlabel('displacement along eigen-direction (% of nominal)')
        ax.set_ylabel('loss'); ax.legend(fontsize=7, ncol=2)
        ax.set_title('(c) Loss along eigen-directions (note the V shape)')
    _despine(ax)
    return fig, (sv, V)


plot_curvature = None   # removed: see plot_geometry


# ----------------------------------------------------------------------------
# 4.  1-D interpolation between start and end (barrier check)
# ----------------------------------------------------------------------------
def plot_line_interpolation(loss_fn, trajs, n=61, extend=0.25, ax=None):
    """Loss along the straight line from each start to its end point (alpha in [-extend, 1+extend]).
    A bump between 0 and 1 means the optimiser had to go around a barrier."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.2, 3.4))
    al = np.linspace(-extend, 1 + extend, n)
    for t, c in zip(trajs, TRAJ_COLORS):
        p0, p1 = t['params'][0], t['params'][int(np.argmin(t['loss']))]
        Ls = [loss_fn(p0 + a * (p1 - p0)) for a in al]
        ax.plot(al, Ls, color=c, lw=1.5, label=t['label'])
    ax.axvline(0, color=GRID, lw=1); ax.axvline(1, color=GRID, lw=1)
    ax.set_yscale('log'); ax.set_xlabel(r'$\alpha$   (0 = start, 1 = optimum)'); ax.set_ylabel('loss')
    ax.set_title('Loss along the straight line start → optimum'); ax.legend(fontsize=8)
    _despine(ax)
    return ax


# ----------------------------------------------------------------------------
# 5.  Animation
# ----------------------------------------------------------------------------
def animate_descent(grid, trajs, filename='descent.gif', n_frames=80, fps=12, figsize=(10, 4.6), dpi=110):
    """GIF of the paths growing over the PCA landscape (log-spaced in iteration) + loss curve."""
    plane = grid['plane']
    N = max(len(t['loss']) for t in trajs)
    frames = np.unique(np.concatenate([[1, 2, 3, 5, 7], np.logspace(1, np.log10(N), n_frames).astype(int)]))
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, width_ratios=[1.25, 1], figure=fig, wspace=0.3)
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
    _landscape(ax, grid['a'], grid['b'], grid['Z'], cbar=False)
    ax.plot(0, 0, '+', color=INK, ms=10, mew=1.2, zorder=20)
    ev = plane['explained'] * 100
    ax.set_xlabel(f'PC1 ({ev[0]:.0f} %)   {_fmt_dir(plane["u"])}'); ax.set_ylabel(f'PC2 ({ev[1]:.0f} %)   {_fmt_dir(plane["v"])}')
    ax.set_xlim(grid['a'][0], grid['a'][-1]); ax.set_ylim(grid['b'][0], grid['b'][-1])
    ax.grid(False)
    proj = [project(plane, t['params'])[:2] for t in trajs]
    lines, heads, lcurves, ldots = [], [], [], []
    for (a, b), t, c in zip(proj, trajs, TRAJ_COLORS):
        l, = ax.plot([], [], color=c, lw=1.6, zorder=5); h, = ax.plot([], [], 'o', color=c, ms=7, mec='white', zorder=6)
        ax.plot(a[0], b[0], 'o', color=c, ms=5, mfc='none', zorder=6)
        lc, = ax2.plot([], [], color=c, lw=1.4, label=t['label']); ld, = ax2.plot([], [], 'o', color=c, ms=6, mec='white')
        lines.append(l); heads.append(h); lcurves.append(lc); ldots.append(ld)
    ax2.set(xscale='log', yscale='log', xlabel='iteration', ylabel='loss', xlim=(1, N * 1.1),
            ylim=(min(t['loss'].min() for t in trajs) * 0.8, max(t['loss'].max() for t in trajs) * 1.3))
    ax2.legend(fontsize=8); _despine(ax2)
    title = ax.set_title('')

    def update(k):
        it = frames[k]
        for (a, b), t, l, h, lc, ld in zip(proj, trajs, lines, heads, lcurves, ldots):
            m = min(it, len(a))
            l.set_data(a[:m], b[:m]); h.set_data([a[m - 1]], [b[m - 1]])
            lc.set_data(np.arange(1, m + 1), t['loss'][:m]); ld.set_data([m], [t['loss'][m - 1]])
        title.set_text(f'Adam descent on the PCA plane — iteration {it:,}')
        return lines + heads + lcurves + ldots + [title]

    ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000 / fps, blit=False)
    ani.save(filename, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return filename


