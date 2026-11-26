"""
paperstyle.py — shared matplotlib style for every figure of the paper.

Single-column Optica page: text width ≈ 13.3 cm.  All figures are drawn at their final
printed size so that no font is smaller than 7 pt after inclusion with \\linewidth.

Colour-blind-safe categorical palette (validated with the dataviz palette checker:
lightness band, chroma floor, CVD separation ≥ 8, normal-vision separation ≥ 15):
    C[0] blue   #2a78d6   C[1] orange #eb6834   C[2] aqua  #1baf7a
    C[3] yellow #eda100   C[4] magenta #e87ba4  C[5] violet #6250d6   C[6] red #e34948
Series are assigned in fixed order (never cycled); any pair that is hard to separate in
print is also separated by line style or marker.

Usage in a notebook
-------------------
    import sys; sys.path.insert(0, '../code')
    import paperstyle as ps
    ps.apply()                       # call AFTER importing loss_landscape_viz (which sets its own rcParams)
    fig, ax = plt.subplots(figsize=ps.figsize(1.0, 0.5))
    ps.panel_label(ax, 'a')
    ps.savefig(fig, '../figs/fig1_five_layer')   # writes .png (300 dpi) and .pdf
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

TEXT_WIDTH_CM = 13.3
TEXT_WIDTH_IN = TEXT_WIDTH_CM / 2.54

# categorical slots, fixed order
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#6250d6', '#e34948']
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, VIOLET, RED = C
INK = '#0b0b0b'       # primary text
INK2 = '#52514e'      # secondary text / axes
GRID = '#d9d8d3'      # hairline grid
NEUTRAL = '#9a9892'   # neutral grey for reference lines / bands
LANDSCAPE_CMAP = 'Greys'   # sequential single hue for loss maps (coloured paths stay legible)

# roles used across the notebooks (same entity -> same colour in every figure)
ROLE = {
    'reference': INK2,   # measured reference pulse: neutral
    'sample':    BLUE,   # measured sample pulse
    'fit':       ORANGE, # fitted / simulated pulse
    'init':      AQUA,   # pulse from the initial guess
    'residual':  VIOLET, # fit residual
    'start1': BLUE, 'start2': ORANGE, 'start3': AQUA,   # the three starts of the quartz/Si fit
    'adam':  BLUE, 'lbfgs': ORANGE,
    'n': BLUE, 'k': ORANGE, 'd': AQUA,                  # parameter families
    'quartz': BLUE, 'Si': ORANGE, 'air': AQUA,
    'stage1': NEUTRAL, 'stage2': BLUE,
}

RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'legend.title_fontsize': 7,
    'legend.frameon': False,
    'legend.handlelength': 1.6,
    'legend.borderaxespad': 0.3,
    'axes.edgecolor': INK2,
    'axes.labelcolor': INK,
    'axes.linewidth': 0.6,
    'axes.titleweight': 'normal',
    'axes.titlepad': 3.0,
    'axes.labelpad': 2.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': False,
    'xtick.color': INK2, 'ytick.color': INK2,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
    'xtick.major.pad': 2, 'ytick.major.pad': 2,
    'lines.linewidth': 0.9,
    'lines.markersize': 3.5,
    'grid.color': GRID, 'grid.linewidth': 0.4,
    'figure.dpi': 110,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'pdf.fonttype': 42, 'ps.fonttype': 42,   # editable text in PDF
    'mathtext.fontset': 'dejavusans',
}


def apply():
    """Install the paper rcParams (call after any module that sets its own rcParams)."""
    mpl.rcParams.update(RC)


def figsize(width=1.0, aspect=0.55, height_cm=None):
    """(w, h) in inches: width as a fraction of the text width; height from aspect (h/w) or absolute cm."""
    w = width * TEXT_WIDTH_IN
    h = (height_cm / 2.54) if height_cm is not None else aspect * w
    return (w, h)


def panel_label(ax, letter, x=-0.18, y=1.02, **kw):
    """Bold (a), (b), ... in the top-left corner outside the axes."""
    ax.text(x, y, f'({letter})', transform=ax.transAxes, fontsize=8, fontweight='bold',
            ha='left', va='bottom', color=INK, **kw)


def despine(ax, grid=False):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    if grid:
        ax.grid(True, color=GRID, lw=0.4)
        ax.set_axisbelow(True)


FONT_LOG = {}          # filled by savefig: per figure, exported width and smallest printed font size


def font_report(fig, target_width_cm=TEXT_WIDTH_CM):
    """Exported width (tight bbox + padding, as written by savefig) and the smallest font of any visible,
    non-empty text.  If the file is placed at `target_width_cm` (\\linewidth = 13.3 cm, or 6.6 cm for a
    half-width figure), every font prints at font_pt * target_width_cm / width_cm."""
    from matplotlib.text import Text
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = fig.get_tightbbox(r)
    pad = mpl.rcParams['savefig.pad_inches']
    w_cm = (bb.width + 2 * pad) * 2.54; h_cm = (bb.height + 2 * pad) * 2.54
    sizes = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    fmin = float(min(sizes)) if sizes else float('nan')
    scale = target_width_cm / w_cm
    return dict(width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), target_width_cm=target_width_cm,
                min_font_pt=fmin, scale_at_target=round(scale, 4), min_printed_font_pt=round(fmin * scale, 3),
                ok=bool(fmin * scale >= 7.0 - 1e-6))


def savefig(fig, path_no_ext, png=True, pdf=True, dpi=300, target_width_cm=TEXT_WIDTH_CM):
    """Save PNG (300 dpi) and PDF; `path_no_ext` without extension.  Also checks the printed font size at the
    final width (see font_report) and records it in FONT_LOG."""
    rep = font_report(fig, target_width_cm)
    for _ in range(4):                       # labels sticking out of the canvas make the tight bbox wider than the
        if rep['width_cm'] <= target_width_cm + 1e-3 and rep['ok']:   # target: shrink until the width fits and fonts print >= 7 pt
            break
        w, h = fig.get_size_inches()
        fig.set_size_inches(w * (target_width_cm - 0.02) / rep['width_cm'], h)
        rep = font_report(fig, target_width_cm)
    FONT_LOG[path_no_ext.split('/')[-1]] = rep
    if png:
        fig.savefig(path_no_ext + '.png', dpi=dpi)
    if pdf:
        fig.savefig(path_no_ext + '.pdf')
    print(f'saved {path_no_ext}.png' + (' + .pdf' if pdf else '') +
          f"  [{rep['width_cm']:.2f} x {rep['height_cm']:.2f} cm; placed at {target_width_cm} cm the smallest font prints at "
          f"{rep['min_printed_font_pt']:.2f} pt{'' if rep['ok'] else '  <-- BELOW 7 pt'}]")


def to_jsonable(obj):
    """Recursively convert numpy / torch scalars and arrays to plain Python for json.dump."""
    import numbers
    try:
        import torch
        if isinstance(obj, torch.Tensor):
            obj = obj.detach().cpu().numpy()
    except ImportError:
        pass
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [to_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, complex):
        return {'re': obj.real, 'im': obj.imag}
    if isinstance(obj, numbers.Number) or obj is None or isinstance(obj, (str, bool)):
        return obj
    return str(obj)
