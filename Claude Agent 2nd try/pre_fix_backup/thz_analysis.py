#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
thz_analysis.py
===============================================================================
THz time-domain spectroscopy (THz-TDS) analysis of a two-layer dielectric stack

        air  |  silicon (d_Si)  |  fused quartz (d_qz)  |  air

measured on a Menlo Systems ScanControl THz-TDS spectrometer under N2 purge.

===============================================================================
PHYSICS / SIGN CONVENTIONS  (read this before touching the code)
===============================================================================
numpy.fft.rfft evaluates   X(f) = sum_n x[n] exp(-2*pi*i*f*t_n),
so a pure time delay tau maps to a factor exp(-i*w*tau).  Consistency therefore
REQUIRES the complex refractive index to be written

        n~ = n - i*k        with  k >= 0                                    (1)

so that propagation across a slab of thickness d contributes

        P(w) = exp(-i*n~*w*d/c) = exp(-i*n*w*d/c) * exp(-k*w*d/c)           (2)
               \______ phase delay ______/  \___ amplitude decay ___/

i.e. a decaying wave.  (The opposite convention n~ = n + i*k belongs with an
exp(+i*w*t) carrier and the opposite FFT sign; mixing the two is the classic
source of "negative absorption" in THz-TDS code.)

Power absorption coefficient:   alpha(w) = 2*w*k/c      [1/m]               (3)

Normal-incidence Fresnel coefficients, medium i -> medium j:

        r_ij = (n~_i - n~_j)/(n~_i + n~_j),   t_ij = 2*n~_i/(n~_i + n~_j)   (4)

===============================================================================
THE INVERSION PROBLEM AND HOW THIS SCRIPT SOLVES IT
===============================================================================
One complex transmission spectrum T(w) = E_sam(w)/E_ref(w) supplies TWO real
numbers per frequency, but the bilayer has FOUR unknowns per frequency
(n_Si, k_Si, n_qz, k_qz) plus two thicknesses.  Pointwise inversion of T(w)
alone is hopelessly under-determined.  The Fabry-Perot echo train supplies the
missing information.

Label the media 0 = air(in), 1 = Si, 2 = quartz, 3 = air(out).  The transmitted
field is a sum of partial waves:

  main pulse   A_0 = t01 t12 t23 * P1 * P2
  Si echo      A_1 = A_0 * (r12 * r10) * P1^2      <- one extra round trip in Si
  quartz echo  A_2 = A_0 * (r23 * r21) * P2^2      <- one extra round trip in qz

Ratios of time-gated sub-pulses cancel every transmission coefficient and the
whole "other" layer exactly:

  R_Si(w) = A_1/A_0 = r12*r10 * exp(-2i*n~_Si*w*d_Si/c)                     (5)
  R_qz(w) = A_2/A_0 = r23*r21 * exp(-2i*n~_qz*w*d_qz/c)                     (6)

Eq. (5) involves silicon only, Eq. (6) quartz only.  That is what makes a
per-layer extraction possible at all.

--- Gating the DECONVOLVED impulse response, not the raw waveform -------------
The raw reference waveform of this instrument carries system ringing (emitter/
detector and purge-box etalons) at +7.2, +7.9 and +8.9 ps after the main pulse,
at up to 2.6% of the peak.  That lands directly on top of the Si echo gate
(+8.78 ps), which is only 14% of the main pulse -> a ~16% amplitude error.
Those artifacts are common to sample and reference, so they cancel exactly in
T(w).  This script therefore band-limits T(w) and inverse-transforms it to the
SAMPLE IMPULSE RESPONSE h(t), and gates the echoes there.  h(t) is clean: the
2nd-order Si echo then appears at exactly 2x the 1st-order delay, and the
Si+quartz cross echo at exactly the sum of the two - see the printed output.

--- Separating n from d (the classic THz-TDS degeneracy) ---------------------
The PHASE of (5)/(6) determines only the optical thickness n*d.  Three measured
delays give three exact, model-free relations:

  c*dt_main = (n1-1)d1 + (n2-1)d2,   c*tau_Si/2 = n1*d1,   c*tau_qz/2 = n2*d2

Eliminating the optical thicknesses gives a thickness sum rule that uses no
material data whatsoever:

  d1 + d2 = c*( tau_Si/2 + tau_qz/2 - dt_main )                             (7)

That is 3 equations for 4 unknowns; exactly ONE external number is needed to
close the system.  This script anchors on n_Si, because high-resistivity float-
zone silicon is the best-characterised THz reference material in existence
(n = 3.4175, flat to <1e-3 from 0.1-4 THz; Grischkowsky et al., JOSA B 7, 2006
(1990)).  Everything else - d_Si, d_qz and n_qz - then follows from measured
delays alone, and d_Si is returned as an independent check against the quoted
wafer thickness.

A fully self-contained alternative exists: |R_Si| and |R_qz| are pure Fresnel
products, independent of thickness, so they can in principle be solved for
(n_Si, n_qz) directly.  That route is computed and printed too, but it is badly
ill-conditioned (the useful combination is |R_Si|/|R_qz| = F(n1)/F(n2) with
F(n) = (n-1)/(n+1), whose sensitivity to n is ~0.6, so it DE-amplifies: a 10%
amplitude error becomes a ~16% index error).  It is reported as a cross-check
only, never used for the quoted numbers.

--- Absorption ---------------------------------------------------------------
k_qz(w) is obtained by inverting the FULL Fabry-Perot transfer matrix against
the measured |T(w)| at each frequency (|T| decreases monotonically in k_qz, so
a vectorised bisection is exact and unconditionally stable).  This is where the
transfer-matrix model does real work rather than merely being plotted.
k_Si(w) comes from |R_Si|; because high-resistivity Si is essentially
transparent, the true value sits below this method's systematic floor, so it is
honestly reported as an UPPER BOUND rather than as a measured value.

===============================================================================
A NOTE ON THE FFT WINDOW  (please read - the brief asked for a Hann window)
===============================================================================
A Hann window spanning the WHOLE record is actively harmful here, for two
independent reasons:

 1. The sample and reference pulses sit at different absolute times (-716.99 ps
    vs -723.29 ps).  A full-record Hann multiplies them by different envelope
    values (0.623 vs 0.377) -> |T| inflated by a factor 1.65 (+4.4 dB), which
    appears as a spurious NEGATIVE absorption coefficient at every frequency.
 2. Even after aligning the pulses, a full-record Hann weights the Fabry-Perot
    echoes differently from the main pulse, distorting the very fringe contrast
    this analysis is built on.

The correct apodisation is a raised-cosine (Hann) TAPER that is flat across the
whole signal region and rolls off only in the noise-only record ends - i.e. a
Tukey window whose tapers are literally Hann half-cosines, reducing exactly to
a textbook Hann at taper fraction 1.0.  That is the default.  Run with
--hann-full to reproduce the naive full-record Hann and watch the absorption go
negative.  For the time-GATED sub-pulse analysis a true Hann gate IS used,
because each sub-pulse is centred in its own gate and so is weighted identically.

===============================================================================
Usage:   python3 thz_analysis.py [--hann-full] [--out FILE]
Outputs: thz_results.png  +  a printed summary table
===============================================================================
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.signal import hilbert, find_peaks

import matplotlib
matplotlib.use("Agg")                       # headless-safe
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Physical constants and analysis settings
# ---------------------------------------------------------------------------
C0 = 299_792_458.0             # speed of light in vacuum [m/s]
PS, UM, THZ = 1e-12, 1e-6, 1e12

NFFT = 2 ** 14                 # zero-padded FFT length (as specified)
TAPER_FRAC = 0.15              # Hann taper fraction; 1.0 -> textbook Hann
GATE_HW_PS = 2.0               # Hann gate half-width for sub-pulses [ps]
IRF_BAND_THZ = (0.15, 3.5)     # band-limit applied before forming h(t)
IRF_EDGE_THZ = 0.15            # raised-cosine edge width of that band-limit
CAL_BAND_THZ = (0.25, 0.60)    # low-loss band used to fix n*d and the Fresnel
FIT_BAND_THZ = (0.20, 3.00)    # reporting / plotting band
TABLE_FREQS_THZ = (0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5)

# The single external anchor: high-resistivity float-zone silicon, the standard
# THz reference material (Grischkowsky et al., JOSA B 7, 2006 (1990)).
N_SI_REF = 3.4175

# Quoted thicknesses - used ONLY for reporting agreement, never in the inversion
D_SI_NOMINAL, D_QZ_NOMINAL = 385.0 * UM, 1000.0 * UM

# --- plotting palette (validated CVD-safe, light mode) ---------------------
CLR_REF, CLR_SAM = "#2a78d6", "#eb6834"     # blue / orange
CLR_SI, CLR_QZ = "#2a78d6", "#eb6834"       # Si -> blue, quartz -> orange
CLR_MODEL, CLR_AUX = "#4a3aa7", "#1baf7a"   # violet model, aqua accent
SURFACE, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
MUTED, GRID, AXIS = "#898781", "#e1e0d9", "#c3c2b7"


# ===========================================================================
# 1.  File I/O
# ===========================================================================
def load_menlo(path):
    """Read a Menlo Systems ScanControl ASCII waveform -> (t [ps], E [a.u.], meta)."""
    meta = {"path": path, "name": os.path.basename(path)}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if (m := re.search(r"Average over\s+(\d+)\s+waveforms", line)):
                meta["averages"] = int(m.group(1))
            if (m := re.search(r"Timestamp:\s*([0-9T:\-]+)", line)):
                meta["timestamp"] = m.group(1)

    data = np.loadtxt(path, comments="#")
    t, e = data[:, 0].astype(float), data[:, 1].astype(float)
    order = np.argsort(t)
    t, e = t[order], e[order]

    dt = np.diff(t)
    if np.ptp(dt) > 1e-6 * np.mean(dt):
        raise ValueError(f"{path}: time axis is not uniformly sampled.")
    meta["dt_ps"], meta["n_pts"] = float(np.mean(dt)), t.size
    return t, e, meta


def find_data_files(folder):
    """Locate the sample and reference waveforms in `folder`."""
    files = sorted(glob.glob(os.path.join(folder, "*.txt")))
    if not files:
        raise FileNotFoundError(f"No .txt data files found in {folder!r}")
    ref = [f for f in files if "ref" in os.path.basename(f).lower()]
    sam = [f for f in files if f not in ref]
    if not ref or not sam:
        raise FileNotFoundError(
            "Could not classify files into sample/reference: "
            f"{[os.path.basename(f) for f in files]}")
    return sam[0], ref[0]


# ===========================================================================
# 2.  Windowing, spectra, band-limiting
# ===========================================================================
def raised_cosine_window(n, taper=1.0):
    """
    Raised-cosine (Hann) tapered window == Tukey window.

    `taper` is the fraction of the record occupied by the two Hann half-cosine
    roll-offs combined:
        taper = 1.0  -> exactly the textbook Hann window 0.5*(1-cos(2*pi*i/(N-1)))
        taper < 1.0  -> unity across the middle, Hann roll-off at the ends only
    """
    if not 0.0 < taper <= 1.0:
        raise ValueError("taper must lie in (0, 1]")
    w = np.ones(n)
    ntap = int(np.floor(taper * (n - 1) / 2.0)) + 1
    ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(ntap) / max(ntap - 1, 1)))
    w[:ntap], w[n - ntap:] = ramp, ramp[::-1]
    return w


def spectrum(t_ps, e, nfft=NFFT, window=None):
    """
    Zero-padded one-sided FFT.  The phase is referenced to t_ps[0]; because
    sample and reference share a time axis, E_sam/E_ref then carries the correct
    physical propagation phase with no further bookkeeping.
    """
    dt = (t_ps[1] - t_ps[0]) * PS
    y = e * (window if window is not None else 1.0)
    return np.fft.rfftfreq(nfft, d=dt), np.fft.rfft(y, n=nfft)


def band_limit(f_thz, lo, hi, edge):
    """Raised-cosine band-pass mask; smooth edges avoid Gibbs ringing in h(t)."""
    b = np.zeros_like(f_thz)
    up = (f_thz > lo - edge) & (f_thz < lo)
    b[up] = 0.5 * (1.0 + np.cos(np.pi * (lo - f_thz[up]) / edge))
    b[(f_thz >= lo) & (f_thz <= hi)] = 1.0
    dn = (f_thz > hi) & (f_thz < hi + edge)
    b[dn] = 0.5 * (1.0 + np.cos(np.pi * (f_thz[dn] - hi) / edge))
    return b


# ===========================================================================
# 3.  Time-domain helpers
# ===========================================================================
def refine_peak(t, y, idx):
    """Sub-sample peak position by parabolic interpolation about index `idx`."""
    if idx <= 0 or idx >= y.size - 1:
        return float(t[idx])
    a, b, c = y[idx - 1], y[idx], y[idx + 1]
    den = a - 2.0 * b + c
    shift = 0.0 if den == 0 else 0.5 * (a - c) / den
    return float(t[idx] + shift * (t[1] - t[0]))


def locate_main_pulse(t, e):
    """Arrival time [ps] of the main pulse, from the analytic-signal envelope."""
    env = np.abs(hilbert(e))
    return refine_peak(t, env, int(np.argmax(env)))


def locate_feature(t, y, t0, window_ps):
    """Locate the strongest |y| feature at a delay inside `window_ps` after t0."""
    mask = (t > t0 + window_ps[0]) & (t < t0 + window_ps[1])
    if not mask.any():
        return np.nan, np.nan
    sub = np.where(mask)[0]
    idx = sub[int(np.argmax(np.abs(y[sub])))]
    return refine_peak(t, np.abs(y), idx) - t0, float(y[idx])


def hann_gate(t, t_centre, half_width):
    """
    A true Hann window of total width 2*half_width centred at t_centre, zero
    elsewhere.  Applied to the FULL-LENGTH record so absolute timing - and hence
    the propagation phase of the gated sub-pulse - is preserved exactly.
    """
    g = np.zeros_like(t)
    m = np.abs(t - t_centre) <= half_width
    g[m] = 0.5 * (1.0 + np.cos(np.pi * (t[m] - t_centre) / half_width))
    return g


# ===========================================================================
# 4.  Optics: Fresnel coefficients and the Fabry-Perot transfer matrix
# ===========================================================================
def r_fresnel(ni, nj):
    """Normal-incidence amplitude reflection coefficient, medium i -> j."""
    return (ni - nj) / (ni + nj)


def t_fresnel(ni, nj):
    """Normal-incidence amplitude transmission coefficient, medium i -> j."""
    return 2.0 * ni / (ni + nj)


def transfer_matrix_stack(freq, n_list, d_list):
    """
    Complex amplitude transmission of an arbitrary stack embedded in air, via
    the 2x2 characteristic-matrix (transfer-matrix) formalism:

        M = D_01 . P_1 . D_12 . P_2 . ... . D_(N,N+1),        t = 1 / M[0,0]

    interface matrices   D_ij = (1/t_ij) [[1, r_ij], [r_ij, 1]]
    propagation matrices P_j  = diag(exp(+i*beta_j), exp(-i*beta_j)),
                         beta_j = n~_j * w * d_j / c.

    This resums ALL Fabry-Perot internal reflections in every layer and all
    coupled-cavity cross terms to infinite order - nothing is truncated.
    """
    w = 2.0 * np.pi * np.asarray(freq, dtype=float)
    one = np.ones_like(w, dtype=complex)
    media = [one] + [np.asarray(n, dtype=complex) * one for n in n_list] + [one]

    m11, m12, m21, m22 = one.copy(), 0 * one, 0 * one, one.copy()

    def mul(a11, a12, a21, a22, b11, b12, b21, b22):
        return (a11 * b11 + a12 * b21, a11 * b12 + a12 * b22,
                a21 * b11 + a22 * b21, a21 * b12 + a22 * b22)

    for j in range(len(n_list) + 1):
        ni, nj = media[j], media[j + 1]
        r, tt = r_fresnel(ni, nj), t_fresnel(ni, nj)
        m11, m12, m21, m22 = mul(m11, m12, m21, m22,
                                 1.0 / tt, r / tt, r / tt, 1.0 / tt)
        if j < len(n_list):                      # propagate inside layer j+1
            beta = media[j + 1] * w * d_list[j] / C0
            m11, m12, m21, m22 = mul(m11, m12, m21, m22,
                                     np.exp(1j * beta), 0.0, 0.0, np.exp(-1j * beta))
    return 1.0 / m11


def model_transmission(freq, n_si, d_si, n_qz, d_qz):
    """
    Measured transfer function T(w) = E_sample/E_reference for the bilayer:
    the stack transmission divided by propagation of the reference pulse over
    the same geometric length of air.

        T(w) = t_stack(w) * exp(+i*w*(d_si + d_qz)/c)

    (Note: for a stack in a symmetric ambient, reciprocity makes t_stack
    invariant under reversing the layer order, so this analysis cannot - and
    need not - distinguish Si-first from quartz-first.)
    """
    t_stack = transfer_matrix_stack(freq, [n_si, n_qz], [d_si, d_qz])
    return t_stack * np.exp(1j * 2.0 * np.pi * freq * (d_si + d_qz) / C0)


# ===========================================================================
# 5.  Phase handling
# ===========================================================================
def guided_unwrap(freq, z, tau_guess_ps, phase0=0.0):
    """
    Unwrap arg(z) robustly using a known approximate group delay.

    Removing the steep linear trend -w*tau first leaves a small, slowly varying
    residual that np.unwrap cannot mis-branch even where the spectrum is noisy;
    the trend is then restored.  A global 2*pi*m branch is finally chosen so the
    extrapolation to w = 0 matches the known zero-frequency phase `phase0`.
    """
    w = 2.0 * np.pi * freq
    tau = tau_guess_ps * PS
    phase = np.unwrap(np.angle(z) + w * tau) - w * tau

    lo = (freq > 0.20 * THZ) & (freq < 0.50 * THZ)
    if lo.sum() > 4:
        intercept = np.polyfit(w[lo], phase[lo], 1)[1]
        phase -= 2.0 * np.pi * np.round((intercept - phase0) / (2.0 * np.pi))
    return phase


def optical_thickness(freq, ratio, tau_ps):
    """
    n*d [m] as a function of frequency from a gated echo ratio
    R(w) = A * exp(-2i*n~*w*d/c) with A real positive:   arg(R) = -2*n*w*d/c.
    """
    w = 2.0 * np.pi * freq
    phase = guided_unwrap(freq, ratio, tau_ps, phase0=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return -C0 * phase / (2.0 * w)


# ===========================================================================
# 6.  Inversions
# ===========================================================================
def solve_indices_from_fresnel(a_si, a_qz, bracket=(1.05, 8.0)):
    """
    Self-contained (thickness-free) index solve from the two Fresnel magnitudes

        a_si = |r12*r10| = F(n1)*G ,    a_qz = |r23*r21| = F(n2)*G
        F(n) = (n-1)/(n+1) ,            G    = (n1-n2)/(n1+n2)

    G cancels in the ratio, so a_si/a_qz = F(n1)/F(n2) expresses n2 through n1
    and reduces the 2x2 system to one well-behaved 1-D root find.

    ILL-CONDITIONED BY CONSTRUCTION: d(F)/F ~ 0.6 * d(n)/n, so amplitude errors
    are AMPLIFIED ~1.6x in the returned indices.  Cross-check only.
    """
    ratio = a_si / a_qz

    def n2_of_n1(n1):
        f2 = np.clip(((n1 - 1.0) / (n1 + 1.0)) / ratio, 1e-9, 1.0 - 1e-9)
        return (1.0 + f2) / (1.0 - f2)

    def residual(n1):
        n2 = n2_of_n1(n1)
        return (n1 - 1.0) / (n1 + 1.0) * (n1 - n2) / (n1 + n2) - a_si

    lo, hi = bracket
    if residual(lo) * residual(hi) > 0:
        return np.nan, np.nan
    n1 = brentq(residual, lo, hi, xtol=1e-10)
    return n1, n2_of_n1(n1)


def solve_k_from_transfer_matrix(freq, t_target, n_si, k_si, d_si, n_qz, d_qz,
                                 k_max=0.5, n_iter=60):
    """
    Invert the FULL Fabry-Perot transfer matrix for k_qz(w), frequency by
    frequency, so that |T_model| == |T_measured|.

    |T| is monotonically decreasing in k_qz, so a vectorised bisection is exact
    and unconditionally stable - no initial guess, no convergence failure.
    """
    lo = np.zeros_like(freq)
    hi = np.full_like(freq, k_max)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        val = np.abs(model_transmission(freq, n_si - 1j * k_si, d_si,
                                        n_qz - 1j * mid, d_qz))
        too_transmissive = val > t_target      # need more absorption
        lo = np.where(too_transmissive, mid, lo)
        hi = np.where(too_transmissive, hi, mid)
    return 0.5 * (lo + hi)


# ===========================================================================
# 7.  Main analysis
# ===========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="THz-TDS analysis of a Si / fused-quartz bilayer.")
    ap.add_argument("--folder", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--hann-full", action="store_true",
                    help="use a full-record Hann window instead of the flat-top "
                         "Hann taper (demonstrates the bias described above)")
    ap.add_argument("--out", default="thz_results.png")
    args = ap.parse_args(argv)

    taper = 1.0 if args.hann_full else TAPER_FRAC
    rule = "=" * 78
    print(rule)
    print(" THz-TDS ANALYSIS  -  silicon / fused-quartz bilayer")
    print(rule)

    # ------------------------------------------------------------- 7.1  I/O
    f_sam, f_ref = find_data_files(args.folder)
    t_s, e_s, meta_s = load_menlo(f_sam)
    t_r, e_r, meta_r = load_menlo(f_ref)
    if t_s.shape != t_r.shape or not np.allclose(t_s, t_r):
        raise ValueError("Sample and reference are not on a common time axis.")

    dt_ps = meta_r["dt_ps"]
    fs_thz = 1.0 / dt_ps
    t_main_r, t_main_s = locate_main_pulse(t_r, e_r), locate_main_pulse(t_s, e_s)

    print("\n[1] Data")
    print(f"    sample     : {meta_s['name']}")
    print(f"    reference  : {meta_r['name']}")
    print(f"    points     : {t_r.size}   span {t_r[0]:.1f} .. {t_r[-1]:.1f} ps"
          f"   dt = {dt_ps*1000:.2f} fs   averages = {meta_r.get('averages','?')}")
    print(f"    sampling   : {fs_thz:.2f} THz   (Nyquist {fs_thz/2:.2f} THz)")
    print(f"    FFT        : zero-padded to 2^14 = {NFFT} pts"
          f"   ->  df = {fs_thz*1e3/NFFT:.2f} GHz")

    # ------------------------------------------------- 7.2  spectra + T(w)
    win = raised_cosine_window(t_r.size, taper)
    freq, E_ref = spectrum(t_r, e_r, NFFT, win)
    _, E_sam = spectrum(t_s, e_s, NFFT, win)
    f_thz = freq / THZ
    w_ang = 2.0 * np.pi * freq
    with np.errstate(divide="ignore", invalid="ignore"):
        T_meas = np.divide(E_sam, E_ref, out=np.zeros_like(E_sam),
                           where=np.abs(E_ref) > 0)

    # Noise floor = high-frequency plateau of the reference spectrum (well above
    # any emitted THz power), which is a like-for-like estimate on the same FFT.
    noise_floor = float(np.median(np.abs(E_ref)[f_thz > 8.0]))
    dyn_db = 20.0 * np.log10(np.abs(E_ref) / noise_floor)
    usable = (f_thz > 0.15) & (dyn_db > 20.0)
    f_max_usable = float(f_thz[usable].max()) if usable.any() else 3.0

    win_name = ("FULL-RECORD HANN" if args.hann_full
                else f"Hann taper, flat top (taper = {taper:.2f})")
    print(f"\n[2] Spectra      window = {win_name}")
    print(f"    peak dynamic range   : {dyn_db.max():.1f} dB")
    print(f"    DR > 20 dB out to    : {f_max_usable:.2f} THz")

    w_full = raised_cosine_window(t_r.size, 1.0)
    i_r = int(np.argmin(np.abs(t_r - t_main_r)))
    i_s = int(np.argmin(np.abs(t_s - t_main_s)))
    bias = w_full[i_s] / w_full[i_r]
    print(f"    [window diagnostic] a full-record Hann weights the sample pulse "
          f"{w_full[i_s]:.3f} vs the")
    print(f"      reference pulse {w_full[i_r]:.3f}  ->  |T| biased x{bias:.3f} "
          f"({20*np.log10(bias):+.1f} dB), i.e. a")
    print(f"      spurious NEGATIVE absorption.  "
          + ("IN USE NOW (--hann-full)!" if args.hann_full
             else "Avoided by the flat-top taper."))

    # ------------------------------- 7.3  deconvolution -> impulse response
    # Dividing out the reference removes all system etalons/ringing common to
    # both traces, leaving the sample's own Fabry-Perot echo train.
    bp = band_limit(f_thz, *IRF_BAND_THZ, IRF_EDGE_THZ)
    h = np.fft.irfft(T_meas * bp, NFFT)
    t_h = np.arange(NFFT) * dt_ps               # delay relative to record start
    t_h0 = refine_peak(t_h, np.abs(h), int(np.argmax(np.abs(h))))

    tau_si, amp_si = locate_feature(t_h, h, t_h0, (6.0, 11.0))
    tau_qz, amp_qz = locate_feature(t_h, h, t_h0, (11.0, 16.0))
    tau_si2, amp_si2 = locate_feature(t_h, h, t_h0, (16.0, 19.5))
    tau_x, amp_x = locate_feature(t_h, h, t_h0, (19.5, 24.0))
    h0 = h[int(np.argmin(np.abs(t_h - t_h0)))]

    print("\n[3] Impulse response h(t) = IFFT[T(w)]  (system artifacts removed)")
    print(f"    main transmitted pulse        : {t_h0:8.4f} ps")
    print(f"    1st Si echo      (2 d_Si)     : {tau_si:8.4f} ps"
          f"   amp/main = {amp_si/h0:+.4f}")
    print(f"    1st quartz echo  (2 d_qz)     : {tau_qz:8.4f} ps"
          f"   amp/main = {amp_qz/h0:+.4f}   (sign inverted)")
    print(f"    2nd Si echo      (4 d_Si)     : {tau_si2:8.4f} ps"
          f"   amp/main = {amp_si2/h0:+.4f}")
    print(f"    Si + quartz cross echo        : {tau_x:8.4f} ps"
          f"   amp/main = {amp_x/h0:+.4f}")
    print(f"    consistency: 2 x tau_Si = {2*tau_si:.3f} ps vs observed "
          f"{tau_si2:.3f} ps   (delta {abs(2*tau_si-tau_si2)*1000:.0f} fs)")
    print(f"                tau_Si + tau_qz = {tau_si+tau_qz:.3f} ps vs observed "
          f"{tau_x:.3f} ps   (delta {abs(tau_si+tau_qz-tau_x)*1000:.0f} fs)")
    print("    Si echo keeps the sign of the main pulse (r_Si,qz > 0); the quartz")
    print("    echo is inverted (r_qz,Si < 0).  Both exactly as Fresnel predicts.")

    # --------------------------------- 7.4  gated echo ratios (per layer)
    g_main = hann_gate(t_h, t_h0, GATE_HW_PS)
    g_si = hann_gate(t_h, t_h0 + tau_si, GATE_HW_PS)
    g_qz = hann_gate(t_h, t_h0 + tau_qz, GATE_HW_PS)
    S_main = np.fft.rfft(h * g_main)
    with np.errstate(divide="ignore", invalid="ignore"):
        R_si = np.fft.rfft(h * g_si) / S_main       # Eq. (5): silicon only
        R_qz = np.fft.rfft(h * g_qz) / S_main       # Eq. (6): quartz only

    cal = (f_thz >= CAL_BAND_THZ[0]) & (f_thz <= CAL_BAND_THZ[1])
    nd_si_f = optical_thickness(freq, R_si, tau_si)
    nd_qz_f = optical_thickness(freq, -R_qz, tau_qz)     # -R_qz: r21 < 0
    nd_si = float(np.nanmedian(nd_si_f[cal]))
    nd_qz = float(np.nanmedian(nd_qz_f[cal]))

    # Main-pulse delay from the phase slope of T(w) (more precise than a peak).
    ph_T = guided_unwrap(freq, T_meas, t_h0)
    slope_band = (f_thz > 0.25) & (f_thz < 1.20)
    dt_main = -np.polyfit(w_ang[slope_band], ph_T[slope_band], 1)[0] / PS

    print(f"\n[4] Optical thicknesses from the echo phase "
          f"(Hann gates, +/-{GATE_HW_PS:.1f} ps)")
    print(f"    n*d (silicon)      = {nd_si/UM:9.2f} um")
    print(f"    n*d (quartz)       = {nd_qz/UM:9.2f} um")
    print(f"    main-pulse delay   = {dt_main:9.4f} ps   (from arg T slope)")

    # ------------------------------------------- 7.5  thickness determination
    d_total = nd_si + nd_qz - C0 * dt_main * PS            # Eq. (7), model-free
    d_si = nd_si / N_SI_REF                                # single external anchor
    d_qz = d_total - d_si
    n_si_c, n_qz_c = N_SI_REF, nd_qz / d_qz

    a1 = float(np.nanmedian(np.abs(R_si)[cal]))
    a2 = float(np.nanmedian(np.abs(R_qz)[cal]))
    n1_fr, n2_fr = solve_indices_from_fresnel(a1, a2)

    print("\n[5] Thickness determination")
    print("    model-free sum rule  d_Si + d_qz = c*(tau_Si/2 + tau_qz/2 - dt_main)")
    print(f"      = {d_total/UM:8.2f} um"
          f"     (quoted total {(D_SI_NOMINAL+D_QZ_NOMINAL)/UM:.0f} um,"
          f" {100*(d_total-D_SI_NOMINAL-D_QZ_NOMINAL)/(D_SI_NOMINAL+D_QZ_NOMINAL):+.2f} %)")
    print(f"    anchoring on n_Si = {N_SI_REF} (float-zone Si, THz standard):")
    print(f"      d_Si    = {d_si/UM:8.2f} um   "
          f"(quoted {D_SI_NOMINAL/UM:.0f} um, {100*(d_si-D_SI_NOMINAL)/D_SI_NOMINAL:+.2f} %)"
          "   <- independent check")
    print(f"      d_quartz= {d_qz/UM:8.2f} um   "
          f"(quoted {D_QZ_NOMINAL/UM:.0f} um, {100*(d_qz-D_QZ_NOMINAL)/D_QZ_NOMINAL:+.2f} %)")
    print(f"      n_quartz= {n_qz_c:8.4f}    (= n*d / d)")
    print(f"    cross-check, thickness-free Fresnel solve from |R| alone:")
    print(f"      |r12 r10| = {a1:.4f},  |r23 r21| = {a2:.4f}"
          f"  ->  n_Si = {n1_fr:.3f}, n_qz = {n2_fr:.3f}")
    print( "      (ill-conditioned by construction - amplitude errors are"
           " amplified ~1.6x; not used)")

    # ------------------------------------------ 7.6  broadband n(w) and k(w)
    band = (f_thz >= FIT_BAND_THZ[0]) & (f_thz <= min(FIT_BAND_THZ[1], f_max_usable))
    n_si_w = nd_si_f / d_si
    n_qz_w = nd_qz_f / d_qz

    # k_Si from |R_Si| against the theoretical Fresnel product.
    nsi_c, nqz_c = complex(n_si_c), complex(n_qz_c)
    A_si = abs(r_fresnel(nsi_c, nqz_c) * r_fresnel(nsi_c, 1.0))
    A_qz = abs(r_fresnel(nqz_c, 1.0) * r_fresnel(nqz_c, nsi_c))
    with np.errstate(divide="ignore", invalid="ignore"):
        k_si_w = -C0 * np.log(np.abs(R_si) / A_si) / (2.0 * w_ang * d_si)
        k_qz_echo = -C0 * np.log(np.abs(R_qz) / A_qz) / (2.0 * w_ang * d_qz)

    # Systematic floor: the offset needed to force alpha = 0 in the calibration
    # band, i.e. this method's absolute accuracy limit for a transparent layer.
    k_si_off = float(np.nanmedian(k_si_w[cal]))
    a_si_off = float(np.nanmedian((2 * w_ang * k_si_w / C0 / 100.0)[cal]))

    # k_qz from a full transfer-matrix inversion of |T(w)| - the TM model doing
    # real work.  k_Si is held at its (vanishing) calibrated value.
    k_si_fit = np.clip(k_si_w - k_si_off, 0.0, None)
    k_si_fit = np.nan_to_num(k_si_fit, nan=0.0, posinf=0.0, neginf=0.0)
    k_qz_w = solve_k_from_transfer_matrix(freq, np.abs(T_meas), n_si_c,
                                          k_si_fit, d_si, n_qz_c, d_qz)

    alpha_si = 2.0 * w_ang * (k_si_w - k_si_off) / C0 / 100.0   # [1/cm], relative
    alpha_qz = 2.0 * w_ang * k_qz_w / C0 / 100.0                # [1/cm], absolute
    alpha_qz_echo = 2.0 * w_ang * k_qz_echo / C0 / 100.0

    print("\n[6] Absorption")
    print( "    k_qz(w): full Fabry-Perot transfer-matrix inversion of |T(w)|"
           " (vectorised bisection)")
    print(f"    k_Si(w): from |R_Si| vs the Fresnel product |r12 r10| = {A_si:.4f}")
    print(f"      systematic floor of the echo-amplitude method: "
          f"|d alpha_Si| ~ {abs(a_si_off):.2f} cm^-1")
    print( "      -> high-resistivity Si absorption lies BELOW this floor;"
           " reported as an upper bound.")

    # ------------------------------------------- 7.7  transfer-matrix check
    n_si_model = n_si_c - 1j * k_si_fit
    n_qz_model = n_qz_c - 1j * np.nan_to_num(k_qz_w)
    T_model = model_transmission(freq, n_si_model, d_si, n_qz_model, d_qz)
    ph_model = guided_unwrap(freq, T_model, t_h0)

    amp_res = float(np.sqrt(np.nanmean(
        (np.abs(T_meas[band]) - np.abs(T_model[band])) ** 2)))
    ph_res = float(np.sqrt(np.nanmean((ph_T[band] - ph_model[band]) ** 2)))
    print("\n[7] Transfer-matrix validation (all Fabry-Perot orders, both layers)")
    print(f"    band {FIT_BAND_THZ[0]:.2f} - {min(FIT_BAND_THZ[1], f_max_usable):.2f} THz")
    print(f"    RMS |T| residual   : {amp_res:.4f}"
          f"   (mean |T| = {float(np.nanmean(np.abs(T_meas[band]))):.4f})")
    print(f"    RMS phase residual : {ph_res:.4f} rad")

    # ----------------------------------------------------- 7.8  summary table
    print("\n" + rule)
    print(" SUMMARY   n = refractive index,  k = extinction coefficient,")
    print("           alpha = 2*w*k/c = power absorption coefficient")
    print(rule)
    hdr = (f"{'f [THz]':>8} | {'n_Si':>7} {'k_Si':>10} {'a_Si':>9} | "
           f"{'n_qz':>7} {'k_qz':>10} {'a_qz':>9}")
    print(hdr)
    print(f"{'':>8} | {'':>7} {'':>10} {'[1/cm]':>9} | {'':>7} {'':>10} {'[1/cm]':>9}")
    print("-" * len(hdr))
    for ft in TABLE_FREQS_THZ:
        if ft > min(FIT_BAND_THZ[1], f_max_usable):
            continue
        i = int(np.argmin(np.abs(f_thz - ft)))
        print(f"{f_thz[i]:8.3f} | {n_si_w[i]:7.4f} {max(k_si_w[i]-k_si_off,0):10.2e} "
              f"{max(alpha_si[i],0):9.3f} | {n_qz_w[i]:7.4f} {k_qz_w[i]:10.2e} "
              f"{alpha_qz[i]:9.3f}")
    print("-" * len(hdr))
    print(f"{'layer':>8} | {'silicon':>28} | {'fused quartz':>28}")
    print(f"{'d':>8} | {d_si/UM:20.2f} um       | {d_qz/UM:20.2f} um")
    print(f"{'(quoted)':>8} | {D_SI_NOMINAL/UM:20.1f} um       |"
          f" {D_QZ_NOMINAL/UM:20.1f} um")
    print(f"{'n*d':>8} | {nd_si/UM:20.2f} um       | {nd_qz/UM:20.2f} um")
    print(f"{'echo':>8} | {tau_si:20.4f} ps       | {tau_qz:20.4f} ps")
    print(f"{'n (mean)':>8} | {float(np.nanmean(n_si_w[band])):20.4f}          |"
          f" {float(np.nanmean(n_qz_w[band])):20.4f}")
    print(rule)
    print(f" alpha_Si  <= {abs(a_si_off):.2f} cm^-1 across the band (upper bound;"
          " Si is transparent here)")
    print(f" alpha_qz  rises from {alpha_qz[np.argmin(np.abs(f_thz-0.5))]:.2f}"
          f" cm^-1 at 0.5 THz to "
          f"{alpha_qz[np.argmin(np.abs(f_thz-2.5))]:.2f} cm^-1 at 2.5 THz")
    print(rule)

    # ---------------------------------------------------------- 7.9  figure
    make_figure(args.out, t_r, e_r, t_s, e_s, t_main_r, t_main_s,
                t_h, h, t_h0, tau_si, tau_qz, tau_si2, tau_x,
                g_main, g_si, g_qz, f_thz, E_ref, E_sam, noise_floor,
                T_meas, T_model, ph_T, ph_model,
                n_si_w, alpha_si, n_qz_w, alpha_qz, alpha_qz_echo,
                band, d_si, d_qz, n_si_c, n_qz_c,
                min(FIT_BAND_THZ[1], f_max_usable), dyn_db,
                meta_s, meta_r, args.hann_full, taper, abs(a_si_off))
    print(f"\nFigure written to {os.path.abspath(args.out)}\n")
    return 0


# ===========================================================================
# 8.  Figure
# ===========================================================================
def _ax_style(ax, title=None):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    for lb in ax.get_xticklabels() + ax.get_yticklabels():
        lb.set_color(INK2)
    if title:
        ax.set_title(title, color=INK, fontsize=11.5, fontweight="bold",
                     loc="left", pad=8)


def make_figure(out, t_r, e_r, t_s, e_s, t_main_r, t_main_s,
                t_h, h, t_h0, tau_si, tau_qz, tau_si2, tau_x,
                g_main, g_si, g_qz, f_thz, E_ref, E_sam, noise_floor,
                T_meas, T_model, ph_T, ph_model,
                n_si_w, alpha_si, n_qz_w, alpha_qz, alpha_qz_echo,
                band, d_si, d_qz, n_si_c, n_qz_c, fmax, dyn_db,
                meta_s, meta_r, hann_full, taper, a_si_bound):
    """Seven-panel summary figure; every axis carries units."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.labelcolor": INK2, "text.color": INK, "legend.frameon": False,
    })
    fig = plt.figure(figsize=(14.0, 15.5))
    gs = GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.22,
                  left=0.068, right=0.975, top=0.918, bottom=0.048)

    fig.suptitle("THz-TDS analysis: silicon / fused-quartz bilayer",
                 x=0.068, ha="left", y=0.975, fontsize=17,
                 fontweight="bold", color=INK)
    fig.text(0.068, 0.949,
             f"Menlo ScanControl, N$_2$ purge, {meta_r.get('averages','?')} averages"
             f"   ·   {meta_s['name']}  /  {meta_r['name']}"
             f"   ·   {'full-record Hann' if hann_full else f'Hann taper, flat top ({taper:.2f})'}"
             f", zero-padded to $2^{{14}}$",
             ha="left", fontsize=9.5, color=MUTED)

    inb = band & (f_thz <= fmax)

    # --- (a) raw time domain -------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    _ax_style(ax, "(a)  Measured waveforms — sample delay and the Fabry–Pérot echo train")
    ax.plot(t_r - t_main_r, e_r, color=CLR_REF, lw=1.5, label="Reference (air)")
    ax.plot(t_s - t_main_r, e_s, color=CLR_SAM, lw=1.5, label="Sample (Si + quartz)")
    ax.set_xlim(-4, 28)
    ax.set_xlabel("Time relative to reference pulse peak  [ps]")
    ax.set_ylabel("THz field amplitude  [a.u.]")
    t0 = t_main_s - t_main_r
    ax.annotate(f"sample delayed {t0:.2f} ps",
                xy=(t0, -2.3), xytext=(t0 + 2.2, -3.1), fontsize=9.5, color=INK2,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
    for dly, lab in ((tau_si, f"Si round trip\n+{tau_si:.2f} ps"),
                     (tau_qz, f"quartz round trip\n+{tau_qz:.2f} ps"),
                     (tau_x, f"Si+quartz cross\n+{tau_x:.2f} ps")):
        ax.axvline(t0 + dly, color=MUTED, lw=0.9, ls=":")
        ax.annotate(lab, xy=(t0 + dly + 0.3, 1.15), fontsize=8.5, color=INK2)
    ax.legend(loc="upper right", fontsize=10, labelcolor=INK2)

    # --- (b) deconvolved impulse response + gates ----------------------
    ax = fig.add_subplot(gs[1, 0])
    _ax_style(ax, "(b)  Impulse response $h(t)=\\mathcal{F}^{-1}[T(\\omega)]$ and the Hann gates")
    hn = h / np.abs(h).max()
    ax.plot(t_h - t_h0, hn, color=MUTED, lw=1.4, label="$h(t)$  (system-artifact free)")
    for g, col, lab in ((g_main, "#6b6a66", "main-pulse gate"),
                        (g_si, CLR_SI, "Si-echo gate"),
                        (g_qz, CLR_QZ, "quartz-echo gate")):
        ax.fill_between(t_h - t_h0, 0, g * 0.55, color=col, alpha=0.20, lw=0)
        ax.plot(t_h - t_h0, g * 0.55, color=col, lw=1.3, ls="--")
    ax.set_xlim(-3, 24)
    ax.set_ylim(-0.45, 1.08)
    ax.set_xlabel("Delay relative to main transmitted pulse  [ps]")
    ax.set_ylabel("Normalised impulse response  [–]")
    ax.text(tau_si - 1.6, 0.62, "Si echo", color=CLR_SI, fontsize=9, fontweight="bold")
    ax.text(tau_qz - 2.2, 0.62, "quartz echo", color=CLR_QZ, fontsize=9, fontweight="bold")
    ax.text(tau_si2 - 1.0, -0.40, f"2×Si\n+{tau_si2:.1f} ps", color=INK2, fontsize=8)
    ax.text(tau_x - 1.0, -0.40, f"Si+qz\n+{tau_x:.1f} ps", color=INK2, fontsize=8)
    ax.legend(loc="upper right", fontsize=9, labelcolor=INK2)

    # --- (c) amplitude spectra -----------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    _ax_style(ax, "(c)  Amplitude spectra and dynamic range")
    ref_db = 20 * np.log10(np.abs(E_ref) / noise_floor + 1e-30)
    sam_db = 20 * np.log10(np.abs(E_sam) / noise_floor + 1e-30)
    m = (f_thz > 0.05) & (f_thz < 6.0)
    ax.plot(f_thz[m], ref_db[m], color=CLR_REF, lw=1.4, label="Reference")
    ax.plot(f_thz[m], sam_db[m], color=CLR_SAM, lw=1.4, label="Sample")
    ax.axhline(0, color=MUTED, lw=1.1, ls="--")
    ax.text(4.6, 2.5, "noise floor", fontsize=8.5, color=MUTED, ha="right")
    ax.axvspan(0.2, fmax, color=CLR_AUX, alpha=0.11, lw=0)
    ax.text(0.28, 4, "analysis band", fontsize=8.5, color="#12805a")
    ax.set_xlabel("Frequency  [THz]")
    ax.set_ylabel("Spectral amplitude  [dB above noise floor]")
    ax.set_xlim(0, 6)
    ax.set_ylim(-8, 95)
    ax.legend(loc="upper right", fontsize=9.5, labelcolor=INK2)

    # --- (d) |T| data vs model -----------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    _ax_style(ax, "(d)  Transmission magnitude — data vs transfer matrix")
    ax.plot(f_thz[inb], np.abs(T_meas[inb]), color=CLR_SAM, lw=1.6,
            label="Measured  $|T(\\omega)|$")
    ax.plot(f_thz[inb], np.abs(T_model[inb]), color=CLR_MODEL, lw=1.4, ls="--",
            label="Fabry–Pérot transfer matrix")
    ax.set_xlabel("Frequency  [THz]")
    ax.set_ylabel("$|T(\\omega)|=|E_\\mathrm{sam}/E_\\mathrm{ref}|$  [–]")
    ax.set_xlim(0.2, fmax)
    ax.legend(loc="lower left", fontsize=9.5, labelcolor=INK2)

    # --- (e) phase data vs model ---------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    _ax_style(ax, "(e)  Transmission phase — data vs transfer matrix")
    ax.plot(f_thz[inb], ph_T[inb], color=CLR_SAM, lw=1.6, label="Measured  $\\arg T$")
    ax.plot(f_thz[inb], ph_model[inb], color=CLR_MODEL, lw=1.4, ls="--",
            label="Fabry–Pérot transfer matrix")
    ax.set_xlabel("Frequency  [THz]")
    ax.set_ylabel("Unwrapped phase  $\\arg T(\\omega)$  [rad]")
    ax.set_xlim(0.2, fmax)
    ax.legend(loc="lower left", fontsize=9.5, labelcolor=INK2)

    # --- (f) refractive index ------------------------------------------
    ax = fig.add_subplot(gs[3, 0])
    _ax_style(ax, "(f)  Refractive index of each layer")
    ax.plot(f_thz[inb], n_si_w[inb], color=CLR_SI, lw=2.0, label="Silicon")
    ax.plot(f_thz[inb], n_qz_w[inb], color=CLR_QZ, lw=2.0, label="Fused quartz")
    ax.set_xlabel("Frequency  [THz]")
    ax.set_ylabel("Refractive index  $n(\\omega)$  [–]")
    ax.set_xlim(0.2, fmax)
    ax.set_ylim(1.7, 3.7)
    ax.text(fmax * 0.50, n_si_c + 0.12, f"Si,  $n = {n_si_c:.4f}$",
            color=CLR_SI, fontsize=10, fontweight="bold")
    ax.text(fmax * 0.50, n_qz_c + 0.12, f"fused quartz,  $n = {n_qz_c:.4f}$",
            color=CLR_QZ, fontsize=10, fontweight="bold")
    ax.legend(loc="center left", fontsize=9.5, labelcolor=INK2)

    # --- (g) absorption -------------------------------------------------
    ax = fig.add_subplot(gs[3, 1])
    _ax_style(ax, "(g)  Power absorption coefficient  $\\alpha = 2\\omega k/c$")
    ax.plot(f_thz[inb], np.clip(alpha_qz[inb], 1e-3, None), color=CLR_QZ, lw=2.0,
            label="Fused quartz  (transfer-matrix inversion of $|T|$)")
    ax.plot(f_thz[inb], np.clip(alpha_qz_echo[inb], 1e-3, None), color=CLR_QZ,
            lw=1.1, ls=":", label="Fused quartz  (quartz-echo ratio, independent)")
    ax.axhspan(1e-3, a_si_bound, color=CLR_SI, alpha=0.16, lw=0)
    ax.axhline(a_si_bound, color=CLR_SI, lw=1.6, ls="--",
               label=f"Silicon: upper bound  $\\alpha \\leq$ {a_si_bound:.1f} cm$^{{-1}}$")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency  [THz]")
    ax.set_ylabel("Absorption coefficient  $\\alpha$  [cm$^{-1}$]")
    ax.set_xlim(0.2, fmax)
    ax.set_ylim(1e-2, 60)
    ax.legend(loc="lower right", fontsize=8.5, labelcolor=INK2)
    ax.text(0.03, 0.93,
            f"$d_\\mathrm{{Si}}$ = {d_si/UM:.1f} µm     "
            f"$d_\\mathrm{{quartz}}$ = {d_qz/UM:.1f} µm",
            transform=ax.transAxes, fontsize=9.5, color=INK2)

    fig.savefig(out, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
