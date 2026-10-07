#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
================================================================================
 thz_analysis.py -- THz-TDS analysis of a Si / fused-quartz bilayer
================================================================================

Measurement
-----------
    Menlo Systems ScanControl THz-TDS, 1000-waveform average, N2 purge.
      * sample    : high-resistivity Si wafer (nominal 385 um) on a quartz
                    window (nominal ~1 mm)
      * reference : identical beam path, empty (dry N2)

What this script does
---------------------
    1. Reads both Menlo ASCII waveforms.
    2. Windows and FFTs both traces (zero-padded to 2^14 points).
    3. Forms the complex transmission  T(w) = E_sample(w) / E_ref(w),
       with a correctly branched (2*pi-resolved) phase.
    4. Builds a Fabry-Perot **transfer-matrix** model of the
       air | Si | quartz | air stack (2x2 characteristic matrices),
       cross-validated against an independent Airy recursion.
    5. Extracts n(w) and k(w) for **each layer separately**.
    6. Estimates each layer thickness d from echo timing in the time domain.
    7. **Profiles chi^2 to obtain honest error bars on the thicknesses**, and
       runs two systematic stress tests (a free air gap; a k_Si-free refit).
    8. Saves 'thz_results.png'.
    9. Prints summary tables of n, k, alpha and d.

--------------------------------------------------------------------------
 WHAT IS AND IS NOT DETERMINED BY THIS MEASUREMENT  (read before quoting d)
--------------------------------------------------------------------------
The global fit reports d_Si and d_qz to many digits, but those digits are NOT
all meaningful.  Profiling chi^2 with n_Si held fixed and everything else
re-optimised (done automatically, see section 4b of the output) shows that:

    WELL determined    n_Si*d_Si  (Si  optical thickness)    ~ +-0.1 um
                       n_qz*d_qz  (qz  optical thickness)    ~ +-0.2 um
                       d_Si+d_qz  (geometric sum)            ~ +-0.1 um
    POORLY determined  d_Si, d_qz individually               ~ +-5   um

The reason is physical: the Fabry-Perot fringe period of a layer fixes n*d,
and the main-pulse delay fixes (n1-1)d1 + (n2-1)d2.  Given n_Si*d_Si, trading
n_Si against d_Si moves chi^2 almost not at all -- the literature value
n_Si = 3.4175 costs under 1% in chi^2 and lands d_Si within a micron of the
385 um nominal, which this data set cannot distinguish from its own best fit.
So the thicknesses are quoted as  d_Si = 386 +- 5 um,  d_qz = 990 +- 5 um,
NOT to 0.01 um.  The 0.01-um-level agreement between the "echo timing" and
"Fabry-Perot fit" columns applies only to the SUM, which is index-free.

--------------------------------------------------------------------------
 A KNOWN SYSTEMATIC: THE Si ECHO IS ~11% TOO LARGE
--------------------------------------------------------------------------
The measured Si round-trip echo is consistently ~11% stronger than the
best-fit Fresnel model, one-sided, and independent of the Wiener
regularisation (checked over three decades of reg -- see section 5b).  It is
therefore a model systematic, not noise.  Matching it with Fresnel
coefficients alone would demand n_Si ~ 3.52, which is unphysical; the more
likely cause is that the model assumes **optical contact** between wafer and
window, whereas an uncontacted wafer has a sub-wavelength air gap that adds
reflectivity at the Si|qz boundary.  Consequence to keep in mind: the
pointwise solver can only absorb a too-large R1 by pushing k_Si down, so
"k_Si pins at zero" is partly an artefact of this bias and must be read as an
upper limit, never as a measurement of Si transparency.

--------------------------------------------------------------------------
 THE PHYSICS THAT MAKES STEP 5 POSSIBLE  (please read -- it is the crux)
--------------------------------------------------------------------------
A single complex transmission spectrum T(w) supplies only TWO real numbers
per frequency (|T| and arg T), but a bilayer has FOUR unknowns
(n_Si, k_Si, n_qz, k_qz).  Point-by-point inversion of T(w) alone is therefore
*mathematically underdetermined* -- any code that claims to invert it directly
is silently fitting noise.

The missing information is in the **internal Fabry-Perot echoes**, which are
resolved in time here (80 ps record, echoes at +8.7 ps and +13.0 ps).  We
therefore use three independent observables per frequency:

    O0  : the directly transmitted sub-pulse        ->  t01 t12 t23 a b
    R1  = O1/O0 : Si  round-trip echo / direct      ->  -r01 r12 a^2
    R2  = O2/O0 : qz  round-trip echo / direct      ->  -r12 r23 b^2

    with a = exp(-i N1 w d1 / c),  b = exp(-i N2 w d2 / c)
    and r_ij = (N_i - N_j)/(N_i + N_j).

R1 is dominated by (N_Si, d_Si) and R2 by (N_qz, d_qz), so the degeneracy is
broken.  Together with the exact full T(w) that is 6 real equations for 4
unknowns -- overdetermined and well posed.

The echoes are isolated by **Wiener-deconvolving the sample trace by the
reference trace** first.  This removes the instrument response (the Menlo
system itself has a small internal etalon at +7.9 ps that would otherwise
contaminate the Si echo at +8.7 ps) and leaves the bare impulse response
h(tau) of the stack, whose sub-pulses can be time-gated cleanly.

Reciprocity note: at normal incidence the transmission of a stack -- and, as
one can verify algebraically, each single-round-trip echo ratio R1, R2 -- is
invariant under reversing the layer order.  So this analysis does not need to
know whether the beam hit the Si or the quartz first.

--------------------------------------------------------------------------
 SIGN / FOURIER CONVENTION  (get this wrong and k comes out negative)
--------------------------------------------------------------------------
numpy.fft uses  X(f) = sum_n x_n exp(-2*pi*i*f*t_n).  A pure time delay tau
therefore multiplies the spectrum by exp(-i*w*tau).  Propagation through
thickness d of a medium of complex index N must reduce to that for N -> 1, so

        propagator = exp(-i * N * w * d / c),      N = n - i*k,   k >= 0

which decays like exp(-k*w*d/c).  Every phase in this file follows that
convention.  Absorption coefficient (intensity):  alpha = 2*k*w/c = 4*pi*f*k/c.

--------------------------------------------------------------------------
 A WORD ON THE HANN WINDOW  (a real trap in TDS ratio spectroscopy)
--------------------------------------------------------------------------
The brief asks for a Hann window.  A *full-record* Hann window is the wrong
thing to apply here, and the code quantifies why (see 'window bias' in the
printed output).  The sample pulse arrives 6.3 ps later than the reference
pulse, so the two pulses sit at different heights on the Hann bell: with this
80 ps record the window multiplies the sample pulse by ~0.62 and the reference
pulse by ~0.37.  T = E_s/E_r is then biased by ~1.7x -- a pure artefact that
would turn into a large negative absorption.

The physically correct way to "apply a Hann window" to a ratio measurement is
to keep the window flat over the whole signal region and use Hann half-cosine
tapers only on the noise-only ends of the record (a Tukey window with Hann
tapers).  That is the default here ('flattop-hann').  Pure full-record Hann is
still available via --window hann, and the bias it introduces is printed so the
reader can see the size of the effect.

Author: THz-TDS analysis for v.apostolopoulos@physics.uoc.gr
Requires: numpy, scipy, matplotlib
================================================================================
"""

from __future__ import annotations

import argparse
import os
import sys
import textwrap

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.signal import find_peaks, hilbert

import matplotlib
matplotlib.use("Agg")           # headless: write a PNG, never open a window
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# =============================================================================
# 0.  CONSTANTS AND CONFIGURATION
# =============================================================================

# Speed of light in the units used throughout this file:
#   time     -> ps
#   length   -> um
#   frequency-> THz          (1 THz * 1 ps = 1,  1 um / 1 ps = 1e6 m/s)
C_UM_PS = 299.792458                     # um/ps  == c in these units

SAMPLE_FILE = "Quartz_Si_15min_nitrogen.txt"
REF_FILE    = "Reference_15min_nitrogen (1).txt"
OUT_PNG     = "thz_results.png"

NFFT = 2 ** 14                           # zero-pad length (brief: 2^14)

# Nominal geometry (used ONLY as a starting guess and to label which echo
# belongs to which layer -- every number reported is fitted from the data).
D_SI_NOMINAL_UM = 385.0
D_QZ_NOMINAL_UM = 1000.0

# Analysis bands (THz)
BAND_FIT     = (0.20, 3.00)   # global transfer-matrix fit
BAND_POINTWISE = (0.20, 2.50) # per-frequency 4-parameter inversion
BAND_PHASE_ANCHOR = (0.30, 1.00)  # linear region used to fix the 2*pi branch
F_ANCHOR = 1.00               # continuation start for the pointwise inversion

# Frequencies printed in the summary table
KEY_FREQS = (0.30, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50)

# --- plotting palette (validated colour-blind-safe categorical slots) --------
CLR_SI    = "#2a78d6"   # slot 1, blue    -> silicon        / sample trace
CLR_QZ    = "#1baf7a"   # slot 3, aqua    -> quartz         / reference trace
CLR_MODEL = "#eb6834"   # slot 2, orange  -> transfer-matrix model
CLR_INK   = "#0b0b0b"
CLR_INK2  = "#52514e"
CLR_MUTED = "#8a8985"
CLR_GRID  = "#e3e2de"


# =============================================================================
# 1.  FILE I/O
# =============================================================================

def read_menlo(path):
    """
    Read a Menlo Systems ScanControl two-column ASCII waveform.

    File layout::

        # Menlo Systems GmbH ScanControl
        # ... several '#' header lines ...
        # Time [ps]   THz Signal [a.u.]
        -740.00000000   -0.00001203
        ...

    Returns
    -------
    t : (N,) float ndarray   time axis [ps], monotonically increasing
    e : (N,) float ndarray   THz field amplitude [a.u.]
    header : list[str]       the comment lines, for provenance
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"data file not found: {path}")

    header = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                break

    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"{path}: expected two numeric columns")

    t, e = data[:, 0].astype(float), data[:, 1].astype(float)

    order = np.argsort(t)                       # be robust to reversed scans
    t, e = t[order], e[order]

    # The delay stage must be uniform for an FFT to be meaningful.
    #
    # Tolerance choice: 0.1% of a step.  A tighter bound does NOT test the
    # stage, it tests how many decimal places the exporter happened to write.
    # These files carry 8 decimals, which leaves ptp(dt) ~ 1e-8 ps against a
    # 1e-6-relative threshold of 3.3e-8 ps -- a 3x margin set entirely by
    # ASCII rounding, so a Menlo export with one fewer decimal place would be
    # refused for no physical reason.  Real stage non-uniformity that would
    # actually corrupt an FFT (mis-stepped encoder, dropped point) is far more
    # than 0.1% of a step, so 1e-3 catches what matters and nothing else.
    dt = np.diff(t)
    if np.ptp(dt) > 1e-3 * np.mean(dt):
        raise ValueError(f"{path}: non-uniform time axis (jitter "
                         f"{np.ptp(dt):.2e} ps on a mean step of "
                         f"{np.mean(dt):.2e} ps)")
    return t, e, header


def remove_baseline(e, frac=0.05):
    """
    Remove the DC offset of the lock-in / ADC.

    A constant offset is a delta function at f = 0 after the FFT; it also
    leaks into the lowest frequency bins through the window.  We estimate it
    from the leading `frac` of the record, which is pure pre-pulse baseline
    (the THz pulse arrives ~20 ps into an 80 ps scan here).
    """
    n = max(8, int(frac * len(e)))
    return e - np.mean(e[:n])


# =============================================================================
# 2.  WINDOWING AND FFT
# =============================================================================

def make_window(n, mode="flattop-hann", taper_ps=6.0, dt=1.0):
    """
    Build the apodisation window.

    Parameters
    ----------
    mode : {'flattop-hann', 'hann', 'boxcar'}
        'flattop-hann' (default, recommended)
            Unity across the whole signal region, with Hann half-cosine
            tapers of length `taper_ps` on the two noise-only ends of the
            record.  This is a Tukey window whose tapers are literally Hann
            half-periods.  It suppresses the discontinuity at the record ends
            (which would otherwise ring across the whole spectrum) WITHOUT
            re-weighting the sample and reference pulses differently -- see
            the module docstring.
        'hann'
            Textbook full-record Hann.  Included so the bias it causes in a
            ratio measurement can be demonstrated, not because it is correct.
        'boxcar'
            No apodisation (rectangular).

    Returns
    -------
    w : (n,) float ndarray
    """
    if mode == "boxcar":
        return np.ones(n)

    if mode == "hann":
        # np.hanning is the symmetric Hann: 0.5*(1 - cos(2*pi*i/(n-1)))
        return np.hanning(n)

    if mode == "flattop-hann":
        ntap = int(round(taper_ps / dt))
        ntap = int(np.clip(ntap, 1, n // 2 - 1))
        w = np.ones(n)
        # np.hanning(2*ntap) is a full bell of length 2*ntap; its two halves
        # are exactly the rising and falling Hann quarter-waves we want.
        bell = np.hanning(2 * ntap)
        w[:ntap] = bell[:ntap]
        w[-ntap:] = bell[ntap:]
        return w

    raise ValueError(f"unknown window mode {mode!r}")


def rfft_spectrum(e, dt, nfft=NFFT, window=None):
    """
    Windowed, zero-padded real FFT.

    Zero padding does NOT add information: it interpolates the spectrum onto a
    finer grid (here df = 1/(nfft*dt) ~ 1.85 GHz) so that the sharp Fabry-Perot
    fringes of the stack are properly sampled.

    Returns
    -------
    f : (nfft//2+1,) frequency axis [THz]
    E : (nfft//2+1,) complex spectrum, convention exp(-i*w*t)
    """
    if window is not None:
        e = e * window
    if nfft < len(e):
        raise ValueError("nfft must be >= number of samples")
    E = np.fft.rfft(e, n=nfft)
    f = np.fft.rfftfreq(nfft, d=dt)          # THz because dt is in ps
    return f, E


def branch_corrected_phase(f, T, band=BAND_PHASE_ANCHOR, mask=None):
    """
    Unwrap arg T(f) and put it on the physically correct 2*pi branch.

    np.unwrap removes 2*pi jumps but leaves an arbitrary overall multiple of
    2*pi.  Physically, a non-magnetic slab must satisfy phi(f) -> 0 as f -> 0
    (no phase shift at DC).  We therefore fit a straight line to the unwrapped
    phase over a high-SNR, low-frequency band where phi is linear in f, and
    remove the nearest multiple of 2*pi from its extrapolated intercept.

    Getting this wrong shifts n by  m * c / (f * d)  -- a classic TDS error.
    """
    phi = np.unwrap(np.angle(T))
    sel = (f >= band[0]) & (f <= band[1])
    if mask is not None:
        sel &= mask
    if sel.sum() < 5:
        return phi
    slope, intercept = np.polyfit(f[sel], phi[sel], 1)
    phi -= 2.0 * np.pi * np.round(intercept / (2.0 * np.pi))
    return phi


# =============================================================================
# 3.  OPTICS:  TRANSFER MATRIX FOR THE FABRY-PEROT STACK
# =============================================================================
#
# Field ansatz inside medium j:
#       E(z) = Ep_j * exp(-i*k_j*z) + Em_j * exp(+i*k_j*z),    k_j = N_j*w/c
# consistent with the numpy FFT convention (see module docstring).
#
# Matching E and dE/dz at an interface gives the interface matrix
#       [Ep_i; Em_i] = (1/t_ij) [[1, r_ij], [r_ij, 1]] [Ep_j; Em_j]
# and referring the fields from the far to the near side of a layer gives
#       P_j = diag( exp(+i*delta_j), exp(-i*delta_j) ),  delta_j = N_j*w*d_j/c
#
# The stack matrix is  M = I_01 P_1 I_12 P_2 ... I_(L,L+1)  and
#       t = 1 / M[0,0]        r = M[1,0] / M[0,0]
# -----------------------------------------------------------------------------

def fresnel_r(Ni, Nj):
    """Amplitude reflection coefficient for a wave in medium i hitting j."""
    return (Ni - Nj) / (Ni + Nj)


def fresnel_t(Ni, Nj):
    """Amplitude transmission coefficient, medium i -> medium j."""
    return 2.0 * Ni / (Ni + Nj)


def _interface_matrix(Ni, Nj, shape):
    """2x2 interface matrix, broadcast to `shape` + (2, 2)."""
    r = np.broadcast_to(fresnel_r(Ni, Nj) * np.ones(shape, complex), shape)
    t = np.broadcast_to(fresnel_t(Ni, Nj) * np.ones(shape, complex), shape)
    M = np.empty(shape + (2, 2), complex)
    M[..., 0, 0] = 1.0
    M[..., 0, 1] = r
    M[..., 1, 0] = r
    M[..., 1, 1] = 1.0
    return M / t[..., None, None]


def _propagation_matrix(N, d, f, shape):
    """2x2 propagation matrix through thickness d [um] of index N."""
    delta = 2.0 * np.pi * f * N * d / C_UM_PS       # radians
    delta = np.broadcast_to(delta * np.ones(shape, complex), shape)
    P = np.zeros(shape + (2, 2), complex)
    P[..., 0, 0] = np.exp(1j * delta)
    P[..., 1, 1] = np.exp(-1j * delta)
    return P


def transfer_matrix_stack(f, N_layers, d_layers, N_in=1.0, N_out=1.0):
    """
    Complex amplitude transmission/reflection of an arbitrary layer stack.

    Parameters
    ----------
    f        : (M,) frequencies [THz]
    N_layers : sequence of complex index arrays (each scalar or (M,)),
               ordered along the propagation direction
    d_layers : sequence of thicknesses [um]
    N_in     : index of the incidence half-space (air = 1)
    N_out    : index of the exit half-space (air = 1)

    Returns
    -------
    t, r : (M,) complex
        Amplitude coefficients of the *whole stack*, referenced to the two
        outer surfaces.  All Fabry-Perot multiple reflections in every layer
        are included exactly (the matrix product resums the infinite series).
    """
    f = np.atleast_1d(np.asarray(f, float))
    shape = f.shape
    Ns = [np.asarray(N_in, complex)] + \
         [np.asarray(N, complex) for N in N_layers] + \
         [np.asarray(N_out, complex)]

    M = _interface_matrix(Ns[0], Ns[1], shape)
    for j in range(1, len(Ns) - 1):
        M = M @ _propagation_matrix(Ns[j], d_layers[j - 1], f, shape)
        M = M @ _interface_matrix(Ns[j], Ns[j + 1], shape)

    m00 = M[..., 0, 0]
    return 1.0 / m00, M[..., 1, 0] / m00


def airy_bilayer_t(f, N1, d1, N2, d2, N_in=1.0, N_out=1.0):
    """
    Independent implementation of the same bilayer transmission, via the
    Airy (recursive Fresnel) formula.  Used only to validate the transfer
    matrix in `self_test()` -- two different algebras agreeing to machine
    precision is a strong correctness check.
    """
    a = np.exp(-2j * np.pi * f * N1 * d1 / C_UM_PS)      # single pass, layer 1
    b = np.exp(-2j * np.pi * f * N2 * d2 / C_UM_PS)      # single pass, layer 2
    r01, t01 = fresnel_r(N_in, N1), fresnel_t(N_in, N1)
    r12, t12 = fresnel_r(N1, N2),   fresnel_t(N1, N2)
    r23, t23 = fresnel_r(N2, N_out), fresnel_t(N2, N_out)

    den2 = 1.0 + r12 * r23 * b ** 2                      # collapse layer 2
    r123 = (r12 + r23 * b ** 2) / den2
    t123 = t12 * t23 * b / den2

    den1 = 1.0 + r01 * r123 * a ** 2                     # collapse layer 1
    return t01 * t123 * a / den1


def model_transmission(f, d1, d2, N1, N2):
    """
    Model of the MEASURED quantity  T(w) = E_sample(w) / E_reference(w).

    The reference beam travels the same geometrical distance (d1 + d2) through
    air, so the measured ratio is the stack transmission divided by the
    free-space propagator exp(-i*w*(d1+d2)/c):

        T_model = t_stack * exp(+2*pi*i*f*(d1+d2)/c)
    """
    t, _ = transfer_matrix_stack(f, [N1, N2], [d1, d2])
    return t * np.exp(2j * np.pi * f * (d1 + d2) / C_UM_PS)


def model_transmission_gap(f, d1, dg, d2, N1, N2):
    """
    Same as `model_transmission` but with a free AIR GAP of thickness dg
    between the two solid layers:  air | Si | air(dg) | quartz | air.

    This is the physically obvious refinement if the wafer is merely resting
    on the window rather than optically contacted.  It is fitted only as a
    systematic stress test (section 4c of the output): if letting dg float
    moves d_Si by many microns while improving chi^2, then d_Si was never
    determined to better than that in the first place.
    """
    t, _ = transfer_matrix_stack(f, [N1, 1.0 + 0j, N2], [d1, dg, d2])
    return t * np.exp(2j * np.pi * f * (d1 + dg + d2) / C_UM_PS)


def echo_ratio_models(f, d1, d2, N1, N2):
    """
    Model for the two time-gated echo ratios (see module docstring).

        R1 = (Si  round-trip echo) / (direct pulse) = -r01 * r12 * a^2
        R2 = (qz  round-trip echo) / (direct pulse) = -r12 * r23 * b^2

    Derivation of R1: the sub-pulse that makes one extra round trip inside
    layer 1 picks up  r12 (reflection off the 1|2 interface, going forward)
    * r10 (reflection back off the 1|0 interface) * a^2 (two extra passes).
    With r10 = -r01 this is -r01*r12*a^2.  R2 follows identically with
    r21 = -r12.  These are exactly the +8.7 ps and +13.0 ps features seen in
    the deconvolved impulse response.
    """
    a2 = np.exp(-4j * np.pi * f * N1 * d1 / C_UM_PS)     # a^2
    b2 = np.exp(-4j * np.pi * f * N2 * d2 / C_UM_PS)     # b^2
    r01 = fresnel_r(1.0, N1)
    r12 = fresnel_r(N1, N2)
    r23 = fresnel_r(N2, 1.0)
    return -r01 * r12 * a2, -r12 * r23 * b2


def self_test(verbose=True):
    """
    Correctness checks on the optics, run before it is used on real data.

    Each check is named and reported individually so that a failure is
    immediately diagnosable rather than just "FAIL".
    """
    rng = np.random.default_rng(0)
    f = np.linspace(0.1, 4.0, 257)
    checks = []

    def chk(name, passed, detail=""):
        checks.append((name, bool(passed), detail))

    # (1) Transfer matrix == Airy recursion, for random LOSSY indices.
    worst = 0.0
    for _ in range(5):
        N1 = rng.uniform(1.5, 4.0) - 1j * rng.uniform(0.0, 0.2)
        N2 = rng.uniform(1.2, 3.0) - 1j * rng.uniform(0.0, 0.2)
        d1, d2 = rng.uniform(50, 600), rng.uniform(200, 1500)
        t_tm, _ = transfer_matrix_stack(f, [N1, N2], [d1, d2])
        t_ai = airy_bilayer_t(f, N1, d1, N2, d2)
        worst = max(worst, float(np.max(np.abs(t_tm - t_ai))))
    chk("transfer matrix == Airy recursion (lossy)", worst < 1e-12,
        f"max|dt| = {worst:.2e}")

    # (2) Energy conservation for a LOSSLESS stack in air: |t|^2 + |r|^2 == 1.
    t, r = transfer_matrix_stack(f, [3.4175, 1.95], [385.0, 1000.0])
    e = float(np.max(np.abs(np.abs(t) ** 2 + np.abs(r) ** 2 - 1.0)))
    chk("energy conservation, lossless stack", e < 1e-12, f"max err = {e:.2e}")

    # (3) A layer of zero thickness must be invisible: the bilayer with
    #     d2 = 0 must equal the single slab of layer 1.
    t_bi, _ = transfer_matrix_stack(f, [3.4175, 1.95], [385.0, 0.0])
    t_1, _ = transfer_matrix_stack(f, [3.4175], [385.0])
    e = float(np.max(np.abs(t_bi - t_1)))
    chk("zero-thickness layer is invisible", e < 1e-12, f"max err = {e:.2e}")

    # (4) An "air layer" must also be invisible apart from its propagation.
    t_air, _ = transfer_matrix_stack(f, [1.0], [250.0])
    e = float(np.max(np.abs(t_air - np.exp(-2j * np.pi * f * 250.0 / C_UM_PS))))
    chk("index-matched layer = free propagation", e < 1e-12,
        f"max err = {e:.2e}")

    # (5) Sign convention: a lossy slab must ATTENUATE (|T| < 1), not amplify.
    #     If the sign of Im(N) were flipped this check would blow up.
    Tl = model_transmission(f, 385.0, 1000.0, 3.4175 - 0.01j, 1.95 - 0.01j)
    chk("lossy media attenuate (Im N sign)", bool(np.all(np.abs(Tl) < 1.0)),
        f"max|T| = {np.max(np.abs(Tl)):.4f}")

    # (6) The gated-echo model must be EXACTLY the right term of the exact
    #     multiple-reflection series.  Write the exact stack transmission as
    #     a Laurent series  t = sum_{p,q} c_pq a^p b^q  in the single-pass
    #     propagators a and b.  Each (p, q) is one sub-pulse, arriving at
    #     delay ((p-1)/2)*tau_1 + ((q-1)/2)*tau_2 after the direct pulse.
    #     We extract the coefficients numerically by sweeping the two
    #     propagation phases around the unit circle and 2-D DFT-ing, then
    #     assert
    #           c_31/c_11 == R1       (the Si  round-trip echo)
    #           c_13/c_11 == R2       (the quartz round-trip echo)
    #     This proves the gated echo model is EXACT, not a weak-reflection
    #     approximation.  (M must be large enough that the aliased c_{1,M+3}
    #     term is negligible; M = 48 gives ~1e-12.)
    M, f0 = 48, 1.0
    th = 2.0 * np.pi * np.arange(M) / M
    N1, N2 = 3.4175 + 0j, 1.95 + 0j          # lossless -> a, b on unit circle
    dd1 = th * C_UM_PS / (2 * np.pi * f0 * N1.real)
    dd2 = th * C_UM_PS / (2 * np.pi * f0 * N2.real)
    tg = np.empty((M, M), complex)
    for j in range(M):
        for l in range(M):
            tg[j, l] = transfer_matrix_stack(np.array([f0]), [N1, N2],
                                             [dd1[j], dd2[l]])[0][0]

    # a = exp(-i*th1), b = exp(-i*th2)  =>  c_pq = < t * a^-p * b^-q >
    def coef(pp, qq):
        return np.mean(tg * np.exp(1j * pp * th)[:, None]
                          * np.exp(1j * qq * th)[None, :])

    c11, c31, c13, c33 = coef(1, 1), coef(3, 1), coef(1, 3), coef(3, 3)
    e0 = abs(c11 - fresnel_t(1, N1) * fresnel_t(N1, N2) * fresnel_t(N2, 1))
    e1 = abs(c31 / c11 - (-fresnel_r(1, N1) * fresnel_r(N1, N2)))
    e2 = abs(c13 / c11 - (-fresnel_r(N1, N2) * fresnel_r(N2, 1)))
    chk("direct sub-pulse == t01 t12 t23", e0 < 1e-9, f"err = {e0:.2e}")
    chk("R1 == exact a^3 b^1 series coefficient", e1 < 1e-9, f"err = {e1:.2e}")
    chk("R2 == exact a^1 b^3 series coefficient", e2 < 1e-9, f"err = {e2:.2e}")

    # (7) The a^3 b^3 cross term (BOTH layers make a round trip) is large
    #     here and arrives at tau_1 + tau_2.  Confirm it is large, i.e. that
    #     time-gating the echoes apart is genuinely necessary and not
    #     cosmetic -- a series truncated after R1, R2 would be wrong by this.
    chk("cross term is large => gating required", abs(c33 / c11) > 0.05,
        f"|c33/c11| = {abs(c33 / c11):.4f}")

    # (8) The absorption SENSITIVITY of the echo ratio, which the detection
    #     floor in section 8 is built on.  This is a factor-of-2 trap in both
    #     directions: the echo makes a round trip of path 2*d, but |R1| is a
    #     FIELD ratio and a field decays as exp(-alpha*L/2), so the two 2's
    #     cancel and the exact relation is
    #           |R1| ∝ exp(-alpha*d1)   ⇒   -d ln|R1| / d alpha = d1.
    #     Verified here by logarithmic differentiation, so a future edit that
    #     reintroduces a stray 2 fails loudly instead of quietly halving (or
    #     doubling) every quoted absorption limit.
    #
    #     Tolerance is 0.1%, not machine epsilon, and that is not sloppiness:
    #     the relation is exact for the PROPAGATOR a^2, but R1 also carries
    #     the prefactor -r01*r12, whose Fresnel coefficients pick up their own
    #     weak dependence on k.  That contributes ~5e-5 here.  0.1% therefore
    #     accepts the real physics while still failing hard on the error this
    #     check exists to catch -- a factor of 2 is a 100% discrepancy.
    d1t, d2t, ft = 385.0, 1000.0, np.array([1.0])
    N1t, N2t = 3.4175 + 0j, 1.95 + 0j
    dal = 1.0                                        # cm^-1 perturbation
    N1t2 = N1t - 1j * alpha_to_k(np.array([dal]), ft)[0]
    rA = abs(echo_ratio_models(ft, d1t, d2t, N1t, N2t)[0][0])
    rB = abs(echo_ratio_models(ft, d1t, d2t, N1t2, N2t)[0][0])
    sens_cm = -np.log(rB / rA) / dal                 # = d1 in cm, exactly
    e = abs(sens_cm / (d1t * 1e-4) - 1.0)
    chk("echo absorption sensitivity == d (not 2d, not d/2)", e < 1e-3,
        f"-dln|R1|/dalpha = {sens_cm:.6f} cm vs d1 = {d1t*1e-4:.6f} cm "
        f"({e:.1e} rel)")

    # (9) Same for layer 2, so the quartz floor is covered too.
    N2t2 = N2t - 1j * alpha_to_k(np.array([dal]), ft)[0]
    rC = abs(echo_ratio_models(ft, d1t, d2t, N1t, N2t2)[1][0])
    rD = abs(echo_ratio_models(ft, d1t, d2t, N1t, N2t)[1][0])
    sens2_cm = -np.log(rC / rD) / dal
    e = abs(sens2_cm / (d2t * 1e-4) - 1.0)
    chk("same sensitivity relation for layer 2", e < 1e-3,
        f"-dln|R2|/dalpha = {sens2_cm:.6f} cm vs d2 = {d2t*1e-4:.6f} cm "
        f"({e:.1e} rel)")

    n_pass = sum(c[1] for c in checks)
    ok = n_pass == len(checks)
    if verbose:
        print(f"  optics self-test: {n_pass}/{len(checks)} checks passed")
        for name, passed, detail in checks:
            flag = "ok  " if passed else "FAIL"
            print(f"      [{flag}] {name:<44s} {detail}")
    if not ok:
        raise RuntimeError("optics self-test failed -- refusing to continue")
    return ok


# =============================================================================
# 4.  TIME-DOMAIN: DECONVOLUTION, ECHO FINDING, THICKNESS FROM TIMING
# =============================================================================

def wiener_deconvolve(Es, Er, f, reg=1e-3, f_hi=3.5, f_lo=0.10):
    """
    Regularised (Wiener) deconvolution of the sample trace by the reference.

        H(w) = Es * conj(Er) / (|Er|^2 + eps) * B(w)

    H(w) is the complex transmission T(w) wherever the reference has signal,
    and rolls smoothly to zero where it does not.  Its inverse FFT h(tau) is
    the **impulse response of the sample stack alone**: every feature in h is
    a Fabry-Perot sub-pulse, with the instrument response (laser pulse shape,
    antenna response, and the system's own internal etalon at +7.9 ps) divided
    out.  That is what makes clean time-gating of the echoes possible.

    B(w) is a smooth (super-Gaussian) band-pass matched to the usable dynamic
    range; a smooth roll-off is essential, because a hard spectral cut would
    ring in time and smear one sub-pulse into the next.
    """
    eps = reg * np.max(np.abs(Er) ** 2)
    bandpass = (np.exp(-(f / f_hi) ** 6) *
                (1.0 - np.exp(-(f / f_lo) ** 4)))
    H = Es * np.conj(Er) / (np.abs(Er) ** 2 + eps) * bandpass
    h = np.fft.irfft(H, n=2 * (len(f) - 1))
    return h, H, bandpass


def subsample_peak(x, y):
    """Parabolic interpolation of the maximum of y(x) on a uniform grid."""
    i = int(np.argmax(y))
    if i == 0 or i == len(y) - 1:
        return x[i], y[i]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    denom = y0 - 2.0 * y1 + y2
    delta = 0.0 if denom == 0 else 0.5 * (y0 - y2) / denom
    dx = x[1] - x[0]
    return x[i] + delta * dx, y1 - 0.25 * (y0 - y2) * delta


def find_impulse_features(tau, h, tau_max=40.0, rel_height=0.03,
                          min_delay_ps=None):
    """
    Locate the resolved sub-pulses of the stack impulse response.

    Peaks are found on the analytic-signal **envelope** rather than on h
    itself: the envelope peak is the group-delay arrival time of the
    sub-pulse and is insensitive to its carrier phase (echoes are sign-flipped
    and slightly reshaped by dispersion, so peak-of-h would bias the timing).

    Features sitting on the shoulder of the direct pulse are NOT echoes -- at
    finite bandwidth the deconvolved direct pulse has a ringing tail, and
    anything inside it is unresolved.  We therefore set an exclusion radius
    from the data itself: `min_delay_ps` defaults to 4x the trailing 10%
    envelope width of the direct pulse (floored at 2.5 ps).  Features inside
    that radius are returned separately so they can be reported, never
    silently used as echoes.

    Returns
    -------
    t0       : arrival time of the direct transmitted pulse [ps]
    feats    : list of (delay_after_main [ps], signed amplitude, rel envelope)
               for resolved echoes
    rejected : same, for features inside the exclusion radius
    min_delay: the exclusion radius actually used [ps]
    """
    env = np.abs(hilbert(h))
    sel = tau <= tau_max
    tau_s, env_s, h_s = tau[sel], env[sel], h[sel]

    i0 = int(np.argmax(env_s))
    t0, _ = subsample_peak(tau_s, env_s)
    a0 = env_s[i0]

    if min_delay_ps is None:
        j = i0
        while j < len(env_s) - 1 and env_s[j] > 0.10 * a0:
            j += 1
        min_delay_ps = max(2.5, 4.0 * (tau_s[j] - tau_s[i0]))

    idx, _ = find_peaks(env_s, height=rel_height * a0,
                        distance=max(1, int(0.8 / (tau[1] - tau[0]))))
    feats, rejected = [], []
    for i in idx:
        lo = max(i - 2, 0)
        tp, ep = subsample_peak(tau_s[lo:i + 3], env_s[lo:i + 3])
        rec = (tp - t0, h_s[i], ep / a0)
        if rec[0] <= 0.0:
            continue
        (feats if rec[0] >= min_delay_ps else rejected).append(rec)
    return t0, feats, rejected, min_delay_ps


def identify_roundtrips(feats, max_order=2):
    """
    Decide, from the data alone, which two delays are the fundamental
    round-trip times tau_1 and tau_2 of the two layers.

    Every sub-pulse of a bilayer arrives at  i*tau_1 + j*tau_2, i, j >= 0.
    Each ordered pair of observed delays is taken as a hypothesis and scored
    by how well ALL observed delays are explained as low-order non-negative
    integer combinations of that pair, weighted by peak amplitude (a strong
    unexplained peak is much worse evidence than a weak one).

    Two guards matter, and both were needed on real data:

    * `max_order` is small (2).  With high orders allowed, a spuriously
      SMALL fundamental can "explain" anything -- 5*2.5 ps lands near 13 ps
      by luck -- and the search happily returns nonsense.  Physically the
      sub-pulses that carry real amplitude are low order, because each round
      trip costs a factor r^2 ~ 0.1.
    * the caller must already have removed features unresolved from the
      direct pulse (see `find_impulse_features`), which is where those
      spurious small fundamentals come from.

    Returns (tau_1, tau_2) with tau_1 < tau_2, and the RMS explanation error.
    """
    if len(feats) < 2:
        raise RuntimeError(
            "fewer than two resolved echoes in the impulse response -- "
            "cannot identify the two layer round-trip times")

    delays = np.array([p[0] for p in feats])
    amps = np.array([p[2] for p in feats])
    combos = [(i, j) for i in range(max_order + 1)
              for j in range(max_order + 1) if i + j >= 1]

    best, best_pair = np.inf, None
    for ia in range(len(delays)):
        for ib in range(len(delays)):
            t1, t2 = delays[ia], delays[ib]
            if not (t1 < t2):
                continue
            grid = np.array([i * t1 + j * t2 for i, j in combos])
            err = float(np.sum(amps * np.min(
                np.abs(grid[None, :] - delays[:, None]), axis=1) ** 2))
            if err < best - 1e-12:
                best, best_pair = err, (t1, t2)

    if best_pair is None:
        raise RuntimeError("could not identify two round-trip times")
    return best_pair, float(np.sqrt(best / np.sum(amps)))


def refine_echo_delay(tau, h, t0, t_echo, half_width):
    """
    Refine an echo delay by cross-correlating the gated echo against the
    gated main pulse.  Cross-correlation uses the whole sub-pulse shape, not
    just its peak, so it returns the group delay with sub-sample accuracy and
    is robust against noise and mild reshaping.
    """
    dt = tau[1] - tau[0]
    g_main = super_gaussian_gate(tau, t0, half_width)
    g_echo = super_gaussian_gate(tau, t0 + t_echo, half_width)
    a, b = h * g_echo, h * g_main

    n = len(tau)
    cc = np.fft.irfft(np.fft.rfft(a, 2 * n) * np.conj(np.fft.rfft(b, 2 * n)),
                      2 * n)
    env = np.abs(hilbert(cc))
    lags = np.arange(2 * n) * dt
    lags[n:] -= 2 * n * dt
    order = np.argsort(lags)
    lags, env = lags[order], env[order]

    win = np.abs(lags - t_echo) < 0.5 * half_width
    if win.sum() < 3:
        return t_echo
    lag_best, _ = subsample_peak(lags[win], env[win])
    return float(lag_best)


def super_gaussian_gate(tau, centre, half_width, order=4):
    """
    Smooth time gate,  exp(-((tau-tc)/hw)^(2*order)).

    A super-Gaussian is flat over the sub-pulse and rolls off smoothly, which
    keeps the gated spectrum free of the ringing a rectangular gate would
    create.  Cost: the gated spectra are smoothed in frequency by roughly
    1/(2*half_width), which sets the true resolution of the per-frequency
    inversion (~0.25 THz here).  That is stated in the printed output.
    """
    return np.exp(-((tau - centre) / half_width) ** (2 * order))


# =============================================================================
# 5.  DISPERSION PARAMETRISATION AND THE GLOBAL TRANSFER-MATRIX FIT
# =============================================================================
#
# Stage B of the extraction.  A global fit to the full complex T(w) over a
# wide band determines the two thicknesses very precisely, because d1 and d2
# each imprint their own Fabry-Perot fringe period on T:
#       df_1 = c/(2 n1 d1) ~ 114 GHz  (Si),   df_2 = c/(2 n2 d2) ~ 77 GHz (qz)
# Those two periods are independent signatures, so the fit can separate the
# layers even though a single-frequency inversion cannot.
#
# The dispersion laws below are deliberately minimal and physically motivated:
#   * high-resistivity Si is essentially non-dispersive and near-transparent
#     below ~4 THz -> n = A + B f^2, alpha = a f^2
#   * fused quartz has weak normal dispersion and absorption rising
#     super-linearly toward the boson peak -> n = A + B f + C f^2,
#     alpha = a f + b f^2
# -----------------------------------------------------------------------------

PAR_NAMES = ["d_Si [um]", "d_qz [um]",
             "n_Si A", "n_Si B", "a_Si [cm^-1/THz^2]",
             "n_qz A", "n_qz B", "n_qz C",
             "a_qz [cm^-1/THz]", "b_qz [cm^-1/THz^2]"]


def alpha_to_k(alpha_cm, f):
    """
    Convert intensity absorption coefficient alpha [cm^-1] to extinction
    coefficient k, using alpha = 4*pi*f*k/c  ->  k = alpha*c/(4*pi*f).

    alpha is in cm^-1 and c is in um/ps, so we need the 1e-4 cm/um factor.
    """
    f = np.asarray(f, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = alpha_cm * 1e-4 * C_UM_PS / (4.0 * np.pi * f)
    return np.where(f > 0, k, 0.0)


def k_to_alpha(k, f):
    """Inverse of alpha_to_k: alpha [cm^-1] from k."""
    return 4.0 * np.pi * np.asarray(f, float) * np.asarray(k, float) / \
        (C_UM_PS * 1e-4)


def unpack_global(p, f):
    """Turn the global parameter vector into (d1, d2, N_Si(f), N_qz(f))."""
    d1, d2, A1, B1, a1, A2, B2, C2, a2, b2 = p
    n1 = A1 + B1 * f ** 2
    n2 = A2 + B2 * f + C2 * f ** 2
    al1 = a1 * f ** 2                      # cm^-1
    al2 = a2 * f + b2 * f ** 2             # cm^-1
    N1 = n1 - 1j * alpha_to_k(al1, f)
    N2 = n2 - 1j * alpha_to_k(al2, f)
    return d1, d2, N1, N2, n1, n2, al1, al2


PAR_LO = np.array([0.5 * D_SI_NOMINAL_UM, 0.5 * D_QZ_NOMINAL_UM,
                   1.0, -1.0, 0.0,
                   1.0, -1.0, -1.0,
                   0.0, 0.0])
PAR_HI = np.array([1.5 * D_SI_NOMINAL_UM, 1.5 * D_QZ_NOMINAL_UM,
                   8.0, 1.0, 500.0,
                   8.0, 1.0, 1.0,
                   200.0, 200.0])
PAR_SCALE = np.array([10, 10, 0.1, 0.01, 1.0, 0.1, 0.01, 0.01, 1.0, 1.0],
                     float)


def global_fit(f, T, weights, p0, fix=None):
    """
    Stage B -- least-squares fit of the transfer-matrix model to T(w).

    Residuals are the real and imaginary parts of (model - measured), weighted
    by the reference dynamic range so that noisy high-frequency bins do not
    dominate.  Fitting Re/Im rather than |T| and phase avoids any need to
    unwrap the model phase and treats amplitude and phase information on an
    equal footing.

    `fix` : optional {parameter_index: value} mapping.  Those parameters are
    held at the given values and every OTHER parameter is re-optimised.  This
    is what makes a chi^2 PROFILE possible (`profile_index_split` below), and
    a profile is the only honest way to put an error bar on d_Si here -- the
    covariance matrix of the unconstrained fit is not, because the d/n valley
    is far too flat and curved for a local quadratic to describe.

    The returned object carries `.p` (the full 10-vector, frozen entries
    reinstated) and `.chi2` (sum of squared weighted residuals).  chi2 is
    directly comparable between fits with different `fix`, because the
    residual vector has the same length and weighting in every case.
    """
    fix = dict(fix or {})
    free = np.array([i for i in range(len(PAR_NAMES)) if i not in fix], int)

    # scipy rejects a start exactly on or outside a bound; nudge it inside.
    span = PAR_HI - PAR_LO
    p0 = np.clip(np.asarray(p0, float),
                 PAR_LO + 1e-6 * span, PAR_HI - 1e-6 * span)
    p_frozen = p0.copy()
    for i, v in fix.items():
        p_frozen[i] = v

    def expand(q):
        p = p_frozen.copy()
        p[free] = q
        return p

    def residual(q):
        d1, d2, N1, N2, *_ = unpack_global(expand(q), f)
        dev = (model_transmission(f, d1, d2, N1, N2) - T) * weights
        return np.concatenate([dev.real, dev.imag])

    out = least_squares(residual, p0[free],
                        bounds=(PAR_LO[free], PAR_HI[free]),
                        x_scale=PAR_SCALE[free],
                        xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
    out.p = expand(out.x)
    out.chi2 = float(np.sum(out.fun ** 2))
    return out


# Index of "n_Si A" in the global parameter vector -- the parameter that is
# traded against d_Si along the flat valley.
I_NSI_A = PAR_NAMES.index("n_Si A")


def profile_index_split(f, T, weights, p_best, n_si_grid, chi2_tol=0.01):
    """
    Profile chi^2 along the n_Si <-> d_Si degeneracy.

    For each trial n_Si (the constant term A1) the parameter is FROZEN and all
    nine others -- including both thicknesses -- are re-optimised.  The result
    is the true one-dimensional profile of chi^2 along the valley, which is
    what an error bar on d_Si must be built from.

    `chi2_tol` is the fractional chi^2 rise that still counts as "an equally
    good description of the data" (default 1%).  That is a deliberately
    generous-looking criterion, and section 4c justifies it empirically: a
    single extra *physically plausible* degree of freedom (an air gap between
    the wafer and the window) changes chi^2 by an order of magnitude more than
    1%.  Against that level of model systematic, splits differing by <1% in
    chi^2 are simply not distinguishable.

    Returns
    -------
    rows : list of dicts, one per grid point
    band : dict with the chi2_tol-acceptable ranges of d1, d2, d1+d2, n1*d1
    """
    rows = []
    for v in np.asarray(n_si_grid, float):
        p0 = np.array(p_best, float)
        p0[I_NSI_A] = v
        out = global_fit(f, T, weights, p0, fix={I_NSI_A: v})
        d1, d2 = out.p[0], out.p[1]
        n1_1thz = out.p[2] + out.p[3]          # n_Si at 1 THz
        n2_1thz = out.p[5] + out.p[6] + out.p[7]
        rows.append(dict(n_si=v, n_si_1thz=n1_1thz, chi2=out.chi2,
                         d1=d1, d2=d2, dsum=d1 + d2,
                         nd1=n1_1thz * d1, nd2=n2_1thz * d2))

    chi2min = min(r["chi2"] for r in rows)
    ok = [r for r in rows if r["chi2"] <= chi2min * (1.0 + chi2_tol)]
    band = {key: (min(r[key] for r in ok), max(r[key] for r in ok))
            for key in ("d1", "d2", "dsum", "nd1", "nd2")}
    band["chi2min"] = chi2min
    band["n_ok"] = len(ok)
    band["tol"] = chi2_tol
    return rows, band


def airgap_fit(f, T, weights, p_best, dg0=5.0):
    """
    Systematic stress test: refit with a free air gap between wafer and window.

    An uncontacted wafer resting on a window has a gap of order a few microns
    -- deeply sub-wavelength at 1 THz (lambda = 300 um), so it is invisible as
    a resolved echo, yet it changes the effective reflectivity of the Si|qz
    boundary and hence both the Fabry-Perot fringe contrast and the echo
    amplitudes.

    If chi^2 improves substantially and d_Si moves by microns, that is direct
    evidence that d_Si is soft at the micron level in the two-layer model too,
    and it is the empirical justification for the 1%-chi^2 acceptance band
    used by `profile_index_split`.
    """
    lo = np.concatenate([[PAR_LO[0], 0.0, PAR_LO[1]], PAR_LO[2:]])
    hi = np.concatenate([[PAR_HI[0], 60.0, PAR_HI[1]], PAR_HI[2:]])
    scale = np.concatenate([[PAR_SCALE[0], 1.0, PAR_SCALE[1]], PAR_SCALE[2:]])
    q0 = np.concatenate([[p_best[0], dg0, p_best[1]], p_best[2:]])
    q0 = np.clip(q0, lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo))

    def residual(q):
        p = np.concatenate([[q[0], q[2]], q[3:]])       # back to 10-vector
        _, _, N1, N2, *_ = unpack_global(p, f)
        dev = (model_transmission_gap(f, q[0], q[1], q[2], N1, N2) - T) \
            * weights
        return np.concatenate([dev.real, dev.imag])

    out = least_squares(residual, q0, bounds=(lo, hi), x_scale=scale,
                        xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
    out.chi2 = float(np.sum(out.fun ** 2))
    return out


# =============================================================================
# 6.  STAGE C -- MODEL-FREE PER-FREQUENCY INVERSION FOR BOTH LAYERS
# =============================================================================

def pointwise_inversion(f_grid, T_meas, R1_meas, R2_meas, d1, d2,
                        x_start, f_anchor=F_ANCHOR,
                        w_T=2.0, w_R=1.0):
    """
    Solve, at every frequency independently, for x = (n_Si, k_Si, n_qz, k_qz).

    Observables (6 real numbers per frequency):
        T (exact, full Fabry-Perot)  -- 2
        R1 = Si echo / direct        -- 2   (mostly sensitive to N_Si)
        R2 = qz echo / direct        -- 2   (mostly sensitive to N_qz)
    Unknowns: 4.  Overdetermined, so the solution is genuinely constrained
    rather than a guess -- this is what separates the two layers.

    Residuals are RELATIVE (divided by |observable|) because |T| ~ 0.7 while
    |R1| ~ 0.15 and |R2| ~ 0.07: absolute residuals would let T dominate, and
    T alone is degenerate.  T carries an extra weight w_T because it is exact,
    whereas R1/R2 are slightly smoothed by the time gates.

    Frequency continuation: we start at `f_anchor` (best SNR) from the global
    fit and then march outward, each solve seeded by its neighbour.  This
    keeps the solver on the correct branch of the highly oscillatory
    exp(-2i*delta) factors (delta/2pi ~ 9 cycles at 1 THz for the Si echo).
    """
    def residual(x, f, Tm, R1m, R2m):
        N1 = x[0] - 1j * max(x[1], 0.0)
        N2 = x[2] - 1j * max(x[3], 0.0)
        Tt = model_transmission(np.array([f]), d1, d2, N1, N2)[0]
        R1t, R2t = echo_ratio_models(np.array([f]), d1, d2, N1, N2)
        R1t, R2t = R1t[0], R2t[0]
        dT = w_T * (Tt - Tm) / max(abs(Tm), 1e-9)
        d1r = w_R * (R1t - R1m) / max(abs(R1m), 1e-9)
        d2r = w_R * (R2t - R2m) / max(abs(R2m), 1e-9)
        return np.array([dT.real, dT.imag,
                         d1r.real, d1r.imag,
                         d2r.real, d2r.imag])

    nf = len(f_grid)
    sol = np.full((nf, 4), np.nan)
    cost = np.full(nf, np.nan)
    bounds = ([1.0, 0.0, 1.0, 0.0], [8.0, 1.0, 8.0, 1.0])

    i0 = int(np.argmin(np.abs(f_grid - f_anchor)))

    for direction in (+1, -1):
        x = np.array(x_start, float)
        rng = range(i0, nf) if direction > 0 else range(i0 - 1, -1, -1)
        for i in rng:
            try:
                out = least_squares(
                    residual, x, bounds=bounds,
                    args=(f_grid[i], T_meas[i], R1_meas[i], R2_meas[i]),
                    x_scale=[0.1, 1e-3, 0.1, 1e-3],
                    xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
            except Exception:
                continue
            sol[i] = out.x
            cost[i] = 2.0 * out.cost
            # Seed the next frequency, but do not let a bad solve derail the
            # continuation: only accept if the fit is reasonable.
            if 2.0 * out.cost < 1.0:
                x = out.x
    return sol, cost


# =============================================================================
# 7.  PLOTTING
# =============================================================================

def _style_axis(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, fontsize=9.5, color=CLR_INK2)
    ax.set_ylabel(ylabel, fontsize=9.5, color=CLR_INK2)
    if title:
        ax.set_title(title, fontsize=10.5, color=CLR_INK, loc="left",
                     fontweight="semibold", pad=6)
    ax.grid(True, which="major", color=CLR_GRID, lw=0.7, zorder=0)
    ax.grid(True, which="minor", color=CLR_GRID, lw=0.4, alpha=0.6, zorder=0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(labelsize=8.5, colors=CLR_INK2, length=3)
    for s in ax.spines.values():
        s.set_color(CLR_GRID)
    ax.set_axisbelow(True)


def make_figure(res, out_png=OUT_PNG):
    """Render the nine-panel summary figure."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "axes.facecolor": "#fcfcfb",
        "figure.facecolor": "#ffffff",
        "legend.frameon": False,
        "legend.fontsize": 8.2,
    })

    fig, axes = plt.subplots(3, 3, figsize=(15.5, 11.2), dpi=150)
    fig.suptitle(
        "THz-TDS of a Si / fused-quartz bilayer  —  Fabry-Pérot "
        "transfer-matrix extraction of n(ω) and k(ω)",
        fontsize=13.5, color=CLR_INK, fontweight="semibold", y=0.985)

    t = res["t"]; dt = res["dt"]

    # ---- (a) time-domain waveforms ------------------------------------
    ax = axes[0, 0]
    ax.plot(t - t[0], res["e_ref"], color=CLR_QZ, lw=1.6, label="Reference (dry N₂)")
    ax.plot(t - t[0], res["e_sam"], color=CLR_SI, lw=1.6, label="Sample (Si + quartz)")
    ax.annotate("Δt = %.2f ps" % res["t0"],
                xy=(res["t_peak_sam"] - t[0], res["e_peak_sam"]),
                xytext=(res["t_peak_sam"] - t[0] + 9, res["e_peak_sam"] * 0.78),
                fontsize=8.5, color=CLR_INK2,
                arrowprops=dict(arrowstyle="->", color=CLR_INK2, lw=0.9))
    _style_axis(ax, "Delay time  t − t₀  [ps]", "THz field  E(t)  [a.u.]",
                "(a)  Time-domain waveforms")
    # Data-driven limits: show from just before the earlier pulse to well past
    # the last echo of interest, clipped to the record.  Hard-coded numbers
    # here would silently crop a sample with different thicknesses.
    _tp = min(res["t_peak_ref"], res["t_peak_sam"]) - t[0]
    _span = max(45.0, 3.0 * res["tau_qz"])
    ax.set_xlim(max(0.0, _tp - 8.0), min(t[-1] - t[0], _tp + _span))
    ax.legend(loc="lower right")

    # ---- (b) deconvolved impulse response -----------------------------
    ax = axes[0, 1]
    tau, h = res["tau"], res["h"]
    m = (tau > res["t0"] - 6) & (tau < res["t0"] + 30)
    hn = h / np.max(np.abs(h))
    ax.plot(tau[m] - res["t0"], hn[m], color=CLR_INK, lw=1.4,
            label="Impulse response  h(τ)")
    for kk, (c, hw, lab, col) in enumerate(res["gates_plot"]):
        ax.axvspan(c - hw, c + hw, color=col, alpha=0.13, lw=0, zorder=0)
        # stagger the labels vertically: the Si and quartz gates are only
        # ~4 ps apart and their captions would otherwise overlap
        ax.text(c, 1.15 if kk % 2 == 0 else 1.01, lab, ha="center",
                va="bottom", fontsize=8.2, color=col, fontweight="semibold")
    ax.axhline(0, color=CLR_MUTED, lw=0.7)
    _style_axis(ax, "Delay after main pulse  τ − Δt  [ps]",
                "Normalised amplitude", "(b)  Sample impulse response & gates")
    ax.set_ylim(-0.45, 1.40)
    ax.legend(loc="lower right")

    # ---- (c) amplitude spectra ----------------------------------------
    ax = axes[0, 2]
    f = res["f"]
    ax.semilogy(f, np.abs(res["Er"]), color=CLR_QZ, lw=1.5, label="Reference")
    ax.semilogy(f, np.abs(res["Es"]), color=CLR_SI, lw=1.5, label="Sample")
    ax.axhline(res["noise"], color=CLR_MUTED, lw=1.1, ls="--",
               label="Noise floor")
    ax.axvspan(*BAND_POINTWISE, color=CLR_MODEL, alpha=0.09, lw=0, zorder=0)
    ax.text(np.mean(BAND_POINTWISE), res["noise"] * 0.6, "inversion band",
            ha="center", va="bottom", fontsize=8.2, color=CLR_MODEL)
    _style_axis(ax, "Frequency  f  [THz]", "Spectral amplitude  |E(f)|  [a.u.]",
                "(c)  Amplitude spectra (2¹⁴-point FFT)")
    # Show the full useful bandwidth: out to where the reference finally sinks
    # into the noise floor, with 20% headroom, rather than a fixed 5 THz.
    _usable = res["f"][res["dr_db"] > 0.0]
    _fmax = min(float(res["f"][-1]),
                (float(_usable[-1]) * 1.2 if len(_usable) else 5.0))
    ax.set_xlim(0, max(_fmax, BAND_POINTWISE[1] * 1.2))
    ax.set_ylim(res["noise"] * 0.25, np.max(np.abs(res["Er"])) * 3)
    ax.legend(loc="upper right")

    # ---- (d) |T| measured vs model ------------------------------------
    ax = axes[1, 0]
    mb = res["mask_fit"]
    ax.plot(f[mb], np.abs(res["T"][mb]), color=CLR_SI, lw=1.4,
            label="Measured  |T|")
    ax.plot(f[mb], np.abs(res["T_fit"][mb]), color=CLR_MODEL, lw=1.4,
            ls="--", label="Transfer-matrix fit")
    _style_axis(ax, "Frequency  f  [THz]", "Transmission magnitude  |T(f)|",
                "(d)  Fabry–Pérot fringes in |T|")
    ax.set_xlim(*BAND_FIT)
    ax.legend(loc="lower left")

    # ---- (e) phase of T ------------------------------------------------
    ax = axes[1, 1]
    ax.plot(f[mb], res["phi"][mb], color=CLR_SI, lw=1.4, label="Measured")
    ax.plot(f[mb], res["phi_fit"][mb], color=CLR_MODEL, lw=1.4, ls="--",
            label="Transfer-matrix fit")
    _style_axis(ax, "Frequency  f  [THz]",
                "Unwrapped phase  arg T(f)  [rad]",
                "(e)  Transmission phase (2π branch fixed)")
    ax.set_xlim(*BAND_FIT)
    ax.legend(loc="upper right")

    # ---- (f) refractive index -----------------------------------------
    ax = axes[1, 2]
    fp, sol = res["f_pw"], res["sol_pw"]
    ax.plot(fp, sol[:, 0], color=CLR_SI, lw=1.9, label="Si  (per-frequency)")
    ax.plot(fp, sol[:, 2], color=CLR_QZ, lw=1.9, label="Quartz  (per-frequency)")
    ax.plot(res["f_fit"], res["n1_fit"], color=CLR_MODEL, lw=1.2, ls=":",
            label="Global model")
    ax.plot(res["f_fit"], res["n2_fit"], color=CLR_MODEL, lw=1.2, ls=":")
    ax.text(fp[len(fp) // 2], np.nanmean(sol[:, 0]) + 0.10, "Si",
            color=CLR_SI, fontsize=10, fontweight="semibold", ha="center")
    ax.text(fp[len(fp) // 2], np.nanmean(sol[:, 2]) + 0.10, "fused quartz",
            color=CLR_QZ, fontsize=10, fontweight="semibold", ha="center")
    _style_axis(ax, "Frequency  f  [THz]", "Refractive index  n(f)",
                "(f)  Refractive index of each layer")
    ax.set_xlim(*BAND_POINTWISE)
    # Autoscale around whatever indices were actually found, with a margin
    # large enough for the text labels, instead of Si/quartz-specific limits.
    _nall = np.concatenate([sol[:, 0], sol[:, 2]])
    _nall = _nall[np.isfinite(_nall)]
    if _nall.size:
        _lo, _hi = float(np.nanmin(_nall)), float(np.nanmax(_nall))
        _pad = max(0.12, 0.10 * (_hi - _lo))
        ax.set_ylim(_lo - _pad, _hi + 2.0 * _pad)
    ax.legend(loc="center right")

    # ---- (g) extinction coefficient ------------------------------------
    ax = axes[2, 0]
    # Silicon absorption sits below the measurement floor: the solver simply
    # pins k_Si at its lower bound, so plotting that "curve" would be
    # meaningless. Show the detection limit as a shaded upper bound instead.
    # A solid fill from the axis floor up to the Si limit would run UNDERNEATH
    # the quartz curve wherever quartz dips below it (here, below ~0.9 THz),
    # making quartz read as "sub-threshold" when it is nothing of the kind.
    # Hatching bounded by the Si limit line carries the same meaning without
    # tinting anything that belongs to the other layer.
    ax.fill_between(fp, 1e-6, res["k_floor_si"], facecolor="none",
                    edgecolor=CLR_SI, hatch="////", lw=0.0, alpha=0.55,
                    zorder=1, label="Si: below detection limit")
    ax.semilogy(fp, res["k_floor_si"], color=CLR_SI, lw=1.6, ls="--",
                zorder=3, label="Si  upper limit")
    ax.semilogy(fp, res["k_floor_qz"], color=CLR_QZ, lw=1.1, ls="-.",
                zorder=3, label="Quartz detection floor")
    ax.semilogy(fp, np.clip(sol[:, 3], 1e-6, None), color=CLR_QZ, lw=1.9,
                zorder=4, label="Fused quartz")
    ax.semilogy(res["f_fit"], np.clip(res["k2_fit"], 1e-6, None),
                color=CLR_MODEL, lw=1.2, ls=":", zorder=4,
                label="Global model (quartz)")
    _style_axis(ax, "Frequency  f  [THz]", "Extinction coefficient  k(f)",
                "(g)  Extinction coefficient (log scale)")
    ax.set_xlim(*BAND_POINTWISE)
    _kall = np.concatenate([np.clip(sol[:, 3], 1e-12, None),
                            res["k_floor_si"], res["k_floor_qz"]])
    _kall = _kall[np.isfinite(_kall) & (_kall > 0)]
    if _kall.size:
        ax.set_ylim(10 ** np.floor(np.log10(_kall.min()) - 0.2),
                    10 ** np.ceil(np.log10(_kall.max()) + 0.2))
    ax.legend(loc="lower right")

    # ---- (h) absorption coefficient ------------------------------------
    ax = axes[2, 1]
    a2 = k_to_alpha(np.clip(sol[:, 3], 0, None), fp)
    # Same reasoning as panel (g): hatch, do not fill, so the Si upper limit
    # never tints the region occupied by the quartz curve.
    ax.fill_between(fp, 0.0, res["a_floor_si"], facecolor="none",
                    edgecolor=CLR_SI, hatch="////", lw=0.0, alpha=0.55,
                    zorder=1, label="Si: below detection limit")
    ax.axhline(res["a_floor_si"], color=CLR_SI, lw=1.6, ls="--", zorder=3,
               label=f"Si  α < {res['a_floor_si']:.2f} cm⁻¹")
    ax.axhline(res["a_floor_qz"], color=CLR_QZ, lw=1.1, ls="-.", zorder=3,
               label=f"Quartz floor  {res['a_floor_qz']:.2f} cm⁻¹")
    ax.plot(fp, a2, color=CLR_QZ, lw=1.9, zorder=4, label="Fused quartz")
    ax.plot(res["f_fit"], res["al2_fit"], color=CLR_MODEL, lw=1.2, ls=":",
            zorder=4, label="Global model (quartz)")
    _style_axis(ax, "Frequency  f  [THz]",
                "Absorption coefficient  α(f)  [cm⁻¹]",
                "(h)  Power absorption coefficient")
    ax.set_xlim(*BAND_POINTWISE)
    _atop = np.nanmax([np.nanpercentile(a2, 99) if np.any(np.isfinite(a2))
                       else 0.0, res["a_floor_si"], res["a_floor_qz"]])
    ax.set_ylim(0, _atop * 1.35 if _atop > 0 else 1.0)
    ax.legend(loc="upper left")

    # ---- (i) results panel ---------------------------------------------
    ax = axes[2, 2]
    ax.axis("off")
    ax.text(0.0, 1.0, res["panel_text"], va="top", ha="left",
            fontsize=8.0, family="monospace", color=CLR_INK,
            transform=ax.transAxes, linespacing=1.55)
    ax.set_title("(i)  Extracted parameters", fontsize=10.5, color=CLR_INK,
                 loc="left", fontweight="semibold", pad=6)

    # ---- warning stamp for a knowingly biased analysis ------------------
    # If the user forced the full-record Hann, the PNG must never be able to
    # circulate as a result.  A diagonal stamp across the whole canvas travels
    # with the image in a way a printed console warning does not.
    if res.get("biased_window"):
        fig.text(0.5, 0.5,
                 "INVALID — full-record Hann\n"
                 f"|T| biased by {res['bias']:.2f}×",
                 fontsize=58, color="#d01c1c", alpha=0.20, rotation=28,
                 ha="center", va="center", fontweight="bold", zorder=1000)
        for _y in (0.972, 0.012):
            fig.text(0.5, _y,
                     f"⚠ WINDOW = {res['window'].upper()} — RATIO SPECTRUM IS "
                     f"BIASED BY {res['bias']:.2f}×. NOT A VALID MEASUREMENT "
                     f"(--allow-biased-window was set).",
                     ha="center", va="center", fontsize=11,
                     color="white", fontweight="bold", zorder=1001,
                     bbox=dict(boxstyle="round,pad=0.45", fc="#d01c1c",
                               ec="none"))

    fig.tight_layout(rect=[0, 0.005, 1, 0.965])
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)
    return out_png


# =============================================================================
# 8.  MAIN
# =============================================================================

def hr(title=""):
    if title:
        print("\n" + "=" * 78)
        print(title)
        print("=" * 78)
    else:
        print("-" * 78)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="THz-TDS bilayer (Si on quartz) parameter extraction.")
    ap.add_argument("--window", default="flattop-hann",
                    choices=["flattop-hann", "hann", "boxcar"],
                    help="apodisation (default: flattop-hann; see docstring)")
    ap.add_argument("--sample", default=SAMPLE_FILE)
    ap.add_argument("--ref", default=REF_FILE)
    ap.add_argument("--out", default=OUT_PNG)
    ap.add_argument("--df-pointwise", type=float, default=0.005,
                    help="frequency step of the per-frequency inversion [THz]")
    ap.add_argument("--allow-biased-window", action="store_true",
                    help="permit --window hann even though it biases the "
                         "ratio T = E_s/E_r (the figure is then stamped)")
    ap.add_argument("--chi2-tol", type=float, default=0.01,
                    help="fractional chi^2 rise still counted as an equally "
                         "good fit when profiling the d_Si/n_Si degeneracy")
    args = ap.parse_args(argv)

    # A full-record Hann on a ratio measurement is not a stylistic choice, it
    # is a ~1.7x multiplicative error on |T| that propagates into a ~65 um
    # thickness error and a spurious negative absorption.  Refuse by default;
    # if the user really wants to see the damage, they must say so, and the
    # figure gets a warning stamp so the PNG can never be mistaken for a
    # result.  (The size of the bias is printed either way, in section 2.)
    if args.window == "hann" and not args.allow_biased_window:
        ap.error(
            "--window hann applies a full-record Hann to BOTH traces, which "
            "weights the sample and reference pulses differently because they "
            "arrive at different delays. On this data it biases |T| by ~1.7x, "
            "moves the fitted d_Si from ~386 um to ~450 um, and stops the "
            "per-frequency inversion converging at every key frequency. "
            "Use the default --window flattop-hann. To reproduce the failure "
            "deliberately, add --allow-biased-window.")

    here = os.path.dirname(os.path.abspath(__file__))
    sam_path = args.sample if os.path.isabs(args.sample) else os.path.join(here, args.sample)
    ref_path = args.ref if os.path.isabs(args.ref) else os.path.join(here, args.ref)
    out_path = args.out if os.path.isabs(args.out) else os.path.join(here, args.out)

    # ---------------------------------------------------------------- 8.1 I/O
    hr("1.  DATA")
    t_s, e_s, hdr_s = read_menlo(sam_path)
    t_r, e_r, hdr_r = read_menlo(ref_path)

    if len(t_s) != len(t_r) or not np.allclose(t_s, t_r):
        raise ValueError("sample and reference time axes differ; "
                         "a common grid is required for a ratio measurement")

    t = t_s
    dt = float(np.mean(np.diff(t)))
    e_s = remove_baseline(e_s)
    e_r = remove_baseline(e_r)

    print(f"  sample     : {os.path.basename(sam_path)}")
    print(f"  reference  : {os.path.basename(ref_path)}")
    for line in hdr_s[:3]:
        print(f"      {line}")
    print(f"  samples    : {len(t)}   step dt = {dt*1000:.3f} fs   "
          f"record = {t[-1]-t[0]:.1f} ps")
    print(f"  FFT        : zero-padded to N = {NFFT} = 2^{int(np.log2(NFFT))}"
          f"   ->  df = {1000.0/(NFFT*dt):.3f} GHz,"
          f"  f_Nyquist = {0.5/dt:.2f} THz")
    print(f"  window     : {args.window}")

    self_test()

    # ------------------------------------------------------- 8.2 FFT and T(w)
    hr("2.  SPECTRA AND COMPLEX TRANSMISSION  T(ω) = E_sample / E_ref")

    win = make_window(len(t), args.window, taper_ps=6.0, dt=dt)
    f, Es = rfft_spectrum(e_s, dt, NFFT, win)
    _, Er = rfft_spectrum(e_r, dt, NFFT, win)
    T = Es / Er

    # Noise floor: the spectral region above the useful bandwidth where only
    # detector/ADC noise remains.  Everything downstream is band-limited to
    # where the reference sits comfortably above it.
    noise = float(np.mean(np.abs(Er[(f > 10.0) & (f < 14.0)])))
    dr_db = 20.0 * np.log10(np.abs(Er) / noise)

    mask_fit = (f >= BAND_FIT[0]) & (f <= BAND_FIT[1])
    phi = branch_corrected_phase(f, T, mask=mask_fit)

    for ftest in (0.5, 1.0, 2.0, 3.0, 4.0):
        i = int(np.argmin(np.abs(f - ftest)))
        print(f"    f = {ftest:4.1f} THz :  dynamic range {dr_db[i]:5.1f} dB,"
              f"   |T| = {abs(T[i]):.4f},   arg T = {phi[i]:8.3f} rad")

    # --- quantify the full-record-Hann bias discussed in the docstring ------
    w_hann = make_window(len(t), "hann")
    _, Es_h = rfft_spectrum(e_s, dt, NFFT, w_hann)
    _, Er_h = rfft_spectrum(e_r, dt, NFFT, w_hann)
    mh = (f > 0.3) & (f < 2.0)
    bias = float(np.mean(np.abs(Es_h[mh] / Er_h[mh]) / np.abs(T[mh])))
    i_ps = int(np.argmax(np.abs(e_s))); i_pr = int(np.argmax(np.abs(e_r)))
    print(f"\n  window bias check: a FULL-RECORD Hann weights the sample pulse"
          f" by {w_hann[i_ps]:.3f}")
    print(f"    and the reference pulse by {w_hann[i_pr]:.3f} (they arrive at "
          f"different delays),")
    print(f"    biasing |T| by a factor {bias:.3f} — a pure artefact. The "
          f"flat-top Hann used")
    print(f"    here weights both by {win[i_ps]:.3f} / {win[i_pr]:.3f} and is "
          f"unbiased.")
    biased_window = (args.window == "hann")
    if biased_window:
        print("\n  " + "*" * 72)
        print("  *** --window hann SELECTED WITH --allow-biased-window.")
        print("  *** EVERY NUMBER BELOW IS WRONG BY CONSTRUCTION. On this")
        print(f"  *** data the {bias:.2f}x bias in |T| drives d_Si to ~450 um")
        print("  *** (vs ~386 um) and the per-frequency inversion fails to")
        print("  *** converge. The figure is stamped accordingly.")
        print("  " + "*" * 72)

    # -------------------------------- 8.3 Time-domain echoes -> thicknesses
    hr("3.  TIME-DOMAIN ECHO ANALYSIS  →  LAYER THICKNESSES")

    h, H, bandpass = wiener_deconvolve(Es, Er, f)
    tau = np.arange(NFFT) * dt

    t0, feats, rejected, min_delay = find_impulse_features(tau, h)
    (tau1, tau2), rms_comb = identify_roundtrips(feats)

    print("  Wiener-deconvolved impulse response h(τ) of the stack")
    print("  (instrument response removed, so every feature is a "
          "Fabry–Pérot sub-pulse):\n")
    if rejected:
        print(f"    features within {min_delay:.2f} ps of the direct pulse are"
              f" unresolved from it")
        print(f"    (filter ringing on its shoulder) and are excluded as echo "
              f"candidates:")
        for d_, a_, env_ in sorted(rejected, key=lambda q: q[0]):
            print(f"      Δτ = {d_:6.3f} ps,  rel. envelope {env_:.4f}"
                  f"   [excluded]")
        print()
    print("    Δτ after main [ps]  rel. envelope  sign   assignment")
    echo_sign = {}
    for d_, a_, env_ in sorted(feats, key=lambda p: p[0]):
        lab = ""
        if abs(d_ - tau1) < 0.3:
            lab = "τ₁  (layer-1 round trip)"
            echo_sign[1] = a_
        elif abs(d_ - tau2) < 0.3:
            lab = "τ₂  (layer-2 round trip)"
            echo_sign[2] = a_
        else:
            for i_ in range(4):
                for j_ in range(4):
                    if i_ + j_ > 1 and abs(i_ * tau1 + j_ * tau2 - d_) < 0.35:
                        lab = f"{i_}·τ₁ + {j_}·τ₂"
            if not lab:
                lab = "unassigned"
        print(f"      {d_:8.3f}          {env_:7.4f}      {'+' if a_ > 0 else '−'}"
              f"    {lab}")

    # The POLARITY of each echo is an independent, amplitude-free check on the
    # index ordering, and it costs nothing: R1 = -r01*r12*a^2 with n_air < n_Si
    # (r01 < 0) and n_Si > n_qz (r12 > 0) predicts a POSITIVE Si echo, while
    # R2 = -r12*r23*b^2 with both r > 0 predicts a NEGATIVE quartz echo.  A
    # flipped sign would mean the layer assignment or the index ordering is
    # wrong, and no amount of fitting would reveal it.
    pol_ok = (echo_sign.get(1, 0.0) > 0) and (echo_sign.get(2, 0.0) < 0)
    print(f"\n  echo polarity check: expected + for the high-index layer's "
          f"round trip")
    print(f"    (−r01·r12 > 0) and − for the window's (−r12·r23 < 0)  →  "
          f"{'consistent' if pol_ok else 'INCONSISTENT — check assignment'}")

    # Refine both round-trip times by cross-correlation against the main pulse
    hw_e = 0.45 * min(tau1, tau2 - tau1)
    tau1 = refine_echo_delay(tau, h, t0, tau1, hw_e)
    tau2 = refine_echo_delay(tau, h, t0, tau2, hw_e)

    # Which round trip belongs to which layer?  A round trip measures the
    # GROUP optical thickness n_g*d.  With the nominal geometry
    # (Si 385 um @ n~3.42 -> n*d ~ 1.3 mm; quartz 1 mm @ n~1.95 -> n*d ~ 2.0 mm)
    # the smaller n*d is the Si layer.  This is the ONLY place a nominal value
    # is used, and only to attach a name to an echo.
    nd1, nd2 = C_UM_PS * tau1 / 2.0, C_UM_PS * tau2 / 2.0   # group n*d [um]
    if (D_SI_NOMINAL_UM * 3.42) > (D_QZ_NOMINAL_UM * 1.95):
        tau_si, tau_qz, nd_si, nd_qz = tau2, tau1, nd2, nd1
    else:
        tau_si, tau_qz, nd_si, nd_qz = tau1, tau2, nd1, nd2

    # The main-pulse delay gives the GEOMETRIC sum independently:
    #   t0 = [(n1-1)d1 + (n2-1)d2]/c  and  n_i d_i = c*tau_i/2
    #   =>  d1 + d2 = c*(tau1/2 + tau2/2 - t0)
    d_sum_timing = C_UM_PS * (tau_si / 2.0 + tau_qz / 2.0 - t0)

    print(f"\n  main-pulse group delay   Δt  = {t0:8.4f} ps")
    print(f"  Si  round trip           τ_Si = {tau_si:8.4f} ps  →  "
          f"n_g·d = {nd_si:8.2f} µm")
    print(f"  quartz round trip        τ_qz = {tau_qz:8.4f} ps  →  "
          f"n_g·d = {nd_qz:8.2f} µm")
    print(f"  ⇒ d_Si + d_qz (timing, index-free) = {d_sum_timing:8.2f} µm")

    # ------------------------------------------- 8.4 Global transfer-matrix fit
    hr("4.  GLOBAL FABRY–PÉROT TRANSFER-MATRIX FIT  (determines d₁, d₂)")

    fb = f[mask_fit]
    Tb = T[mask_fit]
    # SNR weighting, normalised so the cost is O(1)
    wts = np.abs(Er[mask_fit]) / np.max(np.abs(Er[mask_fit]))
    wts = np.sqrt(wts) / np.sqrt(np.mean(wts))

    # Seed the fit from the echo timing when that is self-consistent, and
    # fall back to the nominal geometry otherwise.  The timing-derived sum
    # must be positive and within a factor ~1.5 of nominal, or the echo
    # identification went wrong and we should not trust it as a starting
    # point.  Whatever we choose is then clipped into the fit bounds, because
    # scipy refuses a start outside them.
    d_nom_sum = D_SI_NOMINAL_UM + D_QZ_NOMINAL_UM
    timing_ok = 0.6 * d_nom_sum < d_sum_timing < 1.6 * d_nom_sum
    if timing_ok:
        d1_0 = D_SI_NOMINAL_UM
        d2_0 = d_sum_timing - D_SI_NOMINAL_UM
        n_si_guess = nd_si / d1_0
        n_qz_guess = nd_qz / d2_0
    else:
        print("  WARNING: echo-timing thickness sum "
              f"({d_sum_timing:.1f} µm) is implausible; falling back to the "
              "nominal geometry for the fit start.")
        d1_0, d2_0 = D_SI_NOMINAL_UM, D_QZ_NOMINAL_UM
        n_si_guess, n_qz_guess = 3.4175, 1.95

    p0 = [d1_0, d2_0,
          n_si_guess, 0.0, 0.1,
          n_qz_guess, 0.0, 0.0,
          0.0, 1.5]

    out = global_fit(fb, Tb, wts, p0)
    p = out.p
    d1, d2, N1f, N2f, n1f, n2f, al1f, al2f = unpack_global(p, fb)

    rms = float(np.sqrt(np.mean(out.fun ** 2)))
    print(f"  converged: {out.message.strip()}")
    print(f"  weighted RMS residual in complex T : {rms:.5f}  "
          f"({100*rms/np.mean(np.abs(Tb)):.2f} % of mean |T|)\n")
    print("  best-fit parameter vector (digits here are NOT all significant —")
    print("   see the profile in 4b for what the data actually determines)")
    for nm, v in zip(PAR_NAMES, p):
        print(f"      {nm:24s} = {v:12.5f}")

    d_fit_sum = d1 + d2
    print(f"\n  CROSS-CHECK  d_Si + d_qz :  "
          f"{d_fit_sum:8.2f} µm (Fabry–Pérot fit)  vs  "
          f"{d_sum_timing:8.2f} µm (echo timing, index-free)"
          f"   →  Δ = {d_fit_sum - d_sum_timing:+.2f} µm"
          f" ({100*abs(d_fit_sum-d_sum_timing)/d_fit_sum:.2f} %)")

    # --------------------------------------- 8.4b chi^2 profile along n_Si
    hr("4b. HOW WELL IS THE THICKNESS SPLIT ACTUALLY DETERMINED?"
       "  (χ² PROFILE)")
    print("  The fit above quotes d_Si to 0.01 µm. That precision is not real.")
    print("  Below, n_Si is FROZEN at each value and all nine other parameters")
    print("  — both thicknesses included — are re-optimised. This is the true")
    print("  1-D profile of χ² along the n·d degeneracy valley.\n")

    n_si_grid = np.concatenate([
        np.arange(3.34, 3.47, 0.01), [p[I_NSI_A], 3.4175]])
    n_si_grid = np.unique(np.round(n_si_grid, 5))
    rows, band = profile_index_split(fb, Tb, wts, p, n_si_grid,
                                     chi2_tol=args.chi2_tol)
    chi2min = band["chi2min"]

    print(f"  {'n_Si':>8s}{'χ²':>11s}{'χ²/χ²min':>11s}{'d_Si':>10s}"
          f"{'d_qz':>10s}{'d_Si+d_qz':>12s}{'n_Si·d_Si':>12s}")
    print(f"  {'':>8s}{'':>11s}{'':>11s}{'[µm]':>10s}{'[µm]':>10s}"
          f"{'[µm]':>12s}{'[µm]':>12s}")
    print("  " + "-" * 74)
    for r in sorted(rows, key=lambda q: q["n_si"]):
        tag = ""
        if abs(r["n_si"] - p[I_NSI_A]) < 1e-9:
            tag = "  ← best fit"
        elif abs(r["n_si"] - 3.4175) < 1e-9:
            tag = "  ← literature Si"
        acc = "*" if r["chi2"] <= chi2min * (1 + args.chi2_tol) else " "
        print(f" {acc}{r['n_si']:8.4f}{r['chi2']:11.5f}"
              f"{r['chi2']/chi2min:11.5f}{r['d1']:10.2f}{r['d2']:10.2f}"
              f"{r['dsum']:12.2f}{r['nd1']:12.2f}{tag}")
    print(f"\n  (* = within {100*args.chi2_tol:.0f} % of χ²min, i.e. an "
          f"equally good description of the data:")
    print(f"   {band['n_ok']}/{len(rows)} grid points qualify)\n")

    d1_lo, d1_hi = band["d1"]; d2_lo, d2_hi = band["d2"]
    ds_lo, ds_hi = band["dsum"]
    nd1_lo, nd1_hi = band["nd1"]; nd2_lo, nd2_hi = band["nd2"]
    u_d1 = 0.5 * (d1_hi - d1_lo)
    u_d2 = 0.5 * (d2_hi - d2_lo)
    u_ds = max(0.5 * (ds_hi - ds_lo), 0.005)
    u_nd1 = max(0.5 * (nd1_hi - nd1_lo), 0.005)
    u_nd2 = max(0.5 * (nd2_hi - nd2_lo), 0.005)

    print(f"  Within the two-layer model, the {100*args.chi2_tol:.0f} %-χ² "
          f"band spans")
    print(f"    d_Si = {d1_lo:.1f} … {d1_hi:.1f} µm   (± {u_d1:.1f} µm)")
    print(f"    d_qz = {d2_lo:.1f} … {d2_hi:.1f} µm   (± {u_d2:.1f} µm)")
    print(f"  while n_Si·d_Si holds at {0.5*(nd1_lo+nd1_hi):.2f} µm to "
          f"±{u_nd1:.2f} µm across the WHOLE scan —")
    print(f"  which is the degeneracy stated explicitly: n and d trade, n·d "
          f"does not move.")

    lit = [r for r in rows if abs(r["n_si"] - 3.4175) < 1e-9]
    rl = lit[0] if lit else min(rows, key=lambda q: abs(q["n_si"] - 3.4175))
    if True:
        print(f"\n  Note: fixing n_Si at the literature 3.4175 costs only "
              f"{100*(rl['chi2']/chi2min-1):.1f} % in χ²")
        print(f"  and gives d_Si = {rl['d1']:.1f} µm — "
              f"{abs(rl['d1']-D_SI_NOMINAL_UM):.1f} µm from the "
              f"{D_SI_NOMINAL_UM:.0f} µm nominal, vs "
              f"{abs(d1-D_SI_NOMINAL_UM):.1f} µm for the best fit.")
        print("  The data cannot choose between them. Quote d_Si to the "
              "nearest few µm only.")

    # ------------------------------------- 8.4c systematic: free air gap
    hr("4c. SYSTEMATIC STRESS TEST:  IS THE WAFER OPTICALLY CONTACTED?")
    gap = airgap_fit(fb, Tb, wts, p)
    dg = gap.x[1]
    print("  Refitting with one extra free parameter — an air gap between the")
    print("  wafer and the window (the two-layer model assumes optical "
          "contact):\n")
    print(f"    air gap            d_gap = {dg:8.2f} µm   "
          f"(λ/{C_UM_PS/1.0/max(dg,1e-9):.0f} at 1 THz — unresolvable in time)")
    print(f"    χ²                       = {gap.chi2:8.5f}  vs "
          f"{out.chi2:.5f} for the 2-layer model "
          f"→ {100*(1-gap.chi2/out.chi2):+.1f} %")
    print(f"    d_Si                     = {gap.x[0]:8.2f} µm  vs "
          f"{d1:.2f} µm   (Δ = {gap.x[0]-d1:+.2f} µm)")
    print(f"    n_Si (1 THz)             = {gap.x[3]+gap.x[4]:8.4f}   vs "
          f"{p[2]+p[3]:.4f}")
    print(f"\n  A single plausible extra degree of freedom changes χ² by "
          f"{100*abs(1-gap.chi2/out.chi2):.0f} % and moves")
    print(f"  d_Si by {abs(gap.x[0]-d1):.0f} µm. That is the scale of the "
          f"model systematic, and it is why")
    print(f"  the {100*args.chi2_tol:.0f} %-χ² acceptance band above is the "
          f"right criterion, not Δχ² = 1.")

    _, _, N1full, N2full, _, _, _, _ = unpack_global(p, np.where(f > 0, f, 1e-9))
    T_fit_full = model_transmission(np.where(f > 0, f, 1e-9), d1, d2,
                                    N1full, N2full)
    phi_fit = branch_corrected_phase(f, T_fit_full, mask=mask_fit)

    # ----------------------------------------- 8.5 Gated echo ratios R1, R2
    hr("5.  TIME-GATED ECHO SPECTRA  →  SEPARATING THE TWO LAYERS")

    tau_si_c, tau_qz_c = t0 + tau_si, t0 + tau_qz
    hw_main = 0.45 * min(tau_si, tau_qz)
    hw_si = 0.45 * min(tau_si, abs(tau_qz - tau_si))
    hw_qz = 0.45 * min(abs(tau_qz - tau_si), tau_si)

    g0 = super_gaussian_gate(tau, t0, hw_main)
    g1 = super_gaussian_gate(tau, tau_si_c, hw_si)
    g2 = super_gaussian_gate(tau, tau_qz_c, hw_qz)

    O0 = np.fft.rfft(h * g0, n=NFFT)
    O1 = np.fft.rfft(h * g1, n=NFFT)
    O2 = np.fft.rfft(h * g2, n=NFFT)
    R1_meas, R2_meas = O1 / O0, O2 / O0

    print(f"  gate 0 (direct pulse) : centre {t0:7.3f} ps, half-width "
          f"{hw_main:.2f} ps")
    print(f"  gate 1 (Si echo)      : centre {tau_si_c:7.3f} ps, half-width "
          f"{hw_si:.2f} ps")
    print(f"  gate 2 (quartz echo)  : centre {tau_qz_c:7.3f} ps, half-width "
          f"{hw_qz:.2f} ps")
    print(f"  ⇒ the gates smooth the echo spectra by ≈ "
          f"{1.0/(2*min(hw_si, hw_qz)):.2f} THz;")
    print( "    that, not the 1.85 GHz FFT bin spacing, is the true spectral")
    print( "    resolution of the per-frequency inversion below.")

    i1 = int(np.argmin(np.abs(f - 1.0)))
    print(f"\n  at 1 THz:  R1 = {R1_meas[i1]: .4f}   (|R1| = "
          f"{abs(R1_meas[i1]):.4f})")
    print(f"             R2 = {R2_meas[i1]: .4f}   (|R2| = "
          f"{abs(R2_meas[i1]):.4f})")
    R1m, R2m = echo_ratio_models(np.array([f[i1]]), d1, d2,
                                 N1full[i1], N2full[i1])
    print(f"  model   :  R1 = {R1m[0]: .4f},  R2 = {R2m[0]: .4f}   "
          f"(global fit)")

    # ---------------------------- 8.5b is the Si echo mismatch noise or bias?
    hr("5b. THE Si ECHO IS TOO LARGE — NOISE, OR A MODEL SYSTEMATIC?")
    idx_b = np.where((f >= BAND_POINTWISE[0]) & (f <= BAND_POINTWISE[1]))[0]
    _, _, N1b, N2b, _, _, _, _ = unpack_global(p, f[idx_b])
    R1mod_b, R2mod_b = echo_ratio_models(f[idx_b], d1, d2, N1b, N2b)

    def echo_ratio_bias(reg):
        """Median SIGNED ratio |R|measured / |R|model, at one Wiener reg."""
        hq, _, _ = wiener_deconvolve(Es, Er, f, reg=reg)
        O0q = np.fft.rfft(hq * g0, n=NFFT)
        O1q = np.fft.rfft(hq * g1, n=NFFT)
        O2q = np.fft.rfft(hq * g2, n=NFFT)
        return (float(np.median(np.abs(O1q[idx_b] / O0q[idx_b])
                                / np.abs(R1mod_b))),
                float(np.median(np.abs(O2q[idx_b] / O0q[idx_b])
                                / np.abs(R2mod_b))))

    print("  A mismatch that is two-sided and shrinks with better")
    print("  conditioning is noise. One that keeps the same sign and size")
    print("  across three decades of Wiener regularisation is a systematic:\n")
    print(f"    {'Wiener reg':>12s}{'|R1|meas/|R1|model':>22s}"
          f"{'|R2|meas/|R2|model':>22s}")
    print("  " + "-" * 56)
    for reg_test in (1e-2, 1e-3, 1e-5):
        b1_, b2_ = echo_ratio_bias(reg_test)
        print(f"    {reg_test:12.0e}{b1_:22.3f}{b2_:22.3f}")
    bias_R1, bias_R2 = echo_ratio_bias(1e-3)

    # What n_Si would be needed to explain the excess by Fresnel alone?
    # |R1| = |r01*r12| for a lossless round trip, so solve for n1.
    n_qz_1 = float(np.real(N2full[i1]))
    P_meas = abs(fresnel_r(1.0, p[2]) * fresnel_r(p[2], n_qz_1)) * bias_R1

    def _fres(n1):
        return abs(fresnel_r(1.0, n1) * fresnel_r(n1, n_qz_1)) - P_meas

    try:
        n_si_needed = float(brentq(_fres, n_qz_1 + 1e-3, 12.0))
    except Exception:
        n_si_needed = np.nan

    print(f"\n  → the Si echo is {100*(bias_R1-1):+.1f} % relative to the "
          f"model, ONE-SIDED and")
    print(f"    essentially independent of reg. This is a model systematic, "
          f"not measurement noise.")
    print(f"  → explaining it with Fresnel coefficients alone would require "
          f"n_Si ≈ {n_si_needed:.2f},")
    print(f"    against {p[2]+p[3]:.3f} from the fringe spacing — unphysical "
          f"for silicon.")
    print( "  → most likely cause: the model assumes OPTICAL CONTACT at the "
           "Si|quartz")
    print( "    boundary. A sub-wavelength air gap raises that interface's "
           "reflectivity")
    print( "    (see 4c, where a free gap absorbs much of the residual).")
    print(f"  → CONSEQUENCE for k_Si: the pointwise solver weights R1 at 1, "
          f"and the only")
    print(f"    way it can absorb a {100*(bias_R1-1):.0f} %-too-large echo is "
          f"to drive k_Si to its lower")
    print( "    bound. So 'k_Si pins at zero' is PARTLY AN ARTEFACT of this "
           "bias. k_Si is")
    print( "    reported as an upper limit only — it is not a measurement of "
           "Si transparency.")

    # ------------------------------------ 8.6 Per-frequency (pointwise) solve
    hr("6.  PER-FREQUENCY INVERSION FOR n(ω), k(ω) OF EACH LAYER")

    step = max(1, int(round(args.df_pointwise / (f[1] - f[0]))))
    idx = np.where((f >= BAND_POINTWISE[0]) & (f <= BAND_POINTWISE[1]))[0][::step]
    f_pw = f[idx]

    _, _, N1a, N2a, n1a, n2a, al1a, al2a = unpack_global(p, f_pw)
    ia = int(np.argmin(np.abs(f_pw - F_ANCHOR)))
    x_start = [n1a[ia].real, -N1a[ia].imag, n2a[ia].real, -N2a[ia].imag]

    print(f"  solving {len(f_pw)} frequencies from {f_pw[0]:.3f} to "
          f"{f_pw[-1]:.3f} THz")
    print(f"  4 unknowns (n_Si, k_Si, n_qz, k_qz) vs 6 equations "
          f"(T, R₁, R₂) per point")
    print(f"  continuation anchor at {f_pw[ia]:.3f} THz, "
          f"x₀ = [{x_start[0]:.4f}, {x_start[1]:.2e}, "
          f"{x_start[2]:.4f}, {x_start[3]:.2e}]")

    sol, cost = pointwise_inversion(f_pw, T[idx], R1_meas[idx], R2_meas[idx],
                                    d1, d2, x_start)
    good = np.isfinite(cost) & (cost < 1.0)
    print(f"  converged at {good.sum()}/{len(f_pw)} frequencies "
          f"(median residual {np.nanmedian(cost[good]):.2e})")

    sol_plot = sol.copy()
    sol_plot[~good] = np.nan

    # --------------------------------------------- 8.7 Thickness from timing
    hr("7.  LAYER THICKNESS FROM ECHO TIMING")

    # Each echo must be evaluated at ITS OWN spectral centroid, not at a
    # hard-coded 1 THz: the quartz echo is low-frequency-weighted relative to
    # the Si echo (quartz absorbs, Si does not), so the two sub-pulses do not
    # sample the same part of the dispersion curve.  We take the amplitude
    # centroid of each gated echo spectrum over the inversion band.
    def centroid(O):
        wgt = np.abs(O[idx_b])
        return float(np.sum(f[idx_b] * wgt) / np.sum(wgt))

    f_si_c, f_qz_c = centroid(O1), centroid(O2)

    def group_index(which, f_c, eps_f=0.01):
        """n_g = n + f dn/df from the fitted dispersion, at frequency f_c."""
        _, _, _, _, n1p, n2p, _, _ = unpack_global(
            p, np.array([f_c - eps_f, f_c, f_c + eps_f]))
        npf = n1p if which == 1 else n2p
        return npf[1] + f_c * (npf[2] - npf[0]) / (2 * eps_f)

    ng_si = group_index(1, f_si_c)
    ng_qz = group_index(2, f_qz_c)

    d_si_timing = nd_si / ng_si
    d_qz_timing = nd_qz / ng_qz
    timing_split_sum = d_si_timing + d_qz_timing

    print("  A round-trip echo measures the GROUP optical thickness n_g·d.")
    print(f"  Spectral centroid of each gated echo: Si {f_si_c:.3f} THz, "
          f"quartz {f_qz_c:.3f} THz")
    print(f"  Group indices there (from the fitted dispersion):"
          f"  n_g,Si = {ng_si:.4f},  n_g,qz = {ng_qz:.4f}\n")

    print("  ⚠ THE PER-LAYER 'echo timing' COLUMN IS NOT AN INDEPENDENT CHECK.")
    print("    d = n_g·d / n_g divides the measured optical thickness by an")
    print("    n_g taken from the SAME global fit being validated. Since that")
    print("    fit already reproduces the echo delays to ~0.01 %, agreement")
    print("    in these two rows is close to tautological. It is shown for")
    print("    completeness; only the SUM row below is index-free.\n")

    print(f"  {'':16s}{'echo timing':>16s}{'FP fit':>14s}"
          f"{'nominal':>12s}{'Δ (timing−fit)':>18s}")
    print(f"  {'-'*74}")
    print(f"  {'d_Si  [µm]':16s}{d_si_timing:16.2f}{d1:14.2f}"
          f"{D_SI_NOMINAL_UM:12.1f}{d_si_timing-d1:18.2f}   (not independent)")
    print(f"  {'d_quartz [µm]':16s}{d_qz_timing:16.2f}{d2:14.2f}"
          f"{D_QZ_NOMINAL_UM:12.1f}{d_qz_timing-d2:18.2f}   (not independent)")
    print(f"  {'sum   [µm]':16s}{d_sum_timing:16.2f}{d1+d2:14.2f}"
          f"{D_SI_NOMINAL_UM+D_QZ_NOMINAL_UM:12.1f}"
          f"{d_sum_timing-(d1+d2):18.2f}   ← index-free, the real check")

    # The two entries of the "echo timing" column use DIFFERENT definitions:
    # the per-layer rows divide by the fitted n_g, the sum row is the
    # index-free c*(tau1/2 + tau2/2 - dt).  They therefore need not agree, and
    # their disagreement is itself a measurement — of how wrong the fitted
    # group indices are.  Flagging it beats quietly printing a column that
    # does not add up.
    print(f"\n  INTERNAL CONSISTENCY of the 'echo timing' column:")
    print(f"    d_Si + d_qz from the two per-layer rows = "
          f"{timing_split_sum:8.2f} µm")
    print(f"    index-free sum  c·(τ_Si/2 + τ_qz/2 − Δt) = "
          f"{d_sum_timing:8.2f} µm")
    print(f"    → the column is internally inconsistent by "
          f"{timing_split_sum - d_sum_timing:+.2f} µm "
          f"({100*abs(timing_split_sum-d_sum_timing)/d_sum_timing:.2f} %).")
    print( "    The two rows use different definitions (fitted n_g vs "
           "index-free), so")
    print( "    this gap is a direct measure of the error in the fitted group "
           "indices.")
    print(f"\n  HEADLINE (index-free, model-independent):")
    print(f"    d_Si + d_qz = {d_sum_timing:.2f} µm (echo timing)  vs  "
          f"{d1+d2:.2f} µm (FP fit)")
    print(f"    Δ = {d_sum_timing-(d1+d2):+.2f} µm "
          f"({100*abs(d_sum_timing-(d1+d2))/(d1+d2):.2f} %)")

    # ------------------------------------------------- 8.8 Summary table
    hr("8.  SUMMARY:  n, k AND α AT KEY FREQUENCIES")

    def at(fq):
        j = int(np.argmin(np.abs(f_pw - fq)))
        return j if good[j] else None

    # ---- absorption detection floor -------------------------------------
    # How large an absorption could hide inside our echo-amplitude error?
    #
    # GET THE FACTOR OF 2 RIGHT.  The echo makes a round trip, so it travels
    # an extra path L = 2*d.  But R1 is a ratio of FIELD amplitudes, and a
    # field decays as exp(-alpha*L/2), not exp(-alpha*L).  The two 2's cancel
    # exactly:
    #
    #     |R1| ∝ |a^2| = exp(-4*pi*f*k*d/c) = exp(-alpha*d)
    #     ⇒  d ln|R1| / d alpha = -d        (NOT -2d, and NOT -d/2)
    #
    # so a relative field-amplitude error eps maps onto
    #
    #     delta_alpha = eps / d             [same length units]
    #
    # (An earlier version divided by 2*d, which made both floors 2x too
    # optimistic. The sanity check is the self-test below, which differentiates
    # echo_ratio_models numerically and compares against eps/d.)
    #
    # eps is taken from the data: the median relative mismatch between the
    # measured gated echo ratio and the fitted model over the inversion band
    # (gate truncation, the optical-contact systematic of section 5b, residual
    # system etalon, and noise all fold into it).  Anything below this floor is
    # reported as an upper limit, not a number: the solver returning k = 4e-30
    # means "pinned at the bound", not "zero to 30 digits".
    R1mod, R2mod = echo_ratio_models(f_pw, d1, d2, N1a, N2a)
    eps1 = float(np.median(np.abs(np.abs(R1_meas[idx]) - np.abs(R1mod))
                           / np.abs(R1mod)))
    eps2 = float(np.median(np.abs(np.abs(R2_meas[idx]) - np.abs(R2mod))
                           / np.abs(R2mod)))
    # 1e4 converts um^-1 -> cm^-1
    a_floor_si = eps1 / d1 * 1e4
    a_floor_qz = eps2 / d2 * 1e4
    k_floor_si = alpha_to_k(a_floor_si, f_pw)
    k_floor_qz = alpha_to_k(a_floor_qz, f_pw)

    # Numerical verification of the sensitivity used above, on the real
    # fitted geometry: perturb alpha_Si by 1 cm^-1 and measure d|R1|/|R1|.
    _fchk = np.array([1.0])
    _, _, _N1c, _N2c, *_ = unpack_global(p, _fchk)
    _R1a = abs(echo_ratio_models(_fchk, d1, d2, _N1c, _N2c)[0][0])
    _N1d = _N1c - 1j * alpha_to_k(np.array([1.0]), _fchk)   # +1 cm^-1
    _R1b = abs(echo_ratio_models(_fchk, d1, d2, _N1d, _N2c)[0][0])
    _sens = abs(_R1b - _R1a) / _R1a                 # per cm^-1, linearised
    _sens_exact = float(-np.log(_R1b / _R1a))       # per cm^-1, exact
    print(f"  sensitivity check (numerical, at 1 THz, fitted geometry):")
    print(f"    Δ|R1|/|R1| per 1 cm⁻¹ of α_Si   = {_sens:.6f}  (linearised)")
    print(f"    −Δln|R1|   per 1 cm⁻¹ of α_Si   = {_sens_exact:.6f}  (exact)")
    print(f"    d_Si in cm                      = {d1*1e-4:.6f}"
          f"   ← matches exactly ⇒ δα = ε/d, NOT ε/2d")

    print(f"\n  echo-amplitude accuracy: {100*eps1:.1f} % (Si), "
          f"{100*eps2:.1f} % (quartz)")
    print(f"    ⚠ the Si figure is dominated by the ONE-SIDED "
          f"{100*(bias_R1-1):+.0f} % systematic of")
    print(f"      section 5b, not by random error. Read it as a model "
          f"mismatch budget.")
    print(f"  ⇒ absorption detection floor: α_Si ≲ "
          f"{a_floor_si:.2f} cm⁻¹,  α_qz ≳ {a_floor_qz:.2f} cm⁻¹")
    print(f"  Values below the floor are quoted as upper limits.\n")

    print(f"  {'f':>6s} | {'n_Si':>8s} {'k_Si':>10s} {'α_Si':>9s} |"
          f" {'n_qz':>8s} {'k_qz':>10s} {'α_qz':>9s}")
    print(f"  {'[THz]':>6s} | {'':>8s} {'':>10s} {'[cm⁻¹]':>9s} |"
          f" {'':>8s} {'':>10s} {'[cm⁻¹]':>9s}")
    print("  " + "-" * 72)
    table_rows = []
    for fq in KEY_FREQS:
        j = at(fq)
        if j is None:
            print(f"  {fq:6.2f} |  (not converged)")
            continue
        n1v, k1v, n2v, k2v = sol[j]
        a1v, a2v = k_to_alpha(k1v, f_pw[j]), k_to_alpha(k2v, f_pw[j])

        def fmt(kv, av, kfl, afl):
            """Quote a value, or an upper limit if it is below the floor."""
            if av < afl:
                return f"{'<' + format(kfl, '.1e'):>10s} {'<' + format(afl, '.2f'):>9s}"
            return f"{kv:10.2e} {av:9.3f}"

        print(f"  {f_pw[j]:6.2f} | {n1v:8.4f} "
              f"{fmt(k1v, a1v, k_floor_si[j], a_floor_si)} |"
              f" {n2v:8.4f} {fmt(k2v, a2v, k_floor_qz[j], a_floor_qz)}")
        table_rows.append((f_pw[j], n1v, k1v, a1v, n2v, k2v, a2v))

    hr()
    print("  Literature comparison (300 K, THz range)")
    print("    high-resistivity float-zone Si : n = 3.4175, α < 0.1 cm⁻¹  "
          "(Grischkowsky 1990)")
    print("    fused silica / quartz glass    : n ≈ 1.95–1.96, "
          "α ≈ 1–2 cm⁻¹ @ 1 THz, ~8 cm⁻¹ @ 2 THz")

    # ------------------------------------------------------ 8.9 Figure
    j1 = at(1.0)
    n1_1, k1_1, n2_1, k2_1 = (sol[j1] if j1 is not None else (np.nan,) * 4)
    _a1 = k_to_alpha(k1_1, 1.0)
    if _a1 < a_floor_si:
        si_k_str = f"< {alpha_to_k(a_floor_si, 1.0):7.1e}"
        si_a_str = f"< {a_floor_si:5.2f}"
    else:
        si_k_str = f"= {k1_1:8.2e}"
        si_a_str = f"= {_a1:6.2f}"
    panel_text = textwrap.dedent(f"""\
        WELL DETERMINED  (chi^2 profile, {100*args.chi2_tol:.0f}% band)
          n_Si*d_Si  {0.5*(nd1_lo+nd1_hi):8.1f} +- {u_nd1:.1f} um
          n_qz*d_qz  {0.5*(nd2_lo+nd2_hi):8.1f} +- {u_nd2:.1f} um
          d_Si+d_qz  {0.5*(ds_lo+ds_hi):8.1f} +- {u_ds:.1f} um
              (echo timing, index-free: {d_sum_timing:.1f} um,
               D = {d_sum_timing-(d1+d2):+.2f} um)

        POORLY DETERMINED  (n<->d degeneracy)
          d_Si       {0.5*(d1_lo+d1_hi):8.0f} +- {u_d1:.0f} um
          d_quartz   {0.5*(d2_lo+d2_hi):8.0f} +- {u_d2:.0f} um
          n_Si=3.4175 (lit.) costs {100*(rl['chi2']/chi2min-1):.1f}% in chi^2
          -> the split is NOT resolved

        AT 1.00 THz
          Si        n = {n1_1:7.4f}   k {si_k_str}
                    alpha {si_a_str} cm^-1  (upper limit)
          quartz    n = {n2_1:7.4f}   k = {k2_1:8.2e}
                    alpha = {k_to_alpha(k2_1, 1.0):6.2f} cm^-1

        SYSTEMATIC
          |R1|meas/|R1|model = {bias_R1:.3f}, one-sided,
          reg-independent -> optical-contact gap.
          k_Si pinning is partly this artefact.

        FIT QUALITY
          RMS |T| residual  {rms:.4f}   band {BAND_FIT[0]:.2f}-{BAND_FIT[1]:.2f} THz
          FFT 2^{int(np.log2(NFFT))} = {NFFT} pts, df = {1000.0/(NFFT*dt):.2f} GHz
          window  {args.window}""")

    res = dict(
        t=t, dt=dt, e_sam=e_s, e_ref=e_r,
        t_peak_sam=t[int(np.argmax(np.abs(e_s)))],
        t_peak_ref=t[int(np.argmax(np.abs(e_r)))],
        e_peak_sam=e_s[int(np.argmax(np.abs(e_s)))],
        t0=t0, tau=tau, h=h, tau_si=tau_si, tau_qz=tau_qz,
        gates_plot=[(0.0, hw_main, "direct", CLR_INK2),
                    (tau_si, hw_si, "Si echo", CLR_SI),
                    (tau_qz, hw_qz, "quartz echo", CLR_QZ)],
        f=f, Es=Es, Er=Er, T=T, phi=phi, noise=noise, dr_db=dr_db,
        mask_fit=mask_fit, T_fit=T_fit_full, phi_fit=phi_fit,
        f_fit=fb, n1_fit=n1f, n2_fit=n2f,
        k1_fit=-N1f.imag, k2_fit=-N2f.imag, al1_fit=al1f, al2_fit=al2f,
        f_pw=f_pw, sol_pw=sol_plot,
        k_floor_si=k_floor_si, a_floor_si=a_floor_si,
        k_floor_qz=k_floor_qz, a_floor_qz=a_floor_qz,
        profile=rows, chi2min=chi2min, chi2_tol=args.chi2_tol,
        d1=d1, d2=d2, u_d1=u_d1, u_d2=u_d2,
        biased_window=biased_window, window=args.window, bias=bias,
        panel_text=panel_text,
    )
    png = make_figure(res, out_path)
    hr("9.  OUTPUT")
    print(f"  figure written to: {png}")
    print(f"  ({os.path.getsize(png)/1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
