#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
THz-TDS analysis of a Si-wafer-on-quartz-window stack
Menlo Systems ScanControl ASCII waveforms
================================================================================

Sample   : Si wafer (nominally ~385 um) on a quartz window (nominally ~1 mm),
           measured under a 15 min nitrogen purge.
Reference: empty (purged) beam path.

--------------------------------------------------------------------------------
HOW THE UNDER-DETERMINED PROBLEM IS MADE WELL-POSED -- read before trusting output
--------------------------------------------------------------------------------
A single normal-incidence transmission measurement gives ONE complex number
T(omega) per frequency, i.e. 2 real observables. A two-layer stack has FOUR
unknown optical constants (n_Si, k_Si, n_qz, k_qz) plus two thicknesses. A
point-by-point inversion for "n and k of each layer separately" is therefore
UNDER-DETERMINED: code that claims to do it by brute force at a single frequency
is fitting noise.

What rescues this dataset is that the two layers ring at different periods.
Their internal Fabry-Perot round trips land at different delays and imprint
etalon fringes of different period:

    Si     round trip  2*n*d/c ~  8.8 ps  ->  fringe period ~114 GHz
    quartz round trip  2*n*d/c ~ 13.0 ps  ->  fringe period ~ 77 GHz

That gives FOUR independent observables -- two fringe periods and two fringe
contrasts. The periods fix n*d for each layer; the contrasts are set by the
Fresnel products r(Si|air)*r(Si|qz) and r(qz|air)*r(qz|Si), which depend on the
indices but NOT on the thicknesses. Four observables, four unknowns
(n_Si, n_qz, d_Si, d_qz) -> well-posed, and absorption then follows from the
residual amplitude. This is the physical content of the global fit in STEP 6.

Three estimators are computed and all reported, so they cross-check each other:

  (A) ECHO TIMING -- cross-correlation of each internal echo against the main
      transmitted pulse gives the round-trip delay, hence the optical thickness
      n*d of each layer. This is the most robust number in the analysis: it is
      independent of gate width, of amplitude calibration and of absorption.

  (B) ECHO-RATIO PHASE -- the spectral phase of (echo / main pulse) gives the
      dispersive n(omega)*d of each layer. Layer-selective, because everything
      common to both pulses cancels in the ratio. Also amplitude-insensitive.

  (C) GLOBAL TRANSFER-MATRIX FIT -- the full 4-interface TMM fitted to complex
      T(omega), with n_j(omega) and k_j(omega) as SMOOTH splines whose knot
      spacing (0.2 THz) is far coarser than the fringe periods. The splines
      therefore cannot absorb the fringes, which are left to determine d_Si and
      d_qz. This is the only estimator that yields k, and it is the primary
      source of absolute n and d.

  NOT USED for inversion: the raw echo AMPLITUDE. It is depressed by absorption
  as well as by the Fresnel factors, and on this dataset |echo/main| is not
  monotonic in frequency, so it cannot be extrapolated to omega -> 0 to separate
  the two. It is reported as a diagnostic only. (An earlier version of this
  analysis inverted it for n; that route gave d_Si = 438 um, ~14 % high, and was
  removed.)

Outputs: fig1_time_domain.png, fig2_spectra_transmission.png,
         fig3_optical_constants.png, fig4_tmm_fit.png, thz_results.txt
"""

import os
import numpy as np
from scipy.signal import find_peaks, correlate
from scipy.signal.windows import tukey
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")               # headless: write PNGs, never open a window
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 0. CONFIGURATION
# ------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REF_FILE = os.path.join(HERE, "Reference_15min_nitrogen (1).txt")
SAM_FILE = os.path.join(HERE, "Quartz_Si_15min_nitrogen.txt")

C0 = 299_792_458.0                  # speed of light in vacuum [m/s]

# Nominal geometry quoted by the experimenter. Used ONLY to seed the echo search
# and as a final cross-check -- every thickness reported below is measured.
D_SI_NOM, D_QZ_NOM = 385e-6, 1000e-6
N_SI_GUESS, N_QZ_GUESS = 3.4175, 1.95   # literature seeds (HR-Si, fused quartz)

ZERO_PAD = 1 << 15                  # FFT length after zero padding (32768)

# Usable band, read off the measured dynamic-range curve in fig2: below ~0.25 THz
# the antenna response dies, above ~2.6 THz the SAMPLE spectrum meets the noise
# floor. F_LOCK is the narrower band used while the fit is still locking onto the
# correct fringe order (fewer fringes -> broader, unimodal basin).
F_MIN, F_MAX = 0.25e12, 2.60e12
F_LOCK = 1.20e12

GATE_HW = 2.0e-12                   # echo gate half-width [s]; < half the 4.2 ps
                                    # Si-to-quartz echo separation

# --- Plot styling -------------------------------------------------------------
# Categorical palette validated for colour-vision-deficiency separation: every
# adjacent pair clears the protan/deutan/tritan and normal-vision dE floors.
C_SI, C_QZ, C_REF, C_SAM = "#1B6FE3", "#D1640A", "#158A6E", "#8B5CF6"
C_INK, C_MUTED = "#22252A", "#7A8089"

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 9, "axes.labelsize": 9.5, "axes.titlesize": 10.5,
    "axes.titleweight": "semibold", "axes.edgecolor": "#C9CDD4",
    "axes.labelcolor": C_INK, "text.color": C_INK,
    "xtick.color": C_MUTED, "ytick.color": C_MUTED,
    "axes.grid": True, "grid.color": "#E6E8EC", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def _despine(ax):
    """Recessive chart frame: drop the top and right spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ==============================================================================
# STEP 1 -- READ THE TWO DATA FILES
# ==============================================================================
def read_menlo(path):
    """
    Read a Menlo ScanControl ASCII waveform.

    Layout: '#'-prefixed header lines, then two tab-separated columns,
    'Time [ps]' and 'THz Signal [a.u.]'. The time column is an absolute stage
    coordinate (here starting near -740 ps); only differences are meaningful, so
    it is re-zeroed to the start of the record. Both scans share the same stage
    axis, which is what lets the sample/reference delay be read off directly.
    """
    raw = np.loadtxt(path, comments="#")
    t_ps, amp = raw[:, 0], raw[:, 1]

    # Verify rather than assume uniform sampling: a non-uniform stage axis would
    # silently corrupt every FFT downstream.
    step = np.diff(t_ps)
    if np.ptp(step) > 1e-6 * np.mean(step):
        raise ValueError(f"{os.path.basename(path)}: non-uniform time axis")

    return (t_ps - t_ps[0]) * 1e-12, amp, t_ps[0]


t_ref, E_ref_raw, t0_ref = read_menlo(REF_FILE)
t_sam, E_sam_raw, t0_sam = read_menlo(SAM_FILE)

# The two scans must sit on a common stage grid, else the measured pulse delay --
# and hence every optical thickness -- is offset by the grid mismatch.
if not np.allclose(t_ref, t_sam) or abs(t0_ref - t0_sam) > 1e-9:
    raise ValueError("Reference and sample scans are not on a common time axis")

t = t_ref
dt = float(np.mean(np.diff(t)))
N = len(t)

# Remove DC offset using the quiet pre-pulse stretch (the main pulse arrives
# ~17 ps into the record, so the first 10 ps is signal-free baseline).
n_base = int(10e-12 / dt)
E_ref = E_ref_raw - E_ref_raw[:n_base].mean()
E_sam = E_sam_raw - E_sam_raw[:n_base].mean()

i_pk_ref = int(np.argmax(np.abs(E_ref)))
i_pk_sam = int(np.argmax(np.abs(E_sam)))
t_pk_ref, t_pk_sam = t[i_pk_ref], t[i_pk_sam]
dt_main = t_pk_sam - t_pk_ref       # group delay of the whole stack vs air [s]


# ==============================================================================
# STEP 2 -- FFT BOTH WAVEFORMS WITH ZERO PADDING
# ==============================================================================
# Zero padding adds no true spectral resolution (that is fixed by the 80 ps
# record, ~12.5 GHz). It interpolates onto a finer grid, which matters here
# because we must resolve 77-114 GHz etalon fringes and because phase unwrapping
# is far safer on a densely sampled phase.
#
# A Tukey window is applied first: the record is truncated while ~1 % ringing is
# still present, and an abrupt cut would ring across the whole spectrum. The
# taper is confined to the record edges so the pulse and all echoes pass through
# unmodified.
win = tukey(N, alpha=0.10)

f = np.fft.rfftfreq(ZERO_PAD, dt)
w = 2.0 * np.pi * f
Eref_w = np.fft.rfft(E_ref * win, ZERO_PAD)
Esam_w = np.fft.rfft(E_sam * win, ZERO_PAD)

# Dynamic range: FFT a signal-free stretch of the reference for the noise floor,
# so the usable band is measured rather than guessed.
n_noise = int(12e-12 / dt)
noise_w = np.abs(np.fft.rfft(E_ref[:n_noise] * tukey(n_noise, 0.1), ZERO_PAD))
noise_floor = np.maximum(noise_w, 1e-30)
dr_ref = 20 * np.log10(np.abs(Eref_w) / noise_floor)
dr_sam = 20 * np.log10(np.abs(Esam_w) / noise_floor)

band = (f >= F_MIN) & (f <= F_MAX)


# ==============================================================================
# STEP 3 -- COMPLEX TRANSMISSION  T(omega) = E_sample / E_reference
# ==============================================================================
# Dividing by the reference cancels the emitter spectrum, the detector response
# and the purged beam path, leaving the response of the stack alone.
T = Esam_w / np.where(np.abs(Eref_w) > 0, Eref_w, np.nan)
T_abs = np.abs(T)

# Phase unwrapping with a delay anchor: raw np.unwrap on a steep phase can slip
# 2*pi at any noisy sample, and the error then propagates to every higher
# frequency. Removing the known bulk delay first leaves a small, slowly varying
# residual that unwraps unambiguously; the linear term is restored analytically.
phi_resid = np.unwrap(np.angle(T * np.exp(1j * w * dt_main)))
i_anchor = int(np.argmax(f >= F_MIN))
phi_resid -= 2 * np.pi * np.round(phi_resid[i_anchor] / (2 * np.pi))
T_phase = phi_resid - w * dt_main


# ==============================================================================
# STEP 4 -- TRANSFER-MATRIX (FABRY-PEROT) MODEL OF THE Si/QUARTZ STACK
# ==============================================================================
# Convention: fields vary as exp(+i*omega*t - i*k*z), matching numpy's forward
# rfft kernel exp(-i*omega*t) used above. A complex index n~ = n - i*k then
# propagates as exp(-i*n~*omega*d/c), whose magnitude exp(-k*omega*d/c) decays
# for k > 0, as it must.
#
# Stack (beam side first): air | Si | quartz | air.
#
# LAYER ORDER IS IRRELEVANT HERE: the first-pulse amplitude is the product
# t01*t12*t23, symmetric under reversal, and each echo ratio is the product of
# the two Fresnel coefficients bounding that layer -- also order-independent. So
# nothing below depends on which face the beam strikes first.

def interface_matrix(na, nb):
    """2x2 matrix for a single interface a->b at normal incidence."""
    r = (na - nb) / (na + nb)
    tt = 2.0 * na / (na + nb)
    M = np.empty(na.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = 1.0
    M[..., 0, 1] = r
    M[..., 1, 0] = r
    M[..., 1, 1] = 1.0
    return M / tt[..., None, None]


def propagation_matrix(n_layer, d, omega):
    """
    2x2 matrix for propagation through thickness d of a layer.

    Sign convention: with the interface matrices written as above and the product
    accumulated left-to-right (incident medium first), the layer matrix carries
    exp(+i*delta) in the [0,0] slot. Since transmission is recovered as
    t = 1/M[0,0], this yields the physical exp(-i*delta) on the transmitted
    field, i.e. decay exp(-k*omega*d/c) for k > 0. The opposite choice silently
    produces GAIN instead of loss, so the convention is pinned down at runtime by
    _selftest_tmm() against the analytic Airy formula.
    """
    delta = n_layer * omega * d / C0
    M = np.zeros(delta.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = np.exp(+1j * delta)
    M[..., 1, 1] = np.exp(-1j * delta)
    return M


def tmm_transmission(omega, n_layers, d_layers, n_in=1.0, n_out=1.0):
    """
    Field transmission t(omega) of a stack at normal incidence.

    Includes ALL Fabry-Perot reflections in every layer to infinite order -- that
    infinite resummation is exactly what the matrix product performs.
    """
    omega = np.atleast_1d(omega)
    one = np.ones_like(omega, dtype=complex)
    media = [n_in * one] + [n * one for n in n_layers] + [n_out * one]

    M = np.zeros(omega.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = 1.0
    M[..., 1, 1] = 1.0
    for j in range(len(media) - 1):
        M = M @ interface_matrix(media[j], media[j + 1])
        if j < len(d_layers):
            M = M @ propagation_matrix(media[j + 1], d_layers[j], omega)
    return 1.0 / M[..., 0, 0]


def model_T(omega, n_si, n_qz, d_si, d_qz):
    """
    Modelled MEASURED transmission: stack transmission divided by the reference,
    which is the same geometric length (d_si + d_qz) of air. That division is why
    the free-space factor exp(-i*omega*D/c) is multiplied back out.
    """
    t_stack = tmm_transmission(omega, [n_si, n_qz], [d_si, d_qz])
    return t_stack * np.exp(+1j * omega * (d_si + d_qz) / C0)


def _selftest_tmm():
    """
    Check the matrix code against the closed-form single-layer Airy formula.
    A silent sign or ordering slip in a TMM is the classic way to produce
    plausible-looking but wrong optical constants, so this runs every time.
    """
    om = 2 * np.pi * np.linspace(0.2e12, 3e12, 40)
    n1, d1 = 3.4 - 0.002j, 400e-6
    r01, r12 = (1 - n1) / (1 + n1), (n1 - 1) / (n1 + 1)
    t01, t12 = 2 / (1 + n1), 2 * n1 / (n1 + 1)
    dl = n1 * om * d1 / C0
    analytic = t01 * t12 * np.exp(-1j * dl) / (1 + r01 * r12 * np.exp(-2j * dl))
    numeric = tmm_transmission(om, [n1 * np.ones_like(om)], [d1])
    err = np.max(np.abs(analytic - numeric) / np.abs(analytic))
    assert err < 1e-12, f"TMM self-test failed (rel. err {err:.2e})"


_selftest_tmm()


# ==============================================================================
# STEP 5 -- ECHO TIMING  (estimator A)
# ==============================================================================
def gate(t, E, t_centre, hw=GATE_HW, alpha=0.4):
    """
    Isolate one pulse with a Tukey gate of FIXED width centred on it.

    The width must be identical for the main pulse and for every echo: gating
    convolves the spectrum with the gate's own transform, and only an identical
    gate cancels in the ratio taken later.
    """
    g = np.zeros_like(E)
    m = np.abs(t - t_centre) <= hw
    g[m] = E[m] * tukey(int(m.sum()), alpha)
    return g


def echo_delay(t, E, t_main, tau_guess, search=1.5e-12):
    """
    Round-trip delay of one echo, by cross-correlating the gated main pulse
    against the full waveform.

    Cross-correlation rather than peak-picking, because the transmitted pulse is
    bipolar and an echo whose Fresnel product is NEGATIVE (the quartz one) comes
    back inverted. Peak-picking on |E| then locks onto the wrong lobe and biases
    the delay by half a cycle -- worth ~90 um of optical thickness here.
    Correlating the whole pulse shape is immune to that, and the sign of the
    correlation peak directly measures the sign of the Fresnel product.

    Returns (delay, sign_of_rho).
    """
    m = gate(t, E, t_main)
    c = correlate(E, m, mode="full")
    lags = (np.arange(len(c)) - (len(E) - 1)) * dt

    sel = np.abs(lags - tau_guess) < search
    i = int(np.argmax(np.abs(c) * sel))

    # Parabolic sub-sample interpolation on |correlation|: the 33 fs sampling
    # would otherwise quantise n*d in ~5 um steps.
    y0, y1, y2 = np.abs(c[i - 1]), np.abs(c[i]), np.abs(c[i + 1])
    denom = y0 - 2 * y1 + y2
    shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    return lags[i] + shift * dt, float(np.sign(np.real(c[i])))


tau_si, sgn_si = echo_delay(t, E_sam, t_pk_sam, 2 * N_SI_GUESS * D_SI_NOM / C0)
tau_qz, sgn_qz = echo_delay(t, E_sam, t_pk_sam, 2 * N_QZ_GUESS * D_QZ_NOM / C0)

# Optical thicknesses -- the most robust numbers in this analysis.
nd_si = C0 * tau_si / 2.0
nd_qz = C0 * tau_qz / 2.0

# Independent closure check: the main-pulse delay must equal the sum of the two
# excess optical paths. This ties the echo picture to the bulk measurement and
# catches a mis-assigned echo immediately. Evaluated with the fitted thicknesses
# once those exist (see STEP 8).

# Physics predicts the echo polarities: air < n_qz < n_Si makes the Si echo
# positive (both bounding reflections same sign) and the quartz echo negative
# (opposite signs). Measured polarity is compared against that below.
sgn_si_expected, sgn_qz_expected = +1.0, -1.0
polarity_ok = (sgn_si == sgn_si_expected) and (sgn_qz == sgn_qz_expected)


def find_echoes(t, E, t_main, search=(4e-12, 32e-12), n_report=10):
    """Catalogue the visible echoes, for the annotated time-domain figure."""
    m = (t > t_main + search[0]) & (t < t_main + search[1])
    tt, aa = t[m], np.abs(E[m])
    pk, _ = find_peaks(aa, height=0.02 * np.abs(E).max(),
                       distance=max(1, int(0.8e-12 / dt)))
    order = np.argsort(aa[pk])[::-1][:n_report]
    return tt[pk][np.sort(order)] - t_main, aa[pk][np.sort(order)]


echo_dt, echo_amp = find_echoes(t, E_sam, t_pk_sam)


# ==============================================================================
# STEP 6 -- ECHO-RATIO SPECTROSCOPY  (estimator B): n(omega) per layer
# ==============================================================================
def echo_ratio_spectrum(t, E, t_main, tau):
    """
    Complex spectral ratio (echo / main pulse) for one layer.

    For a stack air|1|2|air, the first transmitted pulse and the pulse that has
    made one extra round trip inside layer j differ by exactly

        A_echo / A_main = rho_j * exp(-2i * n~_j * omega * d_j / c),
        rho_j = r(j|left) * r(j|right)

    Everything else -- emitter spectrum, detector response, the OTHER layer's
    factors -- cancels identically. That is what makes the ratio layer-selective.
    It is taken WITHIN the sample scan, so reference-path artefacts cancel too.
    """
    A_main = np.fft.rfft(gate(t, E, t_main), ZERO_PAD)
    A_echo = np.fft.rfft(gate(t, E, t_main + tau), ZERO_PAD)
    return A_echo / A_main


ratio_si = echo_ratio_spectrum(t, E_sam, t_pk_sam, tau_si)
ratio_qz = echo_ratio_spectrum(t, E_sam, t_pk_sam, tau_qz)


def echo_phase_nd(ratio, tau, sgn):
    """
    Dispersive optical thickness n(omega)*d from the echo-ratio PHASE.

        arg(ratio / sgn(rho)) = -2 * n * omega * d / c

    Dividing out sgn(rho) first is essential: a negative Fresnel product (quartz)
    contributes a flat pi that would otherwise be absorbed into n*d as a ~90 um
    bias. Phase only -- so this is immune to the amplitude problems that make the
    echo MAGNITUDE unusable for inversion (see module docstring).
    """
    r = ratio / sgn
    resid = np.unwrap(np.angle(r * np.exp(1j * w * tau)))
    resid -= 2 * np.pi * np.round(resid[int(np.argmax(f >= F_MIN))] / (2 * np.pi))
    phase = resid - w * tau
    with np.errstate(divide="ignore", invalid="ignore"):
        return -C0 * phase / (2.0 * w)      # = n(omega) * d


nd_si_spec = echo_phase_nd(ratio_si, tau_si, sgn_si)
nd_qz_spec = echo_phase_nd(ratio_qz, tau_qz, sgn_qz)

# Diagnostic only -- NOT inverted for n (see module docstring).
lowf = (f >= 0.30e12) & (f <= 0.80e12)
rho_si_meas = sgn_si * float(np.median(np.abs(ratio_si[lowf])))
rho_qz_meas = sgn_qz * float(np.median(np.abs(ratio_qz[lowf])))


# ==============================================================================
# STEP 7 -- GLOBAL TRANSFER-MATRIX FIT  (estimator C)
# ==============================================================================
# The fit is staged, because the residual landscape is periodic in thickness:
# every fringe order is a local minimum. Stages 1-2 work on the narrow band
# (0.25-1.2 THz), where fewer fringes fit inside the band and the basin around
# the true order is broad and unimodal; stage 3 then opens the full band.
#
#   stage 1 : coarse 2-D grid over (d_Si, d_qz), literature constant n,k
#             -> locks the fringe ORDER
#   stage 2 : continuous refine of d and constant n,k
#   stage 3 : full band, dispersive n(omega), k(omega) as splines + d
#
# Knot spacing (0.2 THz) is deliberately ~2x coarser than the Si fringe period
# (114 GHz) and ~2.6x coarser than the quartz one (77 GHz), so the smooth curves
# CANNOT reproduce the fringes. The fringes are thereby forced to be explained by
# d_Si and d_qz, while the fringe CONTRAST -- which is thickness-independent --
# pins n. That separation is what makes the inversion well-posed.

lock = (f >= F_MIN) & (f <= F_LOCK)
idx_lock = np.where(lock)[0][::3]
om_lock, T_lock = w[idx_lock], T[idx_lock]


def _rms(om, Td, ns, nq, ds, dq):
    return float(np.sqrt(np.mean(np.abs(model_T(om, ns, nq, ds, dq) - Td) ** 2)))


# --- stage 1: coarse grid -----------------------------------------------------
ns0 = N_SI_GUESS - 0.0005j
nq0 = N_QZ_GUESS - 0.006j
d_si_scan = np.arange(0.90, 1.101, 0.0013) * D_SI_NOM
d_qz_scan = np.arange(0.90, 1.101, 0.0013) * D_QZ_NOM

grid = np.empty((len(d_si_scan), len(d_qz_scan)))
for i, ds in enumerate(d_si_scan):
    for j, dq in enumerate(d_qz_scan):
        grid[i, j] = _rms(om_lock, T_lock, ns0, nq0, ds, dq)
i0, j0 = np.unravel_index(np.argmin(grid), grid.shape)
d_si_grid, d_qz_grid = d_si_scan[i0], d_qz_scan[j0]

# --- stage 2: continuous refine, constant n and k -----------------------------
def _resid_const(p):
    n_s, k_s, n_q, k_q, ds, dq = p
    Tm = model_T(om_lock, n_s - 1j * 10 ** k_s, n_q - 1j * 10 ** k_q, ds, dq)
    r = Tm - T_lock
    return np.concatenate([r.real, r.imag])


p_const = least_squares(
    _resid_const,
    [N_SI_GUESS, -3.0, N_QZ_GUESS, -2.3, d_si_grid, d_qz_grid],
    bounds=([2.5, -6, 1.2, -6, 0.9 * d_si_grid, 0.9 * d_qz_grid],
            [4.5, -0.5, 3.0, -0.5, 1.1 * d_si_grid, 1.1 * d_qz_grid]),
    x_scale="jac", max_nfev=600).x

# --- stage 3: full band, dispersive n and k -----------------------------------
KNOT_SPACING = 0.20e12
knots = np.arange(F_MIN, F_MAX + 0.5 * KNOT_SPACING, KNOT_SPACING)
NK = len(knots)

fit_idx = np.where(band)[0][::4]        # still ~7 samples per Si fringe
f_fit, w_fit, T_fit = f[fit_idx], w[fit_idx], T[fit_idx]

# Weight by reference amplitude so high-SNR frequencies dominate rather than
# every point counting equally.
wgt = np.abs(Eref_w[fit_idx])
wgt /= wgt.max()


def unpack(p):
    """Parameter vector -> dispersive complex indices on the fit grid."""
    n_si = np.interp(f_fit, knots, p[0:NK])
    k_si = 10.0 ** np.interp(f_fit, knots, p[NK:2 * NK])
    n_qz = np.interp(f_fit, knots, p[2 * NK:3 * NK])
    k_qz = 10.0 ** np.interp(f_fit, knots, p[3 * NK:4 * NK])
    return n_si - 1j * k_si, n_qz - 1j * k_qz, p[4 * NK], p[4 * NK + 1]


def residuals(p):
    ns, nq, ds, dq = unpack(p)
    r = (model_T(w_fit, ns, nq, ds, dq) - T_fit) * wgt
    return np.concatenate([r.real, r.imag])


n_s0, k_s0, n_q0, k_q0, d_si_2, d_qz_2 = p_const
p0 = np.concatenate([np.full(NK, n_s0), np.full(NK, k_s0),
                     np.full(NK, n_q0), np.full(NK, k_q0), [d_si_2, d_qz_2]])
lo = np.concatenate([np.full(NK, 2.5), np.full(NK, -6.0),
                     np.full(NK, 1.2), np.full(NK, -6.0),
                     [0.95 * d_si_2, 0.95 * d_qz_2]])
hi = np.concatenate([np.full(NK, 4.5), np.full(NK, -0.5),
                     np.full(NK, 3.0), np.full(NK, -0.5),
                     [1.05 * d_si_2, 1.05 * d_qz_2]])

fit = least_squares(residuals, p0, bounds=(lo, hi), method="trf",
                    x_scale="jac", max_nfev=3000)

ns_fit, nq_fit, d_si_fit, d_qz_fit = unpack(fit.x)
T_model_fit = model_T(w_fit, ns_fit, nq_fit, d_si_fit, d_qz_fit)

n_si_tmm, k_si_tmm = ns_fit.real, -ns_fit.imag
n_qz_tmm, k_qz_tmm = nq_fit.real, -nq_fit.imag

# Absorption coefficient (intensity): alpha = 2*omega*k/c, converted to 1/cm.
alpha_si_tmm = 2 * w_fit * k_si_tmm / C0 / 100.0
alpha_qz_tmm = 2 * w_fit * k_qz_tmm / C0 / 100.0

rms_res = float(np.sqrt(np.mean(np.abs(T_model_fit - T_fit) ** 2)))

# Absorption sensitivity floor. A layer of thickness d changes the single-pass
# amplitude by ~exp(-alpha*d/2), so an alpha whose effect is smaller than the fit
# residual is not measurable. Below this line the fitted alpha is an UPPER BOUND,
# not a measurement -- which is exactly the situation for Si here, where the
# optimiser simply sits on its lower bound.
_rel_res = rms_res / float(np.mean(np.abs(T_fit)))
alpha_min_si = 2 * _rel_res / d_si_fit / 100.0      # [1/cm]
alpha_min_qz = 2 * _rel_res / d_qz_fit / 100.0      # [1/cm]

# Echo-phase n(omega), now that a thickness is available to divide by.
n_si_echo = nd_si_spec / d_si_fit
n_qz_echo = nd_qz_spec / d_qz_fit

# Thickness implied by echo timing combined with the fitted index -- an estimate
# of d that uses the fit only for n, not for the fringe counting.
n_si_ref = float(np.nanmean(n_si_tmm[(f_fit >= 0.4e12) & (f_fit <= 1.6e12)]))
n_qz_ref = float(np.nanmean(n_qz_tmm[(f_fit >= 0.4e12) & (f_fit <= 1.6e12)]))
d_si_echo = nd_si / n_si_ref
d_qz_echo = nd_qz / n_qz_ref

# Closure check on the bulk delay.
excess_pred = (nd_si - d_si_fit) + (nd_qz - d_qz_fit)
excess_meas = C0 * dt_main
closure = abs(excess_pred - excess_meas) / excess_meas


# ==============================================================================
# STEP 8 -- FIGURES
# ==============================================================================
fT = f / 1e12

# --- fig1: time domain --------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.4),
                         gridspec_kw={"hspace": 0.32})

ax = axes[0]
ax.plot((t - t_pk_ref) * 1e12, E_ref, color=C_REF, label="Reference (air)")
ax.plot((t - t_pk_ref) * 1e12, E_sam, color=C_SAM, label="Sample (Si + quartz)")
ax.set_xlim(-5, 45)
ax.set_xlabel("Time relative to reference peak  (ps)")
ax.set_ylabel("THz signal  (a.u.)")
ax.set_title("Time-domain waveforms")
ax.legend(loc="upper right")
ax.annotate(f"stack delay {dt_main*1e12:.2f} ps",
            xy=(dt_main * 1e12, E_sam[i_pk_sam]),
            xytext=(26, -10), textcoords="offset points", fontsize=8,
            color=C_MUTED,
            arrowprops=dict(arrowstyle="->", color=C_MUTED, lw=0.9))
_despine(ax)

ax = axes[1]
ax.plot((t - t_pk_sam) * 1e12, E_sam, color=C_SAM, lw=1.3)
ax.set_ylim(-0.45, 0.45)
ax.set_xlim(-4, 36)
# Label heights are staggered because the predicted round trips fall as close as
# 4.2 ps apart and the text blocks would otherwise overlap.
for tau_, lab, col, ytxt in ((tau_si, "Si round trip", C_SI, 0.41),
                             (tau_qz, "quartz round trip", C_QZ, 0.25),
                             (2 * tau_si, "2x Si", C_SI, 0.41),
                             (tau_si + tau_qz, "Si + quartz", C_MUTED, 0.25),
                             (2 * tau_qz, "2x quartz", C_QZ, 0.41)):
    ax.axvline(tau_ * 1e12, color=col, lw=1.1, ls="--", alpha=0.8)
    ax.text(tau_ * 1e12 + 0.25, ytxt, f"{lab}\n{tau_*1e12:.2f} ps",
            fontsize=7, color=col, va="top")
ax.set_xlabel("Delay relative to main transmitted pulse  (ps)")
ax.set_ylabel("THz signal  (a.u.)")
ax.set_title("Internal Fabry-Perot echoes (vertical zoom; labelled by round trip)")
_despine(ax)

fig.suptitle("THz-TDS: Si wafer on quartz window", y=0.98, fontsize=11.5,
             fontweight="semibold")
fig.savefig(os.path.join(HERE, "fig1_time_domain.png"))
plt.close(fig)

# --- fig2: spectra and transmission -------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(8.4, 8.4), sharex=True,
                         gridspec_kw={"hspace": 0.25})

ax = axes[0]
ax.plot(fT, dr_ref, color=C_REF, label="Reference")
ax.plot(fT, dr_sam, color=C_SAM, label="Sample")
ax.axhline(10, color=C_MUTED, lw=0.9, ls=":")
ax.text(3.45, 12, "10 dB", fontsize=7.5, color=C_MUTED, ha="right")
ax.axvspan(F_MIN / 1e12, F_MAX / 1e12, color=C_SI, alpha=0.06, lw=0)
ax.text((F_MIN + F_MAX) / 2e12, 64, "analysis band", fontsize=8,
        color=C_MUTED, ha="center")
ax.set_xlim(0, 3.5)
ax.set_ylim(-5, 74)
ax.set_ylabel("Dynamic range  (dB)")
ax.set_title("Spectral dynamic range (sets the usable band)")
ax.legend(loc="upper right")
_despine(ax)

ax = axes[1]
ax.plot(fT[band], T_abs[band], color=C_INK, lw=1.2)
ax.set_ylabel("|T|")
ax.set_ylim(0, 1.05)
ax.set_title("Transmission magnitude  |T| = |E$_{sample}$/E$_{ref}$|   "
             "(oscillation = the two Fabry-Perot etalons)")
_despine(ax)

ax = axes[2]
ax.plot(fT[band], T_phase[band], color=C_INK, lw=1.2)
ax.set_xlabel("Frequency  (THz)")
ax.set_ylabel("arg T  (rad)")
ax.set_title("Unwrapped transmission phase")
_despine(ax)

fig.suptitle("Complex transmission of the stack", y=0.965, fontsize=11.5,
             fontweight="semibold")
fig.savefig(os.path.join(HERE, "fig2_spectra_transmission.png"))
plt.close(fig)

# --- fig3: layer-resolved optical constants -----------------------------------
# One measure per axis (never a twin y-scale): index and absorption go in
# separate rows, material is the colour series and is also direct-labelled by
# the column title, so identity is never carried by colour alone.
fig, axes = plt.subplots(2, 2, figsize=(9.8, 6.6),
                         gridspec_kw={"hspace": 0.34, "wspace": 0.24})
m = band

ax = axes[0, 0]
ax.plot(fT[m], n_si_echo[m], color=C_SI, lw=1.0, alpha=0.45,
        label="echo-ratio phase")
ax.plot(f_fit / 1e12, n_si_tmm, color=C_SI, lw=2.2, label="TMM global fit")
ax.set_ylim(3.30, 3.55)
ax.set_ylabel("refractive index  n")
ax.set_title("Silicon")
ax.legend(loc="lower right", fontsize=8)
_despine(ax)

ax = axes[0, 1]
ax.plot(fT[m], n_qz_echo[m], color=C_QZ, lw=1.0, alpha=0.45,
        label="echo-ratio phase")
ax.plot(f_fit / 1e12, n_qz_tmm, color=C_QZ, lw=2.2, label="TMM global fit")
ax.set_ylim(1.85, 2.20)
ax.set_title("Quartz")
ax.legend(loc="lower right", fontsize=8)
_despine(ax)

# Absorption panels carry an explicit sensitivity floor. Without it the Si panel
# would draw a confident-looking line at ~1e-3 cm^-1 that is really just the
# optimiser resting on its lower bound -- a measurement of nothing.
ax = axes[1, 0]
ax.axhspan(0, alpha_min_si, color=C_MUTED, alpha=0.14, lw=0)
ax.plot(f_fit / 1e12, alpha_si_tmm, color=C_SI, lw=2.2)
ax.set_ylim(0, max(3.0, 2.2 * alpha_min_si))
ax.text(0.30, alpha_min_si * 0.45,
        f"below detection limit ({alpha_min_si:.1f} cm$^{{-1}}$)",
        fontsize=7.5, color=C_MUTED, va="center")
ax.annotate("fitted $\\alpha_{Si}$ lies here:\nconsistent with zero",
            xy=(1.4, 0.0), xytext=(1.05, 1.28 * alpha_min_si), fontsize=7.5,
            color=C_SI, ha="center",
            arrowprops=dict(arrowstyle="->", color=C_SI, lw=0.9))
ax.set_xlabel("Frequency  (THz)")
ax.set_ylabel(r"absorption  $\alpha$  (cm$^{-1}$)")
_despine(ax)

ax = axes[1, 1]
ax.axhspan(0, alpha_min_qz, color=C_MUTED, alpha=0.14, lw=0)
ax.plot(f_fit / 1e12, alpha_qz_tmm, color=C_QZ, lw=2.2)
ax.set_ylim(bottom=0)
ax.text(2.55, alpha_min_qz * 2.1, "detection limit", fontsize=7.5,
        color=C_MUTED, va="center", ha="right")
ax.set_xlabel("Frequency  (THz)")
_despine(ax)

fig.suptitle("Layer-resolved optical constants", y=0.98, fontsize=11.5,
             fontweight="semibold")
fig.text(0.5, -0.01,
         "Index: thin line = model-light echo-ratio phase, thick line = global "
         "transfer-matrix fit.  Absorption comes from the TMM fit only "
         "(the echo magnitude alone cannot separate loss from Fresnel factors).",
         ha="center", fontsize=8, color=C_MUTED)
fig.savefig(os.path.join(HERE, "fig3_optical_constants.png"))
plt.close(fig)

# --- fig4: fit quality --------------------------------------------------------
fig, axes = plt.subplots(4, 1, figsize=(8.4, 9.6),
                         gridspec_kw={"hspace": 0.42,
                                      "height_ratios": [1, 0.85, 0.85, 0.7]})

ax = axes[0]
ax.plot(f_fit / 1e12, np.abs(T_fit), color=C_MUTED, lw=2.8, alpha=0.5,
        label="measured")
ax.plot(f_fit / 1e12, np.abs(T_model_fit), color=C_SI, lw=1.2, label="TMM fit")
ax.set_ylabel("|T|")
ax.set_ylim(0, 1.05)
ax.set_xlabel("Frequency  (THz)")
ax.set_title("Transfer-matrix fit to the complex transmission")
ax.legend(loc="lower left", ncol=2)
_despine(ax)

# ~2500 fringes across the band collapse into a solid mass at print size, so the
# zoom gets its own row: it carries the actual claim, that every fringe of BOTH
# etalons is matched in position and depth, beat pattern included.
ax = axes[1]
zoom = (f_fit >= 0.95e12) & (f_fit <= 1.20e12)
ax.plot(f_fit[zoom] / 1e12, np.abs(T_fit[zoom]), color=C_MUTED, lw=3.0,
        alpha=0.5, label="measured")
ax.plot(f_fit[zoom] / 1e12, np.abs(T_model_fit[zoom]), color=C_SI, lw=1.4,
        label="TMM fit")
ax.set_xlim(0.95, 1.20)
ax.set_xlabel("Frequency  (THz)")
ax.set_ylabel("|T|")
ax.set_title("Zoom, 0.95 - 1.20 THz: individual fringes of both etalons",
             fontsize=9.5)
_despine(ax)

ax = axes[2]
ax.plot(f_fit / 1e12, np.unwrap(np.angle(T_fit)), color=C_MUTED, lw=2.8,
        alpha=0.5, label="measured")
ax.plot(f_fit / 1e12, np.unwrap(np.angle(T_model_fit)), color=C_SI, lw=1.2,
        label="TMM fit")
ax.set_xlabel("Frequency  (THz)")
ax.set_ylabel("arg T  (rad)")
ax.set_title("Transmission phase", fontsize=9.5)
ax.legend(loc="lower left", ncol=2)
_despine(ax)

ax = axes[3]
ax.plot(f_fit / 1e12, np.abs(T_model_fit - T_fit), color=C_QZ, lw=1.1)
ax.set_xlabel("Frequency  (THz)")
ax.set_ylabel("|residual|")
ax.set_title(f"Fit residual   (RMS = {rms_res:.4f}, i.e. {100*_rel_res:.1f} % of "
             f"mean |T|)", fontsize=9.5)
_despine(ax)

fig.savefig(os.path.join(HERE, "fig4_tmm_fit.png"))
plt.close(fig)


# ==============================================================================
# STEP 9 -- NUMERICAL SUMMARY
# ==============================================================================
def band_avg(x, xf, lo=0.4e12, hi=1.6e12):
    """Mean of a spectrum over a quiet sub-band, ignoring NaNs."""
    sel = (xf >= lo) & (xf <= hi) & np.isfinite(x)
    return float(np.mean(x[sel])) if sel.any() else float("nan")


L = []
A = L.append
A("=" * 78)
A("THz-TDS ANALYSIS -- Si wafer on quartz window")
A("=" * 78)
A("")
A("ACQUISITION")
A(f"  samples / scan          : {N}")
A(f"  sampling interval       : {dt*1e15:.1f} fs   (Nyquist {1/(2*dt)/1e12:.1f} THz)")
A(f"  record length           : {(t[-1]-t[0])*1e12:.1f} ps  "
  f"(raw resolution {1/(t[-1]-t[0])/1e9:.1f} GHz)")
A(f"  zero-padded FFT length  : {ZERO_PAD}  -> grid {f[1]/1e9:.2f} GHz")
A(f"  analysis band           : {F_MIN/1e12:.2f} - {F_MAX/1e12:.2f} THz")
A("")
A("(A) ECHO TIMING  -- most robust quantities here; no model, no calibration")
A(f"  main-pulse delay (stack vs air) : {dt_main*1e12:8.3f} ps")
A(f"  Si     round-trip delay         : {tau_si*1e12:8.3f} ps"
  f"  ->  n*d = {nd_si*1e6:8.2f} um")
A(f"  quartz round-trip delay         : {tau_qz*1e12:8.3f} ps"
  f"  ->  n*d = {nd_qz*1e6:8.2f} um")
A(f"  echo polarities (Si {sgn_si:+.0f}, quartz {sgn_qz:+.0f}) vs Fresnel "
  f"prediction (+1, -1): {'OK' if polarity_ok else 'MISMATCH'}")
A(f"  closure, echoes vs bulk delay   : {closure*100:5.2f} % "
  f"{'OK' if closure < 0.03 else 'SUSPECT'}")
A("")
A("(C) GLOBAL TRANSFER-MATRIX FIT  -- primary source of absolute n, k and d")
A(f"  d_Si                    : {d_si_fit*1e6:8.2f} um"
  f"   (nominal {D_SI_NOM*1e6:.0f}, {100*(d_si_fit/D_SI_NOM-1):+.1f} %)")
A(f"  d_quartz                : {d_qz_fit*1e6:8.2f} um"
  f"   (nominal {D_QZ_NOM*1e6:.0f}, {100*(d_qz_fit/D_QZ_NOM-1):+.1f} %)")
A(f"  RMS residual in T       : {rms_res:.4f}   (|T| ~ {np.mean(np.abs(T_fit)):.2f})")
A(f"  optimiser               : "
  f"{'converged' if fit.success else 'DID NOT CONVERGE'} ({fit.nfev} evals)")
A("")
A("  thickness cross-check -- echo timing divided by the fitted index:")
A(f"    d_Si     = n*d / n_Si  = {d_si_echo*1e6:8.2f} um"
  f"   ({100*(d_si_echo/d_si_fit-1):+.2f} % vs fit)")
A(f"    d_quartz = n*d / n_qz  = {d_qz_echo*1e6:8.2f} um"
  f"   ({100*(d_qz_echo/d_qz_fit-1):+.2f} % vs fit)")
A("")
A("OPTICAL CONSTANTS, averaged over 0.4 - 1.6 THz")
A("                            TMM fit     echo phase")
A(f"  n_Si                   : {band_avg(n_si_tmm, f_fit):9.4f}  "
  f"{band_avg(n_si_echo, f):12.4f}")
A(f"  n_quartz               : {band_avg(n_qz_tmm, f_fit):9.4f}  "
  f"{band_avg(n_qz_echo, f):12.4f}")
A(f"  k_Si                   : {band_avg(k_si_tmm, f_fit):9.5f}")
A(f"  k_quartz               : {band_avg(k_qz_tmm, f_fit):9.5f}")
A(f"  alpha_Si     [1/cm]    : {band_avg(alpha_si_tmm, f_fit):9.2f}"
  f"   <-- below the {alpha_min_si:.1f} cm^-1 detection limit: UPPER BOUND ONLY")
A(f"  alpha_quartz [1/cm]    : {band_avg(alpha_qz_tmm, f_fit):9.2f}"
  f"   (detection limit {alpha_min_qz:.2f} cm^-1)")
A("")
A("  Absorption sensitivity floor = 2 * (RMS residual / |T|) / d. Fitted alpha")
A("  below it means 'not resolvable', not 'zero'.")
A("")
A("DIAGNOSTIC -- echo magnitudes (reported, NOT inverted for n; see docstring)")
A(f"  measured |echo/main|, Si     : {abs(rho_si_meas):.4f}  "
  f"(sign {sgn_si:+.0f})")
A(f"  measured |echo/main|, quartz : {abs(rho_qz_meas):.4f}  "
  f"(sign {sgn_qz:+.0f})")
r_si_pred = ((n_si_ref - 1) / (n_si_ref + 1)) * \
            ((n_si_ref - n_qz_ref) / (n_si_ref + n_qz_ref))
r_qz_pred = ((n_qz_ref - 1) / (n_qz_ref + 1)) * \
            ((n_qz_ref - n_si_ref) / (n_qz_ref + n_si_ref))
A(f"  Fresnel product from fitted n, Si     : {r_si_pred:+.4f}")
A(f"  Fresnel product from fitted n, quartz : {r_qz_pred:+.4f}")
A("  (measured magnitudes sit below these because the echo is also attenuated")
A("   by a double pass through the layer -- that is why they are not inverted.)")
A("")
A("CAVEATS")
A("  * The reference scan carries its own etalon at +1.5 and +3.6 ps (~3 % of")
A("    peak) from optics in the beam path. It sits inside the main-pulse gate")
A("    and is the dominant systematic on the extracted absorption.")
A("  * |echo/main| is non-monotonic in frequency (rises between 0.6 and 1.0 THz),")
A("    which is unphysical for pure absorption. Until that is understood, treat")
A("    alpha as indicative; n is unaffected, being phase-derived.")
A("  * The record is truncated while ~1 % ringing is still present, setting a")
A("    floor on fringe-fitting accuracy and hence on the thicknesses.")
A("  * Si absorption is near the limit this dynamic range can resolve; small")
A("    alpha_Si is an upper bound, not a measurement.")
A("=" * 78)

report = "\n".join(L)
print(report)
with open(os.path.join(HERE, "thz_results.txt"), "w") as fh:
    fh.write(report + "\n")

print("\nFigures written:")
for p in ("fig1_time_domain.png", "fig2_spectra_transmission.png",
          "fig3_optical_constants.png", "fig4_tmm_fit.png"):
    print("  ", p)
