# Coder report: polyimide re-measurement of 2026-11-25 (nb2 Part C)

Part C was added to `notebooks/nb2_thin_film_limits_and_PI.ipynb`, after the loop-1 checks. Parts A and B and the loop-1 checks are
unchanged and reproduce bit-for-bit (Part B n_Si 3.25921, loop-1 fig2b_variant 1.6956σ, C11 values identical). The notebook was
re-executed in place with nbconvert (`RUN_FULL = False`, `NB_THREADS=2`, about 18 min). Every number below is in
`results/nb2.json → pi_2511`. σ always follows the R12 convention: the pre-pulse std over the first 10 ps of the *sample* trace,
1.345e-3 for the new sample and 1.349e-3 for the old one. Sign convention: a positive delay or τ means the sample is later.

## 1. Data checks (`pi_2511/data_checks`)

| trace | time | σ (0–10 ps) | white level | main peak | t_peak (ps) |
|---|---|---|---|---|---|
| ref1 | 14:28 | 6.33e-3 | 1.5e-4 | −3.297 | 12.3618 |
| sample | 14:48 | 1.345e-3 | 1.6e-5 | −2.488 | 15.5148 |
| ref2 | 15:06 | 6.34e-3 | 1.6e-4 | −3.347 | 12.4025 |
| old ref (11-12) | 13:13 | 6.50e-3 | 1.1e-4 | −3.304 | 12.2870 |
| old sample (11-12) | 13:36 | 1.349e-3 | 1.9e-5 | −2.487 | 15.3607 |

- **Time-zero drift, ref2 − ref1 (38.4 min):** +40.0 fs from the cross-spectral phase slope (0.2–2 THz), +40.7 fs from the main
  peak, +39.7 fs from the cross-correlation. That is 1.04 fs/min, with the later scan arriving later.
- **Amplitude ratio ref2/ref1:** 1.015 (peak), 1.010 (rms), 1.012 (|X| over 0.2–2 THz).
- **Sample delay (phase slope / main peak):** vs ref1 3159.1 / 3153.1 fs, vs ref2 3119.1 / 3112.3 fs, vs mean 3138.8 / 3133.5 fs.
  The old pair gave 3080.2 / 3073.7 fs. The new delay against the mean reference is therefore **+58.6 / +59.8 fs longer** than the
  old one, which matches the −65…−70 fs τ that the old pair needed.
- If the drift is linear in time, the sample's timing error is expected to be +21.0 fs against ref1, −19.0 fs against ref2 and
  +1.0 fs against the mean.

## 2–3. Fits (same 8 starts and protocol as Part B; `pi_2511/fits`)

`*` marks κ at the hinge bound (no CR entry). "agree" counts the starts within 0.5 % of the best rMSE. "n·d spread" is the spread of
n·d_PI over the agreeing starts.

| variant | reference | rMSE/σ | n_PI | κ_PI | d_PI (µm) | n·d_PI (µm) | n_Si | κ_Si | d_Si (µm) | τ (fs) | a | agree | n·d spread |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all free, no τ | ref1 | 3.220 | 1.7639 | 7e-6 | 26.74 | 47.16 | 3.4506 | 0* | 377.76 | — | — | 2/8 | 0.09 % |
| all free, no τ | ref2 | 2.020 | 1.7520 | 0.0177 | 25.53 | 44.74 | 3.3551 | 0* | 388.57 | — | — | 8/8 | 0.41 % |
| **all free, no τ** | **mean** | **2.225** | 1.7515 | 0* | **26.42** | 46.28 | **3.4018** | 0* | 383.19 | — | — | 2/8 | 0.15 % |
| n_Si fixed, no τ | mean | 3.030 | 1.7907 | 0* | 25.35 | 45.40 | 3.4175 | 0* | 380.84 | — | — | 3/8 | 0.24 % |
| all free + τ, a | ref1 | 2.385 | 1.7925 | 0.0263 | 25.22 | 45.20 | 3.4502 | 6.5e-4 | 377.83 | +1.7 | 1.029 | 5/8 | 0.58 % |
| all free + τ, a | ref2 | 1.681 | 1.8063 | 0.0455 | 24.61 | 44.45 | 3.4590 | 1.1e-4 | 376.87 | −40.6 | 1.025 | 6/8 | 0.05 % |
| all free + τ, a | mean | 1.900 | 1.7688 | 0.0138 | 25.63 | 45.34 | 3.4070 | 0* | 382.61 | −1.2 | 1.014 | 4/8 | 0.25 % |
| n_Si fixed + τ, a | ref1 | 2.409 | 1.7779 | 0.0207 | 25.42 | 45.20 | 3.4175 | 5.1e-4 | 381.46 | +14.3 | 1.022 | 7/8 | 0.00 % |
| n_Si fixed + τ, a | ref2 | 1.735 | 1.7863 | 0.0375 | 24.88 | 44.43 | 3.4175 | 0* | 381.48 | −24.7 | 1.016 | 8/8 | 0.08 % |
| n_Si fixed + τ, a | mean | 1.905 | 1.7710 | 0.0162 | 25.55 | 45.25 | 3.4175 | 0* | 381.44 | −5.1 | 1.017 | 6/8 | 0.58 % |

**Is n_Si physical without τ?** Nearly. Against the mean reference n_Si = 3.4018, which is −0.46 % from 3.4175, and all 8 starts lie
between 3.3998 and 3.4022. The old pair gave 3.2592 (−4.6 %). The two single references give 3.4506 and 3.3551, so the literature
value lies inside the range set by the timing bracket. From those two fits, ∂n_Si/∂τ = 2.38e-3 per fs. The remaining −0.0157
therefore corresponds to a −6.6 fs time-zero error. By the same scale, the old pair's n_Si corresponds to −66 fs, which is the old τ.
Under the noise-only CR bound (1.7e-4) the deviation is still formally significant. The actual time-zero uncertainty of the bracket
(±20 fs) is far larger than that bound.

**τ:** with n_Si fixed, τ = +14.3 fs against ref1, −24.7 fs against ref2 and −5.1 fs against the mean. The ref1 − ref2 difference is
38.9 fs, matching the measured drift of 40.0 fs. All three values sit about 6 fs below the linear-drift expectation (+21 / −19 / +1 fs).
Against the mean reference |τ| ≤ 5.1 fs, which is within the measured drift and within half of it. The old τ of −65 / −70 fs exceeds
the whole 40 fs drift measured over 38 min on 11-25, and the old reference and sample were 23 min apart. With the mean reference the
τ, a fit gains 15 % in rMSE (2.225 → 1.900σ). That gain comes from the amplitude factor a = 1.014–1.017, which matches the ref1→ref2
power change of 1.0–1.5 %, not from τ.

**CR bounds** (σ = 1.345e-3, exact model; no entry for fixed or at-bound parameters):
- No-τ fit against the mean: n_PI 1.4e-3, d_PI 29 nm, n·d_PI 29 nm, n_Si 1.7e-4, d_Si 28 nm, d_PI + d_Si 49 nm.
- All free + τ, a against the mean: n_PI 1.9e-3, κ_PI 6.9e-4, d_PI 39 nm, n·d_PI 37 nm, n_Si 2.1e-3, d_Si 240 nm, τ 0.83 fs.
- n_Si fixed + τ, a against the mean: d_PI 35 nm, d_Si 12 nm, τ 0.07 fs.

## 4. Old vs new, and the verdict (`pi_2511/comparison_old_vs_new`, `criteria`, `verdict`)

| | old 11-12, no τ | **new 11-25, no τ, mean ref** | old 11-12, τ, a (all free / n_Si fixed) | new 11-25, τ, a, mean ref (all free / n_Si fixed) |
|---|---|---|---|---|
| rMSE/σ | 2.687 | **2.225** | 1.696 / 1.704 | 1.900 / 1.905 |
| n_Si | 3.2592 (−4.6 %) | **3.4018 (−0.46 %)** | 3.4297 / fixed | 3.4070 / fixed |
| d_PI (µm), profilometer 24.7 | 25.84 | 26.42 | 25.24 / 25.25 | 25.63 / 25.55 |
| n·d_PI (µm) | 44.26 | 46.28 | 45.01 / 44.92 | 45.34 / 45.25 |
| τ (fs) | — | — | −70.1 / −65.4 | −1.2 / −5.1 |
| starts within 0.5 %, n·d spread | 6/8, 0.20 % | 2/8, 0.15 % | 2/8 0.22 % / 6/8 0.29 % | 4/8 0.25 % / 6/8 0.58 % |

**Verdict: better**, because the new measurement fixes the problem it was taken for, time zero:
- The drift is now measured (+40 fs in 38 min) and does not need to be fitted.
- Without τ, the paper protocol gives n_Si = 3.402 instead of 3.259, and its misfit falls from 2.69σ to 2.23σ.
- When τ is free, it comes out at −1 to −5 fs against the mean reference, inside the measured drift.

What did not improve:
- d_PI moves further from the profilometer value: +1.7 µm (+7 %) without τ, and +0.85–0.93 µm with τ, compared with +1.1 µm and
  +0.55 µm before.
- n·d_PI without τ rises to 46.3 µm and depends on which reference is used (44.7–47.2 µm). With τ, a it is 44.4–45.3 µm.
- Only 2 of 8 starts agree in the no-τ fit, against 6 of 8 before.
- Both κ go to the bound in the no-τ fit against the mean reference. Without the amplitude factor the model cannot place any absorption.
- The best τ, a misfit (1.90σ) is higher than the old one (1.70σ), although that old value came with an unexplained −70 fs offset.

## 5. Fig. 2(b) (`pi_2511/fig2b_variant`)

n_Si passes the criterion: the literature value lies inside the ref1/ref2 bracket and is less than 1 % away. Panel (b) therefore shows
the **paper-protocol fit without τ against the mean reference** (legend "fit"). Its numbers:
- n_PI 1.7515, κ_PI at bound (no CR entry), d_PI 26.42 µm, n·d_PI 46.28 µm.
- n_Si 3.4018, κ_Si at bound, d_Si 383.19 µm; τ and a not fitted (a = 1).
- rMSE 2.225σ, MAE 1.12σ.
- CR bounds as in section 3.
- n·d spread: 2.65 % over all 8 starts, 0.15 % over the 2 agreeing starts. Film echo marked at +0.31 ps.

The smallest-|τ| alternative (all free + τ, a against the mean; τ = −1.2 fs, rMSE 1.90σ) is stored as
`pi_2511/fig2b_alternative_smallest_abs_tau` and can be selected by changing `FIG2B_C` in the Part C figure cell.

Panel (a) and the layout are unchanged: the drawing code is now `draw_fig2`, with an identical export (13.05 cm wide, smallest
printed font 7.13 pt at 13.3 cm). The loop-1 version is kept as `figs_old/fig2_thin_film_loop1.*`.

**Needed in the paper:**
- The Fig. 2 caption ("fitted with a free time offset τ") must change.
- Table `tab:pi` and the Sec. 4 PI paragraph still quote the 11-12 numbers. The τ narrative should become: drift measured with a
  bracket, τ not needed.
