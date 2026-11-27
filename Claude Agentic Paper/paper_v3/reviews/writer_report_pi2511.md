# Writer report: PI re-measurement (2026-11-25 bracket) in the paper

Compile (`bash briefs/compile.sh`): **8 pages**, 0 errors, 0 undefined references/citations, 0 overfull boxes. Abstract 198 words; its number count is unchanged (the new sentence has no numbers).
Backup of the pre-edit manuscript: `/tmp/main_before_pi2511.tex`.

## What changed (main.tex line numbers after the edit)

| line | where | change |
|---|---|---|
| 5 | header comment | number source for PI/Si is now `nb2 pi_2511` (the 2026-11-25 bracket) |
| 47 | abstract | "a time offset … must be fitted, so absolute thicknesses rest on time zero" → "the absolute thicknesses rest on time zero, whose drift a reference–sample–reference bracket measures" |
| 58 | Fig. 1 caption (b) | "fitted with a free time offset τ (σ = 1.35e-3)" → "fit without τ against the mean reference (MAE 1.12σ, σ = 1.34e-3); residual and ±σ band ×10" |
| 93 | Sec 4 PI paragraph | rewritten as the bracket story (see numbers below). Rendered length goes from 11 to 12 lines. The τ, a model is now defined in words, not by an equation, to save space |
| 94–95 | LaTeX comment | the old "AUTHOR NOTE … will be re-measured … (R27)" comment is removed. New record comment: the old pair's numbers are in results/nb2.json (loop1_checks; also pi_2511/comparison_old_vs_new) |
| 98–106 | Table tab:pi | rows: mean ref. (no τ), CR bound, ref. 1 only, ref. 2 only, independent (unchanged). Caption gives the times 14:28/14:48/15:06, 10 000 waveforms each, and the 0.15 % / 2.7 % n·d spread |
| 155 | Sec 6.2 (iv) | parenthetical "(we fit τ only where, as for PI, it brings a fixed-index substrate into agreement …)" → "(for PI the drift is measured instead, Section 4)". Not in the brief: the old text contradicted the new no-τ protocol |
| 206 | conclusions | "given time zero, the total thickness" → "given a bracketed time zero, the total thickness". The R33/R29 quartz sentence is unchanged |

## Deviation: τ, a row and columns dropped from Table 1
With the row "τ, a free (mean ref)" plus the τ and a columns, the paper compiles to **9 pages** (tested). Removing them restores 8. The τ, a result stays in the prose: "τ of −1 to −5 fs, inside the drift, the gain coming from a = 1.014". Its full row is in `pi_2511/fits/all_free_tau_amp|mean` (n_PI 1.769, κ_PI 0.014, d_PI 25.63, n_Si 3.407, d_Si 382.6, n·d 45.34, τ −1.2, a 1.014, rMSE 1.90σ) if the director wants it back and room is found.
Some prose numbers were also left out for space: the ref1/ref2 n·d range 44.7–47.2 µm is shown only in the table, and the 10–15 µm PI-like threshold is gone from the text.

## Every number and its JSON key (all under results/nb2.json → pi_2511)

| number in paper | key | JSON value |
|---|---|---|
| 23 min (old ref→sample) | data_checks/sample_delays/old_pair/minutes_ref_to_sample | 22.88 |
| −65 to −70 fs (old τ) | data_checks/old_pair_tau_fs | −65.44 / −70.14 |
| 38 min | data_checks/time_zero_drift/minutes_ref1_to_ref2 | 38.43 |
| +40 fs drift | data_checks/time_zero_drift/ref2_minus_ref1_phase_slope_fs | 40.03 (peak 40.73) |
| 1.0–1.5 % amplitude | data_checks/amplitude_ratio_ref2_over_ref1 (rms / peak) | 1.0100 / 1.0151 |
| n_Si 0.46 % below | criteria/n_Si_without_tau/new_minus_lit_percent | −0.459 |
| −6.6 fs | criteria/n_Si_without_tau/equivalent_time_zero_error_fs | −6.58 |
| literature value between single-reference fits | criteria/n_Si_without_tau/literature_inside_ref1_ref2_bracket | true |
| both κ at bound | fits/all_free\|mean/at_bound | true, true |
| τ −1 to −5 fs | fits/all_free_tau_amp\|mean, nSi_fixed_tau_amp\|mean tau_fs | −1.22 / −5.12 |
| a = 1.014 | fits/all_free_tau_amp\|mean/best/amplitude | 1.0144 |
| 38.9 fs | criteria/tau/nSi_fixed_tau_ref1_minus_ref2_fs | 38.95 |
| d_PI 7 % above profilometry | fig2b_variant/d_PI_minus_profilometer_um ÷ 24.7 | 1.722 µm → 6.97 % |
| CR bound 29 nm | fig2b_variant/cramer_rao/d_PI_nm (nd_PI_nm) | 29.45 (28.77) |
| 0.16 % PI-like sweep, 10 fs ≈ 3 µm | unchanged from loop-1 text (C11c, R27) | — |
| Fig 1(b): MAE 1.12σ, σ 1.34e-3 | fig2b_variant/mae_over_sigma, sigma | 1.1201, 1.3445e-3 |
| Table, mean ref.: 1.751, at bound, 26.42, 3.402, at bound, 383.2, 46.28, 2.23 | fits/all_free\|mean/best (n_PI, at_bound, d_PI_um, n_Si, d_Si_um, nd_PI_um, rmse_over_sigma) | 1.75147, –, 26.422, 3.40182, –, 383.195, 46.278, 2.2253 |
| Table, CR: 0.0014, 0.029, 0.00017, 0.028, 0.029 | fits/all_free\|mean/cramer_rao (n_PI, d_PI_nm, n_Si, d_Si_nm, nd_PI_nm) | 1.424e-3, 29.45, 1.721e-4, 27.52, 28.77 |
| Table, ref 1: 1.764, 7e-6, 26.74, 3.451, at bound, 377.8, 47.16, 3.22 | fits/all_free\|ref1/best | 1.76390, 7.4e-6, 26.737, 3.45056, κ_Si at bound, 377.763, 47.161, 3.2199 |
| Table, ref 2: 1.752, 0.018, 25.53, 3.355, at bound, 388.6, 44.74, 2.02 | fits/all_free\|ref2/best | 1.75205, 0.01771, 25.535, 3.35513, κ_Si at bound, 388.573, 44.738, 2.0198 |
| 2 of 8 starts, 0.15 %, 2.7 % | fig2b_variant/n_starts_agreeing_within_0p5pct, nd_PI_spread_starts_within_0p5pct, nd_PI_spread_all_8_starts | 2, 0.151 %, 2.653 % |
| 14:28 / 14:48 / 15:06, 10 000 waveforms | timestamps; data_checks/noise/*/n_waveforms | 14:28:22 / 14:48:33 / 15:06:48; 10000 |
| independent 1.726, 24.7, 3.4175, 42.63 | unchanged | — |

Note: the 38.9 fs comes from the n_Si-fixed τ fits, and the text says so ("with n_Si fixed"). The all-free single-reference τ values (+1.7 / −40.6 fs) differ by 42.3 fs.
