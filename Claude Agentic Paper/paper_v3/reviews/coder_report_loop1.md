# Coder report — loop 1 (2026-11-24)

All checks C1–C14 of `reviews/loop1_domain.md` §5 were run **live inside the notebooks** (appended "Loop-1 checks" part of
each notebook; the long loop-0 fits are reloaded from `results/*_cache.npz`, `RUN_FULL = False`). Every number below is in
`results/nbX.json → loop1_checks` (nb6: top level of `results/nb6.json`). New notebook `nb6_bilayer_model_variants.ipynb`
holds C3 and C4 (they need their own models; nb3 was already long). New helper module `code/loop1_tools.py`;
`code/paperstyle.py` gained a font/width check (`font_report`, automatic shrink in `savefig`). All six notebooks were
re-executed in place with nbconvert (outputs embedded); `notebooks/README.md` updated.

σ below is always the R12 convention (pre-pulse std, first 10 ps) unless stated otherwise.

---

## C14 — noise table of the six traces (nb2 `C14_noise_table`, quartz/Si part also nb3 `C14_noise_quartzSi`)

| trace | std 0–10 ps (= σ, R12) | std 0–5 / 5–10 ps | white (first-diff., quietest 2-ps seg.) | σ / white | max \|segment mean\| | ACF 1/e (ps) |
|---|---|---|---|---|---|---|
| quartz/Si reference | 1.58e-3 | 1.11e-3 / 1.93e-3 | 4.5e-5 | 35 | 5.8e-4 | 0.26 |
| quartz/Si sample | 2.90e-4 | 2.23e-4 / 2.18e-4 | 4.5e-5 | 6.4 | 5.5e-4 | 0.76 |
| PI reference | 6.50e-3 | 2.12e-3 / 8.94e-3 | 1.1e-4 | 57 | 1.9e-3 | 0.26 |
| PI sample | 1.35e-3 | 1.03e-3 / 1.60e-3 | 1.9e-5 | 72 | 5.7e-4 | 0.26 |
| run-3 reference | 2.40e-2 | 1.82e-3 / 3.38e-2 | 2.8e-5 | 848 | 1.5e-2 | 0.30 |
| run-3 sample | 9.04e-4 | 4.72e-4 / 1.19e-3 | 2.7e-5 | 34 | 3.5e-4 | 0.23 |

**Meaning.** The R12 σ is dominated by bursts and offsets; the white component is 6–850× smaller (6–72× for the sample
traces). Confirms domain 2.8 / R19: σ is a bookkeeping convention, not a white-noise level.

## C1 — PI/Si model variants (nb2 `C1_PI_model_variants`; same 8 starts, coarse-to-fine L-BFGS, rMSE + hinge)

| variant | n_PI | κ_PI | d_PI (µm) | n·d_PI (µm) | n_Si | κ_Si | d_Si (µm) | τ (fs) | a | rMSE/σ | vs loop 0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| loop 0, all free | 1.7125 | 0 (bound) | 25.84 | 44.26 | 3.2592 | 5.8e-4 | 400.19 | — | — | 2.687 | — |
| (a) n_Si = 3.4175 | 2.998 | 2.2e-3 | 117.0 | 350.9 | 3.4175 | 5e-5 | 283.6 | — | — | 18.06 | +572 % |
| (b) n_Si = 3.4175, τ, a free | 1.7792 | 0.0262 | **25.25** | 44.92 | 3.4175 | 0 (bound) | 381.60 | **−65.4** | 1.022 | **1.704** | −36.6 % |
| (c) all free + τ, a | 1.7836 | 0.0286 | 25.24 | 45.01 | **3.4297** | 0 (bound) | 380.24 | −70.1 | 1.024 | 1.696 | −36.9 % |
| (d) as (b), d_PI = 24.7 fixed | 1.8013 | 0.0305 | 24.70 | 44.49 | 3.4175 | 0 | 381.51 | −65.0 | 1.024 | 1.749 | −34.9 % |

6/8 starts agree within 0.5 % in (b), 2/8 in (c), 5/8 in (d).
**C6 (PI part):** CR at the loop-0 point without/with τ, a: d_PI 62/64 nm, d_Si 30/251 nm, d_PI + d_Si 73/273 nm,
n_PI·d_PI 52/52 nm, n_Si 1.7e-4/2.0e-3, τ ± 0.9 fs.
**Meaning.** Yes: a free time-zero offset of −65 to −70 fs (sample/reference) removes the impossible n_Si. With n_Si free it
goes to 3.430 (+0.36 %), with it fixed d_PI = 25.25 µm (+2.2 % vs profilometer), and the misfit falls by 37 %. Forcing
d_PI = 24.7 µm costs only +2.6 % in rMSE over (b). The absorption moves from Si to the film (κ_PI ≈ 0.026–0.03, κ_Si at the
bound), which is the physically plausible assignment. Without τ, fixing n_Si is catastrophic (a), so the loop-0 n_Si = 3.26
was absorbing a ~65 fs timing error. Suggested framing: "the PI extraction exposes a time-zero systematic; once τ is free,
d_PI = 25.2 µm and n_Si = 3.43". n·d_PI moves 44.3 → 44.9 µm. The visible-index × profilometer product 1.726 × 24.7 = 42.63 µm
is now in nb2.json (`independent_nd_visible_index_times_profilometer_um`).

## C11 — thin-film sweep (nb2 `C11_sweep`)

**(a) CR bounds at the truth (sweep geometry, σ = 1.58e-3).** Reference pulse FWHM: field peak 0.355 ps, intensity envelope
0.511 ps (film round trip at 10 µm: 0.24 ps).

| d (µm) | 5 | 7 | 10 | 15 | 20 | 30 | 50 | 70 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| CR n (%) | 2.21 | 0.90 | 0.34 | 0.105 | 0.055 | 0.044 | 0.047 | 0.049 | 0.049 |
| CR κ (%) | 46 | 28 | 16 | 10 | 7.7 | 5.8 | 4.0 | 3.1 | 2.3 |
| CR d (%) | 5.9 | 2.4 | 0.94 | 0.30 | 0.13 | 0.068 | 0.053 | 0.053 | 0.051 |
| loop-0 median\|err\| / (0.674 CR), d | 1.22 | 1.06 | 1.55 | 1.81 | 0.93 | 1.00 | 0.33 | 0.51 | 1.21 |

**(b) 10 more draws at 5–20 µm** (5 starts each; best start re-polished with 1500 full-band L-BFGS iterations).

| d (µm) | median\|err\| n / κ / d (%) | std(err)/CR n / κ / d | mean signed d err (%) ± s.e. | max polish shift of d | 13-draw median\|err\| d (%) |
|---|---|---|---|---|---|
| 5 | 0.64 / 20 / 1.35 | 0.77 / 0.82 / 0.80 | +1.3 ± 1.5 | 1.7e-7 µm | 2.87 |
| 7 | 0.64 / 11 / 1.64 | 0.99 / 1.04 / 0.98 | −1.1 ± 0.75 | 8e-9 µm | 1.66 |
| 10 | 0.30 / 17 / 0.86 | 0.94 / 1.50 / 0.98 | −0.21 ± 0.29 | 2e-4 µm | 0.88 |
| 15 | 0.093 / 11 / 0.28 | 1.26 / 1.19 / 1.14 | −0.13 ± 0.11 | 3e-4 µm | 0.31 |
| 20 | 0.060 / 6.5 / 0.12 | 1.12 / 1.06 / 1.15 | +0.05 ± 0.05 | 2e-8 µm | 0.12 |

**Meaning.** The scatter equals the CR bound to within ±25 % at every thickness (std/CR 0.77–1.26 for n and d), there is no
significant bias, and re-polishing moves nothing (≤ 3e-4 µm). So the sweep errors **are** noise-limited (white noise, exact
model). Domain 2.21's worry about optimiser-limited scatter is not borne out: at 5 µm the starts at equal loss spread by only
0.08 µm (median) after the coarse-to-fine schedule. The median d error crosses 1 % between 7 and 10 µm (13 draws: 1.66 % → 0.88 %),
so "the median error in d exceeds 1 % below 10 µm (round trip 0.24 ps, pulse FWHM 0.36 ps)" stands.

**(c) PI-like sweep** (film 1.72 / 0.01 on Si 3.4175 / 1e-4 / 400 µm, PI reference, σ = 1.35e-3, 5 draws × 5 starts).
Echo Fresnel product 0.087 vs 0.135 for the sweep.

| d (µm) | 10 | 15 | 20 | 25 | 30 | 40 |
|---|---|---|---|---|---|---|
| round trip (ps) | 0.115 | 0.172 | 0.229 | 0.287 | 0.344 | 0.459 |
| median\|err\| n / κ / d (%) | 1.2 / 28 / 1.85 | 0.35 / 30 / 0.33 | 0.20 / 14 / 0.28 | **0.043 / 17 / 0.16** | 0.020 / 11 / 0.073 | 0.044 / 7.9 / 0.057 |
| median\|err\| n·d (%) | 0.74 | 0.11 | 0.10 | 0.13 | 0.053 | 0.022 |
| CR n / κ / d (%) | 2.3 / 105 / 3.2 | 0.71 / 52 / 1.1 | 0.30 / 29 / 0.47 | 0.16 / 19 / 0.26 | 0.096 / 13 / 0.17 | 0.055 / 7.8 / 0.10 |

**Meaning.** In an exact-model world a 25-µm PI film on Si is recovered to 0.16 % in d and 0.04 % in n. It is **not** "at the
threshold": the PI-like 1 % crossing in d lies between 10 and 15 µm. In thickness terms that is slightly above the sweep's
7–10 µm (consistent in direction with D1), but in round-trip terms it is lower (≈ 0.13 ps vs 0.2 ps). The PI discrepancies are
systematic (C1), not a resolution limit. Fig. 2(a) shows the 25-µm PI-like median and range as hollow markers.

## C5 — |κ₂| against noise models (nb3 `C5_kappa2_vs_noise_models`)

| noise model | Adam best (κ₂ = −1.46e-3): \|κ₂\|/bound | (n₂ − 3.4175)/bound | rMSE optimum (κ₂ = −7.27e-3): \|κ₂\|/bound | (n₂ − 3.4175)/bound |
|---|---|---|---|---|
| canonical σ = 2.9e-4 (R12) | 77 | 108 | 372 | −425 |
| white σ from first differences (4.5e-5) | 493 | 693 | 2390 | −2729 |
| GLS, Toeplitz pre-pulse covariance | 91 | 89 | 447 | −350 |
| residual-RMS-scaled (RMS = 41.3σ / 38.7σ) | **1.9** | 2.6 | 9.6 | −11 |
| sandwich HAC, Bartlett 1 ps | **0.7** | 0.9 | 3.8 | −3.8 |
| sandwich HAC, Bartlett 3 ps | 0.7 | 0.7 | 3.6 | −3.3 |
| block bootstrap (1 ps blocks, 100×, rMSE) | — | — | **4.4** (std κ₂ 1.7e-3; κ₂ < 0 in 99 %) | −4.0 (std n₂ 0.028) |

The rMSE optimum (warm-started L-BFGS, 2173 evaluations, |∇| 4e-6) is n₁ 1.99678, κ₁ 0.00778, d₁ 977.709 µm, n₂ 3.30372,
κ₂ −0.00727, d₂ 398.136 µm, identical to nb4. Bootstrap std: d₁ 3.3 µm, d₂ 3.3 µm.
**Meaning.** "77×" holds only under the pre-pulse noise model; any model-error-aware uncertainty (residual-scaled, sandwich)
makes κ₂ < 0 at the Adam point insignificant (0.7–1.9 bounds), and only moderately significant at the least-squares point
(3.6–9.6). Step (iv) should rest on the residual RMS (41σ) and on the estimator dependence of κ₂ (−0.0073 LS, −0.0015 Adam,
+0.0044 L-BFGS on Eq. 4), not on 77. The same holds for n_Si.

## C6 — τ/a columns and the valley directions (nb3 `C6_tau_amp_and_valley_directions`)

| CR (σ = 2.9e-4) | d₁ | d₂ | d₁ + d₂ | n₁d₁ | n₂d₂ | n₁ | n₂ | \|κ₂\|/bound |
|---|---|---|---|---|---|---|---|---|
| without τ, a | 31.6 nm | 31.1 nm | 7.5 nm | 8.3 nm | 5.9 nm | 6.1e-5 | 2.8e-4 | 77 |
| with τ, a | 137 nm | 53 nm | 172 nm | 8.3 nm | 5.9 nm | 2.7e-4 | 4.8e-4 | 71 |

τ is bound to 0.57 fs. The thickness-transfer direction (d₁↑ d₂↓ at constant n₁d₁, n₂d₂) is eigen-direction s₄
(cos 0.985, stiffness 0.072 vs s₄ = 0.073). The κ split (κ₁d₁ + κ₂d₂ constant) is s₆ (cos 0.999). The (n, d) part of the
displacements Adam → softplus, Adam → rMSE and fit 1 → fit 2 lies along the transfer direction with cos 0.9999, 0.9997, 0.997
(Δd₁ = +2.4 / −16.5 / +0.12 µm, Δd₂ equal and opposite). **Meaning:** this confirms R17 quantitatively. Freeing τ leaves n·d
untouched but makes d₁ + d₂ 23× less certain, i.e. absolute thicknesses rest on time-zero.
**Invariants over all 16 end points** (nb4 `end_points_timing_invariants`): d₁ + d₂ 1375.5–1375.8 µm (0.026 %), Σ(n−1)d 0.014 %,
n₁d₁ 0.031 %, n₂d₂ 0.057 %, κ₁d₁ + κ₂d₂ 4.25–4.90 µm (14 %, 4.67–4.90 without the unconverged log-cosh start 1); d₁ 977.7–1000.0,
d₂ 375.7–398.1 µm, n₂ 3.304–3.499, κ₂ −0.0073 … +0.0044.

**Singular-value mini-table** (nb3 `singular_value_table`): singular values s_k of J/√L (√eig G), J = ∂y/∂x in
nominal-scaled coordinates x_i = Δp_i/p_i^nom:

| k | s_k | dominant combination |
|---|---|---|
| 1 | 8.665 | +0.69 n₁ +0.49 n₂ +0.39 d₂ +0.36 d₁ |
| 2 | 1.541 | +0.59 d₂ +0.51 n₂ −0.50 n₁ −0.37 d₁ |
| 3 | 0.6895 | +0.81 d₁ −0.45 n₁ +0.31 d₂ |
| 4 | 0.0733 | +0.66 n₂ −0.63 d₂ +0.27 d₁ (= thickness transfer) |
| 5 | 0.0192 | +0.93 κ₁ +0.32 κ₂ |
| 6 | 0.00159 | +0.95 κ₂ −0.32 κ₁ (= κ split) |

## C8 — Adam to 6×10⁴ iterations from start 1 (nb3 `C8_adam_to_60000`)

| iteration | 19 999 (loop 0) | 2×10⁴ | 3×10⁴ | 4×10⁴ | 5×10⁴ | 6×10⁴ |
|---|---|---|---|---|---|---|
| loss | 4.45223e-3 | 4.45220e-3 | 4.45165e-3 | 4.45173e-3 | 4.45177e-3 | 4.45305e-3 |
| κ₂ | −1.444e-3 | −1.445e-3 | −1.482e-3 | −1.475e-3 | −1.476e-3 | −1.480e-3 |
| d₁ / d₂ (µm) | 994.28 / 381.41 | 994.28 / 381.41 | 994.20 / 381.50 | 994.21 / 381.51 | 994.22 / 381.50 | 994.19 / 381.49 |
| \|∇_u L\| | 3.5e-2 | 3.5e-2 | 3.7e-3 | 3.9e-3 | 9.2e-3 | 8.9e-2 |

Iterate 19 999 reproduces loop-0 fit 1 exactly. **Meaning:** after 2×10⁴ iterations Adam hovers on a plateau. The loss is
within 0.013 % and κ₂ within 3 % of their later values, and the loss fluctuates (the 6×10⁴ value is the highest). It does not
converge to the least-squares point (different loss). "Budget iterate on a plateau; follows the same valley" is the accurate
statement.

## C7 — Adam learning rate (nb4 `C7_adam_learning_rate`)

L-BFGS reaches its stall level (4.6172–4.6178e-3) after 123 / 152 / 154 evaluations (first iterate within 1e-5); total
evaluations including the confirming restarts are 250 / 272 / 299.

| η | Adam iterations to the L-BFGS stall level (starts 1/2/3) | ratio to L-BFGS evals (to stall / total) | iterations to 4.5e-3 | diverged |
|---|---|---|---|---|
| 1e-4 (nb3) | 5223 / 5929 / 4998 | 33–43 / 17–22 | 9688 / 10449 / 9019 | no |
| 3e-4 | 2218 / 2494 / 2058 | 13–18 / 6.9–9.2 | 4091 / 4243 / 3582 | no |
| 1e-3 | 838 / 945 / 812 | 5.3–6.8 / 2.7–3.5 | 1480 / 1505 / 1269 | no |
| 3e-3 | 357 / 437 / 460 | 2.9–3.0 / 1.4–1.6 | 621 / 662 / 622 | no |

**Meaning:** "twenty times cheaper" is specific to η = 1e-4. At η = 3e-3 Adam needs only ~1.5× L-BFGS's evaluations and then
continues *below* the L-BFGS stall: it reaches 4.5e-3 in ~620 iterations, a level L-BFGS never reaches on this loss. One
Adam iteration = one closure (forward + backward) throughout.

## C9 — the kink (nb4 `C9_kink_probe`)

At the L-BFGS end point (start 1, loss 4.617188e-3, |∇| 3.65e-2):

| direction | ∇L·p | one-sided FD slope, t = 1e-9 … 1e-5 | residual samples changing sign (t ≤ 1e-6) |
|---|---|---|---|
| −∇L | −3.65e-2 | **+2.3e-2 … +5.2e-2** | 2 |
| last L-BFGS step | +3.43e-2 | +3.43e-2 … +3.8e-2 | 0 |
| toward Adam end point | −2.46e-4 | −2.46e-4 (all t) | 0 |

Gradient norms: rMSE end points 1.0–2.7e-6 (vs 2.6–6.0 at the starts; d₁ 977.7087–977.7089 µm); log-cosh 0.9e-7–1.1e-5
(start 1 capped at 400 iterations, d₁ 987.70); plain MAE 1.7–4.5e-2 (stalled like Eq. 4). **Meaning:** a textbook
non-differentiable point. Along −∇L the loss rises at once (slope +0.023 instead of −0.036, from residual samples crossing
zero), but a smooth descent direction to Adam's point exists (slope −2.5e-4). The rMSE "convergence" is real: tiny gradient,
three starts agree to 0.2 nm.

## C10 — L-BFGS iterates and the Adam planes (nb4 `C10_out_of_plane`)

| start | share of \|Δx\|² in own Phase-I plane | end point out of / in that plane | share in the Phase-II plane | end out of / in Phase-II plane |
|---|---|---|---|---|
| 1 | 12 % | 0.14 | 99.8 % | 0.050 / 1.13 |
| 2 | 55 % | 0.15 | 97.2 % | 0.051 / 1.13 |
| 3 | 25 % | 0.14 | 97.6 % | 0.050 / 1.13 |

**Meaning — note for the director:** in the coordinates of the shared Phase-II plane (nominal-scaled), the L-BFGS
displacement lies 97–99.8 % *in* that plane. It is dominated by κ₁ 0.0071 → 0.0031, which is the plane's v axis. So R21's
reason for dropping the Phase-II panel ("iterates lie off that plane") does not hold as stated. Conversely, only 12–55 % of
each L-BFGS path lies in its own Phase-I plane, so the Phase-I panels of Fig. 5 show projections. I followed R21 (panel
dropped); the caption should say "projections; 45–88 % of the L-BFGS displacement lies out of the Phase-I plane".

## C3 — six knots on the Sec. 5 bilayer, λ_κ = 0 (nb6 `C3_six_knots_bilayer`)

| | stage 1 / tied (flat) | stage 2, six knots |
|---|---|---|
| (i) Adam protocol (MAE + log-cosh) | κ_Si −1.83e-3, κ_q 5.41e-3, d 993.1 / 382.7 µm, loss 4.4535e-3 | Si κ knots −0.0028, −0.0029, −0.0022, −0.0036, −0.0026, +0.0006 (97 % of band < 0); Si n 3.413–3.430; d 991.3 / 385.1 µm; loss 2.7893e-3 (−37 %) |
| (ii) converged rMSE (L-BFGS) | κ_Si −7.27e-3, n_Si 3.304, rMSE 1.124e-2 | Si κ knots −0.0060 … −0.0030 (all < 0); Si n 3.392–3.402; d 989.3 / 387.3 µm; rMSE 4.840e-3 (−57 %, \|∇\| 1e-7) |

Protocol deviation: at lr 1e-4 Adam moves κ by ≤ 5e-7 per iteration, so each stage starts with lr 1e-3 (stage 1: 6000 +
4000 it.; stage 2: 6000 + 14 000 it. at 1e-4). **Meaning:** dispersion does **not** remove κ_Si < 0 on the sample where the
artefact occurs. It removes a large misfit (−37 to −57 %) and moves n_Si to 3.39–3.43, near the literature value, but κ_Si
stays negative over the band under both estimators.

## C4 — interface layer, free τ/a (nb6 `C4_interface_layer`; rMSE L-BFGS from the three starts)

| model (best valid run) | d_gap (µm) | n_gap | κ_Si | n_Si | d₁ / d₂ (µm) | rMSE change | κ_Si > 0 in valid runs |
|---|---|---|---|---|---|---|---|
| two layers (baseline) | — | — | −7.27e-3 | 3.304 | 977.7 / 398.1 | 0 | 0/3 |
| air gap | 8.00 (all starts) | 1 | −1.80e-3 | 3.238 | 958.6 / 409.1 | −5.0 % | 0/7 |
| adhesive, n_g ∈ [1.4, 1.8] | 9.68 | 1.40 (bound) | −2.29e-3 | 3.244 | 958.2 / 407.8 | −4.3 % | 0/6 |
| two layers + τ, a | — | — | +2.72e-3 | **3.644** | 942.1 / 360.9 | −30 % (τ −242 fs, a 1.14) | 2/3 |
| six knots + air gap (C3 + C4) | 6.39 (3 starts) | 1 | −0.0126 (0.2 THz) … +0.005 (3 THz) | 3.316–3.322 | 971.8 / 398.3 | −4.0 % vs six knots | — |

Two air-gap and three adhesive runs (start 2) diverged to unphysical values; they are flagged `valid = false`.
**Meaning:** no interface layer removes the sign either. A gap of 8 µm is preferred, but it lowers rMSE only 5 % and leaves
κ_Si < 0 with n_Si further from 3.4175. A free τ, a gives κ_Si > 0 only at an unphysical n_Si = 3.64 (τ and d₁ + d₂ are
degenerate). Neither "dispersion removes it" nor "an unmodelled interface layer removes it" is supported. The defensible
statement is that κ_Si (and n_Si) of the Sec. 5 bilayer are not determined: their values and signs depend on the model and
the estimator.

## C2 — run 3 without hinge (nb5 `C2_run3_no_hinge`; 10 000 + 15 000 iterations)

Stage 1: κ_Si = **+1.27e-3**, κ_q 4.51e-3, loss 3.6484e-3. Stage 2: Si κ knots −0.00015, −0.00122, +0.00083, +0.00006,
+0.00058, +0.00603 (45 % of the band < 0); loss 2.9187e-3 (with hinge 2.9328e-3, −0.5 %). **Meaning:** on run 3 the
frequency-independent model gives κ_Si > 0 without any penalty, as domain 2.15 suspected. With free knots and no hinge,
κ_Si dips slightly negative at 0.2–0.76 THz. The run-3 result therefore owes its positive κ_Si to the mounting/gap
(stage 1), not to dispersion, and the hinge only cleans up the low-frequency knots.

## C12 — 27-parameter CR bounds and true profiles (nb5 `C12_cramer_rao_and_profiles`)

CR 1σ (σ = 9.04e-4; residual RMS = 6.39σ, multiply by 6.4 for residual-scaled):
quartz n knots 4.2–7.5e-4, quartz κ 3.0e-5–6.0e-4, Si n 1.9–2.2e-3, Si κ 4.4e-5–1.2e-3 (largest at the 0.2/3 THz edges).
**Thicknesses:** d_q 221 nm, **d_air 27 nm**, d_Si 220 nm, d_q + d_Si 39 nm, corr(d_q, d_Si) −0.984. All other parameters
fixed: d_air 15 nm, d_q 6 nm, d_Si 1.9 nm, thickness transfer 96 nm. Mean n_Si 3.3915 ± 1.9e-3 (14 bounds below 3.4175).
**Profiles** (GN warm start + 1500 projected Adam it.; control at Δ = 0: 2.93258e-3):

| mean n_Si | 3.36 | 3.37 | 3.38 | 3.39 | 3.40 | 3.41 | 3.4175 | 3.42 | 3.43 |
|---|---|---|---|---|---|---|---|---|---|
| Δloss vs control | +0.45 % | +0.25 % | +0.09 % | +0.02 % | −0.03 % | −0.04 % | −0.03 % | −0.04 % | 0.00 % |

| d_air offset (µm) | −10 | −6 | −2 | +2 | +6 | +10 |
|---|---|---|---|---|---|---|
| Δloss vs control | +62 % | +27 % | +3.7 % | +5.3 % | +38 % | +79 % |

**Meaning:** the air gap is far better determined than the per-layer split: 27 nm vs 220 nm (marginal), and ±2 µm costs
4–5 % in loss. That supports R15/R22's sentence with a real number. The n_Si profile is flat (within ±0.05 %, the Adam noise)
from 3.38 to 3.43, so the literature value 3.4175 lies on the valley floor. The old "valley floor near 3.37, +4 %" is
wrong and should be deleted; the CR bound of the mean n_Si (1.9e-3) understates this flatness under a non-quadratic,
MAE-type loss.

## C13 — Si knots tied (nb5 `C13_Si_tied`; 15 000 it. from the stage-1 state)

Si n 3.3862, κ 1.4e-6 (at the hinge), d_Si 389.71 µm; loss 2.9534e-3 vs 2.9328e-3 (+0.70 %); MAE/σ 3.25 vs 3.23; 17 vs 27
parameters. ΔAIC / ΔBIC (tied − full): +93 / +31 with L = 3486 white samples; **−15 / −45** with the residual correlation
length (24 samples, n_eff = 146). **Meaning:** once residual correlation is accounted for, the information criteria prefer
non-dispersive Si. Silicon can be presented as non-dispersive at a cost of 0.7 % in loss.

---

## Figures (all re-exported; widths and smallest printed font in `loop1_checks.figure_export`)

| figure | exported width | placed at | smallest printed font | changes |
|---|---|---|---|---|
| fig1_five_layer | 6.58 cm | 6.6 cm (0.5\linewidth) | 7.0 pt | redrawn half width, panels stacked |
| fig2_thin_film | 13.05 cm | 13.3 cm | 7.1 pt | (a) 13 draws at 5–20 µm, 3 at ≥ 30 µm, CR 1σ dashed per parameter, PI-like 25 µm (hollow, range bars); (b) "residual ×10" kept, ±σ (×10) band |
| fig3_exp_fit | 13.20 cm | 13.3 cm | 7.05 pt | (b) σ = 2.9×10⁻⁴ written in the panel, "thinner than the line", residual RMS = 41σ; units on axes |
| fig4_adam_two_phase | 13.30 cm | 13.3 cm | 7.0 pt | fonts only (and two-term axis loadings in (d)) |
| fig5_lbfgs_trajectory | 13.29 cm | 13.3 cm | 7.0 pt | three Phase-I panels (L-BFGS iterates, Adam faint) + loss vs evaluations with linear inset (3.7 % gap); Phase-II panel dropped |
| fig6_qas_fit | 13.26 cm | 13.3 cm | 7.0 pt | fonts, "(arb. u.)"; known cosmetic issue: the (b) legend still touches the −0.1 tick and the stage-1 spike (fix in loop 2, needs a nb5 re-run) |
| fig7_qas_nk | 13.27 cm | 13.3 cm | 7.0 pt | CR 1σ bands + knot error bars (C12), HR-Si 3.4175 line, edge knots hollow, shared legend |

Printed size check: points = font_pt × 13.3 cm / width_cm (6.6 cm for fig1). Every figure is ≥ 7.0 pt.
`fig_sensitivity` and `fig8_qas_scans` were also re-exported at ≥ 7 pt (σ_k → s_k) but are not paper figures (R15).

## Timings, budget, deviations

- Check run times in the delivered notebooks (one torch thread each; several notebooks shared the two cores, so values are up to
  2× inflated): C1 59 s; C11b 322 s; C11c 212 s; C5 rMSE optimum 23 s + bootstrap 236 s; C8 334 s (5.6 ms/it; 4.3 ms alone);
  C7 80 s; C12 Jacobian 0.5 s + profiles 541 s; C2 451 s; C13 153 s; C3 333 s (Adam) + 97 s (L-BFGS); C4 445 s; six knots + gap
  118 s. C6, C9, C10 and C14 take seconds.
- Total CPU ≈ 100 min, over the 60–75 min target. The checks themselves took ≈ 55 min. The rest went on smoke tests (~8 min)
  and on re-executing nb2 (twice) and nb3 (once) to fix figure layout. The re-runs gave identical numbers: everything is seeded,
  and only the timings differ.
- Cuts: C2 10 000 + 15 000 iterations (not 15k + 30k); C3 learning-rate schedule (see C3); C13 15 000 iterations; profiles
  1500 Adam iterations (control drift ±0.05 %).
- **nb1 thread sensitivity:** executed with 1 thread, the five-layer fit ends at rMSE 1.6e-11 instead of 5.1e-8 (round-off
  changes the L-BFGS path). nb1 was re-executed with `NB_THREADS=2` to reproduce the loop-0 numbers the text uses. The 1-thread
  result is in this report only; the floor is the stopping criterion, not the method.
- All other loop-0 numbers (nb2 part B, nb3, nb4, nb5) reproduce bit-for-bit.

## Follow-up (R27, R30): three figure changes, re-executed nb2, nb4, nb5

- **fig5:** the shared Phase-II panel is back, per R30. Layout: (a)–(c) Phase-I planes, (d) shared Phase-II plane with the L-BFGS and Adam Phase-II paths, (e) loss vs evaluations with the linear inset showing the 3.7 % gap. Colours in (e) follow the legend of (d). Smallest printed font 7.0 pt.
- **fig6(b):** shorter legend labels and a lower y-limit; the legend no longer touches the tick.
- **fig2(b):** now shows the C1 variant with the lower rMSE, **(c) everything free + τ, a**: 1.696σ, against 1.704σ for (b).
  - Fitted values: n_PI 1.7836, κ_PI 0.0286 (not at bound), d_PI 25.24 µm, n·d_PI 45.01 µm, n_Si 3.4297, κ_Si 3e-14 (at bound, no CR entry), d_Si 380.24 µm, τ −70.1 fs, a 1.024.
  - CR bounds (τ and a in J): d_PI 38 nm, n·d_PI 37 nm, n_PI 1.9e-3, n_Si 2.2e-3, d_Si 240 nm, τ 0.83 fs.
  - n·d_PI spread: 2.97 % across all 8 starts; 0.22 % across the 2 starts within 0.5 % of the best rMSE.
  - Everything is in `nb2.json → loop1_checks/fig2b_variant`.
- No other numbers changed; the re-runs differ only in timings.
