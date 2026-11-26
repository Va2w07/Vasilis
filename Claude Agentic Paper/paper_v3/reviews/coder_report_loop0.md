# Coder report — loop 0 (2026-11-23)

All five notebooks are in `notebooks/`, executed in place with `jupyter nbconvert --execute` (outputs embedded),
and regenerate every figure (`figs/*.png` 300 dpi + `figs/*.pdf`) and every number (`results/nb*.json`).
Shared style: `code/paperstyle.py` (13.3 cm text width, fonts 7–8 pt, colour-blind-safe palette validated with the
dataviz checker). `notebooks/README.md` maps notebook → figures → JSON. The cell-scripts from which the notebooks
were built, the builder and the smoke tests are in `scratch/` (not deliverables).

## 1. Files produced

| notebook | figures | results |
|---|---|---|
| `nb1_five_layer_simulation.ipynb` | `fig1_five_layer` | `nb1.json` |
| `nb2_thin_film_limits_and_PI.ipynb` | `fig2_thin_film` | `nb2.json`, `nb2_cache.npz` |
| `nb3_quartz_si_adam_landscape.ipynb` | `fig3_exp_fit`, `fig4_adam_two_phase`, `fig_sensitivity` | `nb3.json`, `nb3_trajectories.npz`, `nb3_cache.npz` |
| `nb4_lbfgs_vs_adam.ipynb` | `fig5_lbfgs_trajectory` | `nb4.json` |
| `nb5_quartz_air_si_six_knots.ipynb` | `fig6_qas_fit`, `fig7_qas_nk`, `fig8_qas_scans` (loss scans; optional panel) | `nb5.json`, `nb5_cache.npz` |

New code in `code/`: `paperstyle.py`; `quartzsi_model.py` (the quartz/Si problem shared by nb3 and nb4: data,
float64 forward, loss, three starts, Jacobian); `dispersive_grid64.py` (float64 transcription of
`dispersive_multilayer_grid.py`, same model/loss/fit loop). Existing modules were not modified.

## 2. Protocol deviations from the reference notebooks

1. **Precision.** All fits are float64 (brief). The reference quartz/Si Adam notebook and the quartz/air/Si grid
   notebook are float32. Result: identical to 4–5 significant digits (see §4), so this is a non-issue for the numbers.
2. **nb3 iteration count.** The reference notebook's cells were executed with 50 000 iterations; the brief and the
   old manuscript say 20 000, which is what nb3 uses. The 20 000-iteration values match the old Table 2 exactly.
3. **nb1 / nb2 fits** use the coarse-to-fine L-BFGS of `RobustExtractor` (rMSE loss, low-pass 0.25 / 0.5 / 1 THz
   then full band, hinge against gain) with random multi-starts instead of the Bayesian + Adam pipeline. Starts for
   the sweep are uniform in n ± 0.5, κ ± 0.01 (clipped ≥ 0), d × [0.5, 1.5] around the truth (5 starts); nb1 uses the
   `RobustExtractor` jitter (n ± 0.15, d × e^{±0.35}, κ ∈ [0, 0.01]) around the hand guess (4 starts).
   The PI/Si reference notebook used Adam + `scipy.least_squares` (soft-L1): nb2 uses the same multi-start L-BFGS as
   the sweep instead, as the brief asks.
4. **nb4 L-BFGS mechanics.** `torch.optim.LBFGS` is stepped one iteration at a time (`max_iter=1`) so every accepted
   iterate is recorded; `max_eval` must then be set explicitly (torch otherwise gives the line search zero
   evaluations and the step is always 0 — this cost an hour). When the strong-Wolfe line search can no longer
   decrease the loss, the history is cleared and L-BFGS restarted from the current point (up to 3 restarts); the
   stop reason is in the JSON. A run that diverges (non-finite loss) is repeated once with lr = 0.1 (only rMSE from
   start 2 needed this).
5. **Noise σ definition** (used for every "MAE/σ" and Cramér–Rao number): std of the first 10 ps of the trace in
   question — the *sample* trace for quartz/Si (nb3/nb4: σ = 2.90e-4; the old manuscript said 2.4e-4, window
   unknown), for run 3 (nb5: 9.04e-4) and for PI/Si (nb2 B: 1.35e-3; the reference notebook used t < 8 ps,
   1.47e-3); the *measured reference* for the sweep (nb2 A: 1.58e-3, old manuscript 1.6e-3). The pre-pulse noise
   of these ScanControl traces is bursty, so the value depends on the window by up to ±30 %; the window is in each JSON.
6. **nb5 comparison tolerance.** The brief's expected κ knots have 2 significant digits; I compare against the
   5-decimal printout of the reference notebook and flag a knot only if both the relative deviation exceeds 0.2 %
   and the absolute deviation exceeds 1e-5 (otherwise κ_Si ≈ −1e-5 is "41 % off" by rounding alone).
7. **PI/Si secondary minimum.** The reference notebook's "film collapses to ~16 µm" minimum belongs to its Adam /
   least-squares protocol. With coarse-to-fine L-BFGS none of 8 random starts ends below 20 µm, and a start placed
   at 16 µm converges to the global optimum. What the multi-start exposes instead is (i) a κ-driven local minimum
   (2 of 8 starts: n_PI 1.61, κ_PI 0.0107, d_PI 28.5 µm, rMSE +24 %; and n_PI 1.72, κ 0.0106, d 25.5 µm, +0.9 %) and
   (ii) with a *plain* full-band L-BFGS from 16 µm, a minimum at d_PI 22.3 µm, n_PI 1.93, rMSE +48 %. All in `nb2.json`.

## 3. Timings (two CPU cores, float64, each notebook executed alone)

- nb1: 38 s for 4 starts (2751 evaluations); the hand start needs 361 evaluations / 5.5 s.
- nb2: sweep 27 fits × 5 starts = 159 s; PI/Si 8 starts 12 s. Whole notebook ≈ 4 min.
- nb3: 4.1–4.4 ms per Adam iteration (old manuscript: 8 ms) → 82–89 s per 20 000-iteration start; four 61×61 planes
  ≈ 40 s; whole notebook ≈ 9 min.
- nb4: 4.5–4.6 ms per L-BFGS evaluation (forward + backward); 250–300 evaluations per start; notebook ≈ 3 min.
- nb5: 10 ms per iteration → stage 1 (15 000 it.) 151 s, stage 2 (30 000 it.) 301 s; scans ≈ 1 min; notebook ≈ 9 min.

## 4. Numerical results, and where they differ from the old manuscript

**nb3 (quartz/Si, Adam).** Reproduces the expected optimum: fit 1 n₁ 1.96387, κ₁ 5.41e-3, d₁ 994.28 µm, n₂ 3.4478,
κ₂ −1.44e-3, d₂ 381.41 µm, loss 4.4522e-3 at iteration 20 000 (best point of the three runs, start 3: loss 4.45158e-3,
MAE 4.38e-3 = 15.1 σ). Spread between starts: n₁ 2e-4, d₁ 0.15 µm, n₂ 1.5e-3, d₂ 0.16 µm. Phase boundary at
iteration 1615 (identical to the reference), 92 % of the iterations in Phase II; Phase II lowers the loss by only
22–29 %. Singular values 8.66, 1.54, 0.69, 0.073, 0.019, 0.0016 (σ₁/σ₆ = 5456); ρ(n₁,d₁) 0.973, ρ(n₂,d₂) 0.993,
ρ(κ₁,κ₂) 0.951; κ₁ 277× and κ₂ 885× less visible than n₁ (old text: 250× / 1000×).
Cramér–Rao with σ = 2.90e-4: n₁ 6.1e-5, κ₁ 8.7e-6 (0.16 % of its value), d₁ 31.6 nm, n₂ 2.8e-4, κ₂ 1.9e-5, d₂ 31.1 nm
(old: 5e-5, 7e-6, 30 nm, 2e-4, 1.6e-5, 30 nm — the difference is the larger σ). |κ₂|/bound = 77 (old text: "ninety").
Softplus-constrained fit (β = 5, start 1): n₁ 1.9595, κ₁ 4.78e-3, d₁ 996.6, n₂ 3.4681, κ₂ 1.05e-4, d₂ 379.1 µm, loss
4.4646e-3 = +0.29 % (old softplus row: 1.9595, 0.00478, 996.6, 3.4680, 0.00011, 379.1, 4.4649e-3 — same).
Residual MAE/σ = 15.1 (old text "18σ", because of the σ difference).

**nb4 (L-BFGS vs Adam).** This is the result the writer must handle with care. With the *same* MAE + log-cosh loss,
L-BFGS from all three starts reaches loss 4.617e-3 in 250–300 evaluations and then **stalls** — the strong-Wolfe
line search cannot decrease the loss along the (sub)gradient direction even for steps of 1e-8 (the MAE part is a
V-shaped kink), although the loss decreases *monotonically* along the straight line from the L-BFGS end point to the
Adam end point (both in `nb4.json`). So L-BFGS does **not** reach the 20 000-iteration Adam level (4.4517e-3); its
end point is 3.7 % higher in loss. Adam needs 5000–5900 iterations to reach the L-BFGS end level. The L-BFGS end
point is a different point on the n·d valley: d₁ 999.6–1000.0 µm, d₂ 375.7–376.1 µm, κ₂ = +4.4e-3 (vs Adam 994.3 /
381.4 / −1.5e-3); n₁d₁ and n₂d₂ agree to 0.012 % and 0.034 % over all six end points — consistent with the old
Table 3 row "time only, L-BFGS (320 evaluations): 1000.0 / 375.7", so the old L-BFGS evidently stalled at the same
point; the old claim that it "reached the loss level of 2×10⁴ Adam iterations" is not supported (the old table's
"same loss to within 1 %" vs 3.7 % here). Robustness runs (JSON only): plain MAE stalls the same way (d₁ 997.5–1000);
rMSE converges (smooth) to d₁ 977.7, d₂ 398.1 (start 2 diverged at lr 1, converged with lr 0.1); log-cosh alone
converges to the same 977.7 / 398.1 point from starts 2–3 and to 987.7 / 388.1 from start 1 (400-iteration cap, still
creeping). The valley spread over all fits: d₁ 977.7–1000.0 µm.

**nb5 (quartz/air/Si run 3, six knots).** Reproduces the reference run to ≤ 0.004 % in every n knot and thickness
(quartz d 976.63 µm, air 125.87 µm, Si 388.91 µm; loss 2.9328e-3; κ knots within 1e-5 absolute); nothing flagged.
Stage 1 → stage 2: loss 3.6223e-3 → 2.9328e-3 (−19 %), MAE/σ 3.97 → 3.23. d_air scan: minimum at 126.0 µm, loss at
50 µm 1.62e-2 (5.5× the optimum).

**nb2 A (sweep).** σ = 1.58e-3 (power SNR 44 dB on the reference). Median |error| of the film (n / κ / d, %):
5 µm 1.9 / 25 / 4.8; 7 µm 0.62 / 9.5 / 1.8; 10 µm 0.29 / 8.3 / 0.99; 15 µm 0.11 / 5.7 / 0.37; 20 µm 0.025 / 6.3 / 0.08;
30–100 µm n 0.01–0.05, κ 0.9–4.9, d 0.01–0.05. Converged loss = (0.98–1.00) σ at every thickness, so the errors are
noise-limited. Compared with the old text: same shape (threshold ≈ 10 µm = 0.24 ps round trip; d < 1 % median from
10 µm, n < 1 % from 7 µm); κ at 5 µm is 25 % (old 76 %), and the old "better than 0.1 % from 20 µm" holds for n and d.
**nb2 B (PI/Si).** Best of 8 starts: n_PI 1.712, κ_PI ≈ 0, d_PI 25.84 µm (n·d 44.26 µm), n_Si 3.2592, κ_Si 5.8e-4,
d_Si 400.19 µm (n·d 1304.3 µm), rMSE 2.69 σ, MAE 1.40 σ; 6 of 8 starts agree to 0.1 %. Cramér–Rao: n_PI 2.4e-3,
d_PI 62 nm, n·d_PI 52 nm, ρ(n_PI,d_PI) = −0.93. d_PI is 1.14 µm (+4.6 %, 18 CR-σ) above the profilometer 24.7 µm;
n_PI(THz) is 0.8 % below n∥(633 nm) = 1.726; n·d / d_prof = 1.79. Film round trip 0.295 ps, at the sweep threshold.

**nb1.** All 15 parameters recovered: max error n 1.1e-3 %, d 2.1e-3 %, κ 0.27 %; rMSE 5.1e-8 (float64 floor region);
hand start alone reaches 7.9e-8 in 361 evaluations.

## 5. Open questions / things a reviewer should know

- nb4's L-BFGS stall on the non-smooth loss changes the Sec. 7 narrative: the fair statement is "L-BFGS reaches in
  ~250 evaluations what Adam reaches after ~5500 iterations, then stops at a kink 3.7 % above the Adam end point;
  where it stops on the n·d valley depends on the loss (MAE/log-cosh: d₁ ≈ 1000; rMSE: 978)". A smooth loss or a
  non-smooth-aware line search would fix the stall; not implemented (not in the brief).
- The sample-trace σ (2.9e-4) differs from the old manuscript's 2.4e-4; all σ-relative numbers (MAE/σ, CR bounds)
  move accordingly. The writer should take them from `nb3.json`.
- `fig_sensitivity` is delivered as a separate 13.3 × 4.6 cm figure; it can be dropped or merged if space is short.
- `fig8_qas_scans` (nb5 loss scans) is an optional figure/panel; not requested as a numbered paper figure.
- The reference run 3 notebook used `n_scale = 0.1` (not the module default 0.05); nb5 does the same.
