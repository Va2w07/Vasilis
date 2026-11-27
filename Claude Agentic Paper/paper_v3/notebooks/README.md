# Notebooks that regenerate every figure of the paper

Run from this directory (`notebooks/`); every notebook adds `../code` to `sys.path` and uses the relative paths
`../data`, `../figs`, `../results`.  Execute all of them with

```
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 nb*.ipynb
```

Order: nb3 before nb4 (nb4 reads `results/nb3_trajectories.npz` and `results/nb3.json`); nb6 reads `results/nb3.json`; the others are independent.
Every fit runs in float64.  No Bayesian optimisation anywhere; multi-start around the nominal guess instead.
Figures are written as PNG (300 dpi) **and** PDF at their final single-column size (13.3 cm wide), with the shared
style of `code/paperstyle.py`.  Every number that can appear in the paper is in the `results/*.json` files.

| notebook | content | figures | results |
|---|---|---|---|
| `nb1_five_layer_simulation.ipynb` | synthetic pulse, five-layer stack, 15 parameters, coarse-to-fine L-BFGS from a hand guess (4 starts) | `fig1_five_layer` | `nb1.json` |
| `nb2_thin_film_limits_and_PI.ipynb` | A: film-thickness sweep 5–100 µm with the measured reference pulse and measured noise (3 draws, 5 starts); B: polyimide on Si, multi-start fit, secondary minimum, Cramér–Rao bounds, comparison with profilometer / visible index; C (2026-11-27): PI re-measured 2026-11-25 as a reference–sample–reference bracket | `fig2_thin_film` (panel b from Part C) | `nb2.json` (Part C under `pi_2511`), `nb2_cache.npz` |
| `nb3_quartz_si_adam_landscape.ipynb` | quartz/Si Adam fit (3 starts, 20 000 it.), sensitivity matrix, Cramér–Rao bounds and the κ₂-sign argument, softplus-constrained fit, two-phase PCA landscape | `fig3_exp_fit`, `fig4_adam_two_phase`, `fig_sensitivity` | `nb3.json`, `nb3_trajectories.npz`, `nb3_cache.npz` |
| `nb4_lbfgs_vs_adam.ipynb` | L-BFGS (strong Wolfe, history 20) on the same parameters / loss / starts, every iterate recorded, projected on the planes of Fig. 4; MAE, rMSE and log-cosh losses for robustness | `fig5_lbfgs_trajectory` | `nb4.json` |
| `nb5_quartz_air_si_six_knots.ipynb` | quartz / air / Si run 3, six-knot n(f), κ(f) per layer, two-stage Adam fit (15 000 + 30 000 it.), loss scans | `fig6_qas_fit`, `fig7_qas_nk`, `fig8_qas_scans` (optional panel) | `nb5.json`, `nb5_cache.npz` |
| `nb6_bilayer_model_variants.ipynb` | (loop 1) model variants on the Sec. 5 quartz/Si data: six knots without hinge (C3), interface gap / adhesive layer, free τ and a (C4), six knots + gap | — | `nb6.json` |

Run times on two CPU cores (float64): nb1 ≈ 1 min, nb2 ≈ 4 min, nb3 ≈ 9 min, nb4 ≈ 3 min, nb5 ≈ 9 min.
nb2, nb3 and nb5 have a `RUN_FULL` flag (default `True`); set it to `False` to reload the cached fits from
`results/*_cache.npz` and only redo the analysis and figures.

Code in `code/` used by the notebooks: `Simulate.py` (synthetic pulse), `RobustExtractor.py` (float64 forward model,
scaled parameters, coarse-to-fine L-BFGS), `quartzsi_model.py` (quartz/Si problem shared by nb3 and nb4),
`loss_landscape_viz.py` (PCA planes, phase split, landscape drawing), `dispersive_grid64.py` (float64 six-knot
model, transcription of `dispersive_multilayer_grid.py`), `paperstyle.py` (figure style).

## Loop 1 (2026-11-24): review checks C1–C14 and final-size figures

Every notebook ends with a "Loop-1 checks" part (after the loop-0 analysis) that runs the checks of
`reviews/loop1_domain.md` §5 live and writes them into its JSON under the key `loop1_checks`:

| notebook | checks (key in `results/nbX.json → loop1_checks`) |
|---|---|
| nb1 | figure export only (fig1 redrawn as a half-width 6.6 cm figure, panels stacked) |
| nb2 | C14 noise table of all six traces; C1 PI/Si model variants (n_Si fixed, free τ and a, d_PI fixed); C6 PI/Si CR with τ, a; C11 (a) sweep CR bounds and pulse FWHM, (b) 10 more draws at 5–20 µm with re-polish, (c) PI-like sweep |
| nb3 | C14 (quartz/Si traces); C5 κ₂ against several noise models + block bootstrap; C6 CR with τ, a and the thickness-transfer / κ-split directions; singular-value mini-table; C8 Adam to 6×10⁴ iterations |
| nb4 | C7 Adam learning-rate scan; C9 kink probe and gradient norms; C10 out-of-plane fractions; timing invariants over the 16 end points; new Fig. 5 |
| nb5 | C12 27-parameter CR bounds and true profiles (n_Si, d_air); C2 run 3 without hinge; C13 Si knots tied, AIC/BIC; new Fig. 7 with CR bands |
| nb6 | C3 six knots on the bilayer (λ_κ = 0, Adam protocol and converged rMSE); C4 interface gap / adhesive / τ, a; six knots + gap |

Figures: every figure is exported at its final physical width — 13.3 cm (fig2…fig7, placed at `\linewidth`) and 6.6 cm
(fig1, placed at `0.5\linewidth`) — and `paperstyle.savefig` checks that the smallest font prints at ≥ 7 pt at that width
(it shrinks the canvas if labels stick out); the measured widths and smallest printed font sizes are in each JSON under
`loop1_checks.figure_export`.  The loop-0 notebooks executed with `RUN_FULL = True`; in loop 1 `RUN_FULL = False` (the long
loop-0 fits are reloaded from `results/*_cache.npz`, the new checks run live).  `NB_THREADS` (environment) sets the torch
threads (default 2; loop 1 ran two notebooks in parallel with 1 thread each).  **nb1 must be run with 2 threads** to reproduce
the numbers in `nb1.json`: with 1 thread, round-off changes the L-BFGS path and the noise-free five-layer fit stops at rMSE
1.6e-11 instead of 5.1e-8 (both exact recoveries; the floor is the stopping criterion).  The loop-0 results of nb2–nb5 reproduced bit-for-bit with
1 thread; the loop-1 checks in the delivered notebooks were run with `NB_THREADS=1`.

Loop-1 run times (one thread per notebook, two notebooks in parallel on two cores): see `reviews/coder_report_loop1.md`.
New helper module: `code/loop1_tools.py` (noise statistics, τ/a extractor, CR / GLS / sandwich covariances, block bootstrap, L-BFGS).

## Part C of nb2 (2026-11-27): polyimide re-measured as a reference–sample–reference bracket

Data `data/PI_silicon_ref1_15min_2511.txt` (14:28), `PI_silicon_15min_2511.txt` (14:48), `PI_silicon_ref2_15min_2511.txt` (15:06), 10 000
waveforms each. Part C (appended after the loop-1 checks; Parts A/B and the loop-1 checks unchanged and reproduced bit-for-bit) writes
`results/nb2.json → pi_2511`:

| key | content |
|---|---|
| `data_checks` | pre-pulse σ / white level / peak of all five PI traces; time-zero drift ref2 − ref1 (phase slope 0.2–2 THz, main peak, cross-correlation); amplitude ratio; sample delays vs ref1, ref2, mean and for the old pair |
| `fits` | `<variant>\|<reference>` for variants `all_free` (paper protocol, no τ), `nSi_fixed`, `all_free_tau_amp`, `nSi_fixed_tau_amp` × references `ref1`, `ref2`, `mean`: all 8 starts, best, CR bounds, at-bound flags, n·d spread, starts within 0.5 % |
| `comparison_old_vs_new`, `criteria`, `better`, `verdict` | old (2026-11-12 pair, Part B / C1) vs new, the decision criteria and the verdict |
| `fig2b_variant` | every number for Fig. 2(b) (now: paper protocol, no τ, against the mean reference); `fig2b_alternative_smallest_abs_tau` for the record |

Figure handling: the loop-1 drawing code is now the function `draw_fig2(panel_b, path)`; the loop-1 version (2026-11-12 pair, free τ) is
written to `figs_old/fig2_thin_film_loop1.*`, and Part C writes `figs/fig2_thin_film.*` (panel (a) and layout unchanged; smallest printed
font 7.13 pt at 13.3 cm). If `pi_2511/better` were false, Part C would rewrite the loop-1 panel (b) instead.
Run time of the whole notebook with `RUN_FULL = False`, `NB_THREADS=2`: ≈ 18 min (Part C fits 226 s; the rest is the loop-1 C11 checks,
which are not cached). Report: `reviews/coder_report_pi2511.md`.
