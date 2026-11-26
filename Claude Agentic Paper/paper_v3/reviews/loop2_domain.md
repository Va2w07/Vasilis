# Domain review — loop 2

Reviewer: domain reviewer (THz-TDS parameter extraction, inverse problems, optimisation).

**Material read:**
- `DECISIONS_LOG.md` (brief, R1–R32, D1–D5)
- `reviews/loop1_domain.md`, `loop1_outside.md`, `coder_report_loop1.md`, `writer_report_loop1b.md`
- `main.tex` and `reviews/main_loop2.pdf` (8 pp., read page by page)
- `figs/*.png`
- `results/nb1–nb6.json` (flattened and searched)
- `notebooks/README.md`, the nb5 code cells, `code/dispersive_grid64.py`, `code/Simulate.py`

**Conventions:**
- Line numbers refer to `main.tex`.
- Figure numbers follow the PDF: Fig. 1 = `fig2_thin_film`, Fig. 2 = `fig3_exp_fit`, Fig. 3 = `fig4_adam_two_phase`, Fig. 4 = `fig5_lbfgs_trajectory`, Fig. 5 = `fig6_qas_fit`, Fig. 6 = `fig7_qas_nk`.
- Tags: **[MAJOR]** a claim is wrong or over-stated in a way a referee would reject; **[MINOR]** should be fixed but does not change a conclusion; **[CHECK]** a bounded computation (seconds to minutes).

**Bottom line.** The revision is a large improvement. Every loop-1 MAJOR item is resolved or reduced to a residual minor point, except 2.19, which is now addressed with the wrong mechanism (item 3.1). The numbers are clean: I found no number that is not in a JSON, and only a handful of rounding or labelling slips. The physics is now mostly correct and honestly framed. One substantive problem remains, and it is new evidence against a ruling:
- The paper attributes the 1.2 % difference in quartz n·d between the two mountings to time zero.
- Time zero cannot change an echo spacing. The paper's own counting argument and its own τ-variant (C4) show this.
- The interface-gap model that the coder already ran (C4) reproduces the difference.

So the abstract's "optical thickness to better than 0.1 %" is a within-model precision, not an accuracy.

**Verdict: minor revision**, conditional on item 3.1 (two sentences in Sec. 9 and one clause in the abstract).

---

## 1. Loop-1 MAJOR items — status

| loop-1 item | topic | status | where / remark |
|---|---|---|---|
| 2.1 | PI: impossible n_Si = 3.26 | **resolved** | Sec. 4 l.94, Table 1: stated; timing-invariant argument; τ, a fitted. |
| 2.2 | "18 CR-σ" | **resolved** (small residue) | Deleted. n·d agreement and CR bounds with τ, a in J are given. Residue: "best starts" means 2 of 8, and the profilometer uncertainty is still not given (item 3.12). |
| 2.3 | sweep vs PI not comparable | **resolved** | PI-like sweep (C11c) in Fig. 1(a) and l.94. Small residue on the top axis (item 4.3). |
| 2.4 | κ split PI film/substrate | **resolved** | "split undetermined"; κ_Si "at bound", no CR entry; "local minimum" row gone; footnote d on 42.63. |
| 2.5 | quartz/Si n_Si also biased | **resolved** | l.113 (+0.9 %); step (iv) l.155 (≈100 bounds); Sec. 9 estimator spread. |
| 2.6 | κ "visibility" scaling artefact | **resolved** | l.144 holomorphy sentence; s_k of J/√L; relative-precision wording. |
| 2.7 | valley structure / counting argument | **resolved** | l.141, as R17. A caveat on its scope is needed (item 3.1). |
| 2.8 | non-white noise, CR argument | **resolved** | l.85 (σ is a bookkeeping scale, white part 6–72× smaller); step (iv) rests on 41σ and the C5 bounds. One range slip (item 2.2). |
| 2.9 | step (vi) a hypothesis | **resolved** | Rewritten with C3, C4 and the τ variant. One wording issue (item 3.3). |
| 2.10 | Adam not converged | **resolved** | l.113 "budget iterate on a plateau", "follow the same valley". |
| 2.12 | "twenty times" depends on η | **resolved** | l.161, range over η. One omission (item 3.5). |
| 2.15 | six knots credited with the sign | **resolved** | Withdrawn; C2 in l.202. |
| 2.16 | band-edge knots, ripple | **resolved** (interpretation over-stated) | Fig. 6 hollow edge knots, HR-Si line, ripple sentence; Si non-dispersive (C13). The ripple interpretation is not significant at model-error scale and sits at the knot-grid resolution (item 3.7, C15). |
| 2.17 | "valley floor near 3.37" grid artefact | **resolved** | Replaced by the C12 profile, l.202. |
| 2.18 | air-gap conditioning | **resolved** | Table 3 (27 vs 221/220 nm, ρ = −0.98, 39 nm); "eight times" l.202. "Order of magnitude" in the abstract (item 2.7). |
| 2.19 | sample identity / accuracy of n·d | **partly — wrong mechanism** | l.206 attributes the 1.2 % to time zero. Time zero cannot change n_j d_j; the C4 interface gap reproduces the difference (item 3.1; disagreement with R23/R29). |
| 2.21 | "noise-limited" non sequitur | **resolved** | C11b; l.91; CR curves in Fig. 1(a). |
| 2.23 | protocol of Sec. 4 | **partly** | Scaling, coarse-to-fine and hinge are in l.78. Missing: starts and start box, hinge weight, quartz/Si starts and iteration budget, which were dropped from Table 2 in this loop (item 3.11). |
| 3.1 | fonts below 7 pt | **resolved** | All figures at \linewidth, smallest printed font 7.0 pt (`figure_export` in the JSONs); legible in the PDF. |
| 3.2 | no uncertainties in Sec. 8 | **resolved** (interpretation note) | CR bands and error bars in Fig. 6, CR column in Table 3, residual factor 6.4 in the caption. The bands of n_Si are common-mode (item 3.8). |

Loop-1 MINOR items (spot-checked):
- **Fixed:** 1.1, 1.3, 1.4, 1.11, 1.12, 1.14, 2.11, 2.13, 2.20, 2.22, 2.25, 2.26, 3.3–3.9.
- **Dropped with R32's consent:** 1.8.
- **1.9:** acceptable ("<10⁻⁴").

---

## 2. Number audit

**Method.**
- **Script** (`/tmp/a2/audit.py`): it flattens nb1–nb6.json (7417 leaves) and extracts every number from the abstract-to-backmatter body of `main.tex`. Each number is matched at its printed precision, with unit scalings. The automatic pass is permissive, as in loop 1.
- **Hand check:** every non-trivial number in the prose, every caption and all three tables, against the named JSON field. This includes:
  - the PI τ-variant numbers (`nb2 loop1_checks/fig2b_variant`, `C1_PI_model_variants/b_…`);
  - the C5 ratios (`nb3 C5_kappa2_vs_noise_models`);
  - the C7 η ranges (`nb4 C7_adam_learning_rate`);
  - the C9 kink numbers;
  - the C10 fractions;
  - the C11 sweep numbers;
  - the C12 bounds and profile;
  - C13;
  - the C3/C4 rows (`nb6`);
  - the 16-end-point invariants (`nb4 end_points_timing_invariants`).

**Verified correct**, in brief:
- **Fig. 1 caption:** 1.58e-3; 13 and 3 draws; 2.9 / 0.88 / 0.12 % (13-draw medians); PI-like 0.043 / 0.16 / 17 %; 1.35e-3.
- **Sec. 4:** five-layer 0.002 / 0.3 %; std/CR 0.77–1.26 ("within 25 %"); 1 % crossing between 7 and 10 µm; FWHM 0.355 ps; round trip at 10 µm 0.240 ps (= ⅔ × 0.36).
- **PI numbers:** n_Si 3.26; τ-fit misfit −37 %; n_Si 3.4297 (+0.36 %); d_PI 25.24 (+2.2 %); n_Si-fixed variant 25.25 / −65.4 / 1.704σ; n·d spread 0.217 % (2 starts) and 2.97 % (8 starts).
- **Table 1:** all CR entries (fig2b_variant `cramer_rao`).
- **Sec. 5:** <0.02 % over 6×10⁴ iterations (max deviation +0.019 %); +0.9 %; 41σ (41.3).
- **Fig. 2 caption:** σ 2.90e-4; peak/σ 7.9×10³ (7911); MAE 15σ (15.1).
- **Table 2:** all rows match:
  - Adam = `nb3 best`, which is the lowest-loss iterate of start 3;
  - softplus and L-BFGS start 1;
  - rMSE, gap and adhesive from `nb6 C4 summary`;
  - CR row = `nb3 cramer_rao/bound_physical`.
- **Table 2 caption:**
  - 0.28 % (4.46456/4.45223);
  - −5.0 / −4.3 %;
  - six-knot rMSE 989.3 / 387.3, n_Si 3.392–3.402, κ_Si −0.0060 … −0.0030, −57 %;
  - six knots + gap 6.4 µm, −0.0126 / +0.0050, −58.6 %.
- **Sec. 6.1:** invariants 0.026 / 0.031 / 0.057 % (< 0.06 %); d_q spread 2.25 %; κ_Si sign change.
- **Fig. 3 caption:** Phase I 1–1615; −81.1 / 87.7 / 91.5 %; Phase II 99.7 % (92.4 + 7.3); −22–29 %; all six s_k and their loadings (`singular_value_table`).
- **Sec. 6.2:** 77 (76.8); ≈100 (107.8); below 2 at Adam (1.86 / 0.71 / 0.66); κ_Si −0.007 … +0.004; τ-variant n_Si 3.64; softplus 0.3 %.
- **Sec. 7:**
  - L-BFGS: 123–154 evaluations to the stall, 250–299 in total;
  - Adam: 4998–5929 iterations at η = 10⁻⁴ and 357–460 at η = 3×10⁻³, none diverging;
  - 3.7 % gap (3.71 %);
  - kink: slope +0.023 against −0.036, two sign changes;
  - rMSE: gradient norm 2.6–6.0 → 1.0–2.7e-6; d₁ spread 0.21 nm.
- **Fig. 4 caption:** 12–55 %; 97.2–99.8 %; 0.050 / 1.13.
- **Fig. 5 caption:** σ 9.04e-4; −19.0 %; 3.62e-3 → 2.93e-3; MAE 3.97σ → 3.23σ.
- **Fig. 6 caption:** factor 6.39.
- **Table 3:** all 24 knots, the stage-1 column, n·d at 1 THz (1926.0 / 1319.3), ρ = −0.984, 39.1 nm, and all CR ranges. The ranges are quartz n 4.2–7.5e-4, κ 0.30–6.0e-4, Si n 1.85–2.22e-3, Si κ 0.44–11.6e-4, and d 221 / 27 / 220 nm.
- **Sec. 8:** 8.2× ("eight times"); c/2d_air = 1.191 THz; +0.70 % (C13); C2 stage-1 κ_Si +1.27e-3; profile within ±0.05 % over 3.38–3.43.
- **Sec. 9:** 2.25 % / 5.8 %.

Slips found:

**2.1 [MINOR] Table 2 caption (l.122): "others within 2×10⁻⁴ in n_q, 0.16 µm in d".**
- **What the row is.** The Adam row is the *lowest-loss iterate* of start 3 (`nb3/best`, d_q 994.215). It is not start 3's end point (994.248).
- **Distance of the others from it:** start 2's end point lies 0.19 µm away in d_q (994.401) and 1.9×10⁻⁴ in n_q.
- **What 0.16 µm is:** the max–min spread of the three end points; their n_q spread is 2.3×10⁻⁴.
- **Fix:** "Adam: lowest-loss iterate (start 3); the three end points spread by 2.3×10⁻⁴ in n_q, 0.16 µm in d, 1.5×10⁻³ in n_Si".

**2.2 [MINOR] Sec. 6.2(iv) (l.155): "4–10 at the least-squares point".** The C5 values at the rMSE optimum are:
- residual-scaled: 9.6;
- HAC sandwich: 3.8 (1 ps) and 3.6 (3 ps);
- block bootstrap: 4.4.

Write "4–10" as "3.6–9.6", or "about 4 to 10".

**2.3 [MINOR] Fig. 4(e) y-axis label reads "loss, Eq. (4)".** After the merge of Eqs. (1)–(2) the loss is Eq. (3). The label is hard-coded in nb4; re-export the figure. The nb4 C7 description string also says "Eq. (4)", which is cosmetic.

**2.4 [MINOR] Fig. 4 caption (l.167): "the iterates are offset in κ_q from this plane".** This is not what C10 shows.
- The κ_q difference between L-BFGS and Adam (0.00305 against 0.00541) lies *in* the Phase-II plane.
- I recomputed the L-BFGS end point's in-plane coordinates from the κ's alone, with the plane loadings u = 0.97κ_Si + 0.24κ_q and v = −0.97κ_q + 0.24κ_Si about Adam's optimum: (0.966, 0.581), |·| = 1.127. This reproduces C10's 1.13.
- The 0.05 out-of-plane offset comes from the n, d components: thickness transfer and n_Si contribute ≈0.02 in nominal-scaled units. The rest is rounding of the κ loadings.

Fix: "the end points lie 0.05 off this plane (n, d along the thickness-transfer valley), 1.13 in it, so the contours are Adam's".

**2.5 [MINOR] Sec. 7 (l.161): "crosses the Phase-I bowl in about thirty iterations" is not in any JSON (R10).** It is inherited from R11. Either the coder writes the count to nb4.json (CHECK C17), or the phrase is replaced by the JSON number that exists ("reaches its final level in 123–154 evaluations").

**2.6 [MINOR] Fig. 1 caption: "SNR 44 dB".**
- **What it is:** `power_snr_dB` = 10 log₁₀(⟨x²⟩/σ²), averaged over the whole 80-ps window. That is not a standard THz-TDS SNR and it depends on the window length.
- **What the paper uses elsewhere:** Fig. 2 quotes peak/σ.
- **Fix:** "peak/σ = 2.1×10³" (`partA/peak_over_sigma` = 2107). That keeps the two captions consistent.

**2.7 [MINOR] Abstract (l.47) "an order of magnitude better" against Sec. 8 (l.202) "eight times better" (221/27 = 8.2).** Use "eight times" in both.

**2.8 [MINOR] Sec. 8 (l.202): "n_Si is 0.8 % low".**
- Correct for the knot mean (3.3915, −0.76 %).
- The value at 1 THz is 3.3922 (−0.74 %) and the tied-Si fit gives 3.3862 (−0.92 %).
- Fix: "the mean Si knot is 0.8 % low".

No other discrepancies. The literature constants (3.4175, 1.726, 24.7) and the instrument metadata (1000 waveforms, 15 min, PI-2525) are acceptable, as in loop 1.

---

## 3. Physics and presentation — fresh-referee reading

**3.1 [MAJOR] Sec. 9 (l.206) and abstract (l.47): the ~1 % accuracy of n·d is attributed to time zero, but time zero cannot change an optical thickness. The interface model can, and the paper's own data show it.**

*The argument.*
1. **Time zero cannot move n·d.** By the counting argument of Sec. 6.1, each n_j d_j is fixed by an *echo spacing* inside the sample trace, and that spacing is independent of the reference time zero. A time-zero error moves only Σ(n−1)d, and hence only d_q + d_Si. The paper's own runs confirm this:

   | run | n_q d_q (µm) | n_Si d_Si (µm) | d_q + d_Si (µm) |
   |---|---|---|---|
   | quartz/Si, τ free (`nb6 C4 tau_amp`, τ = −242 fs) | 1952.5 → 1952.2 | 1315.0 → 1315.0 | 1375.6 → 1303.0 |
   | PI, τ free | — | 1304.3 → 1304.1 | — |

2. **The interface model does move n·d.** In the C4 runs of the same Sec. 5 data (`nb6 C4 summary`), adding an interface layer shifts both optical thicknesses:

   | model | n_q d_q (µm) | n_Si d_Si (µm) | d_q + d_Si (µm) |
   |---|---|---|---|
   | 8.0-µm air gap | 2.0148 × 958.59 = 1931.4 | 3.2382 × 409.10 = 1324.8 | 1367.7 (+ 8.0 gap = 1375.7, the two-layer total) |
   | 9.7-µm adhesive | 1929.1 | 1322.9 | 1366.0 |

3. **The run-3 group optical thickness of the quartz** (my CHECK C15) is n_g d = 1928.2 µm at 1 THz (n = 1.9721, n_g = 1.9743). Phase n d is 1926.0.

*What follows.*
- The interface gap alone brings the Sec. 5 quartz n·d from 1.2 % above the Sec. 8 value to within 0.05–0.17 %.
- It brings d_q + d_Si to within 2.2 µm (0.16 %). Only this last residue is time-zero-like (≈7 fs).
- For Si, the gap over-corrects: +0.3–0.4 % against −0.33 % without it.
- The physical reason is plain. A sub-wavelength interface layer changes the phase of the quartz back-face reflection, and that phase depends on frequency. Its sign is opposite for quartz→Si and quartz→air. A phase change of this kind shifts the apparent echo delay. The counting argument assumes real, instantaneous Fresnel coefficients.

*Consequences.*
- The precision of 0.01–0.06 % across estimators (16 end points) is real.
- The accuracy of about 1 % is set by the interface model, not by time zero.
- The abstract's "pulse timings fix each layer's optical thickness and the total thickness to better than 0.1 %" reads as an accuracy claim, and so over-claims.
- The total thickness also needs "given time zero" (C6: freeing τ makes d_q + d_Si 23× less certain). The conclusions say this; the abstract does not.

*Fix.* Approximately 50 words, no new figure:
- **Abstract:** "…fix each layer's optical thickness and, given time zero, the total thickness reproducibly to 0.1 %; an unmodelled interface gap shifts them by about 1 %."
- **Sec. 9, replacing "likely set by time zero":** "The same plates measured twice differ in quartz n d by 1.2 %. An 8-µm interface gap in the bilayer model removes most of this (to 0.2 %), so the accuracy, about 1 %, is set by the interface model. Time zero (here ≈ 7 fs) limits only the total thickness."
- **Sec. 6.1:** add "for interfaces without intermediate layers" to the counting argument.
- **Traceability:** the coder writes these products to nb6.json and nb5.json (CHECK C16).

This conflicts with R23/R29; see §6.

**3.2 [MINOR] The κ_Si wording over-states in two places (abstract l.47, Sec. 6.2 l.155–157).**
- **(a) "A softplus constraint … resolves it."**
  - A constraint *imposes* the physical sign at a cost of 0.3 % in loss. It does not resolve the model error that drives it.
  - The constrained fit also moves n_Si *further* from the literature value: 3.4681 (+1.5 %) against 3.4476 (+0.9 %) unconstrained (Table 2). The text does not say so.
  - Write "a softplus constraint imposes the physical sign at 0.3 % in loss (n, d move along the thickness-transfer valley, n_Si to +1.5 %); richer models do not". Use "imposes" in the abstract.
- **(b) "it lies along a direction the data do not determine."**
  - At the least-squares point, κ_Si < 0 survives every model-error-aware uncertainty at 3.6–4.4 bounds, and bootstrap gives κ₂ < 0 in 99 % of resamples (C5).
  - The accurate statement: it lies along the most weakly determined direction (s₆), and its *estimate* depends on how the estimator weights the model error.
  - Replace "the data do not determine" with "the data determine only weakly".

**3.3 [MINOR] Sec. 6.2(vi) (l.155): "Richer models keep the sign" is not true of the richest model run.**
- Six knots + gap (Table 2 caption, `nb6 C3_C4`) gives κ_Si from −0.0126 at 0.2 THz to +0.0050 at 3 THz, with an interior mean of +3×10⁻⁴.
- The sign therefore *changes across the band*. That is equally unphysical, but it is not "kept".
- Write "Richer models do not make κ_Si physical: six knots leave it negative across the band, adding a gap makes it swing from −0.013 to +0.005, and a free τ gives κ_Si > 0 only with n_Si = 3.64".

**3.4 [MINOR] Free τ is accepted for PI (Sec. 4) and rejected for quartz/Si (Sec. 6.2), with no stated criterion.**
- A referee will ask why the same nuisance parameter is physical in one case and not in the other.
- Its behaviour differs between the two samples:

  | sample | τ | a | effect on n_Si |
  |---|---|---|---|
  | PI | −65 to −70 fs | 1.02 | restores 3.43 |
  | quartz/Si | −242 fs | 1.14 | moves to 3.64 |

- On quartz/Si, τ trades against d_q + d_Si (C6), and the misfit there is not a timing error.
- Add half a sentence to (vi): "…unlike for the PI sample, where τ is tens of fs and restores the substrate index".
- In Sec. 4, the τ nuisance parameter is defined (formula, sign convention consistent with ñ = n − iκ) and motivated (scans 23 min apart). Good.

**3.5 [MINOR] Sec. 7 (l.161): the Adam/L-BFGS comparison omits the result most relevant to "which to use".**
- C7: at η = 3×10⁻³, Adam reaches loss 4.5×10⁻³ in 621–662 iterations. That level is *below* the L-BFGS stall, and L-BFGS never reaches it on Eq. (3).
- Also say what "Adam needs 5000–5900 iterations" is needed *for* ("to reach the L-BFGS final level").
- Suggested: "…360–460 at 3×10⁻³, none diverging, and continues below the L-BFGS stall within ~650."
- The framing "the advantage of L-BFGS is not a fixed speed-up but that no step size needs tuning" then stands. Note that Adam did not diverge across the 30× range of η tested, so the tuning burden is modest here. That is fair to state.

**3.6 [MINOR] Sec. 7 (l.163): "the known failure of BFGS on non-smooth functions [Lewis2013]" slightly mis-cites Lewis & Overton.**
- Their message is that BFGS with an *inexact weak-Wolfe* line search is often effective on non-smooth functions.
- Breakdowns occur when the line search cannot satisfy its conditions at a kink. They argue that the strong-Wolfe curvature condition is inappropriate for non-smooth functions.
- PyTorch's L-BFGS uses strong Wolfe (nb4 `line_search = strong_wolfe`).
- Write "…the line-search breakdown expected of a strong-Wolfe search at a non-smooth point [Lewis2013]".
- This also tells the reader what to try (weak Wolfe, or a smooth loss).

**3.7 [MINOR] Sec. 8 (l.202): "its knots ripple with the gap's free spectral range, c/2d_air ≈ 1.19 THz (etalon mismodelling)" is stated as fact. The data cannot support it at that strength.**
- **(a) Resolution.** The knot spacing is 0.56 THz, so the shortest period the spline can represent is 1.12 THz. A "ripple of period 1.19 THz" cannot be told apart from knot-to-knot alternation, the generic behaviour of an under-determined spline.
- **(b) Significance (CHECK C15, done).** Knot-difference contrasts from the 27-parameter covariance:

  | contrast | value | pre-pulse bounds | residual-scaled bounds (×6.4) |
  |---|---|---|---|
  | interior n_Si alternation | (5.4 ± 0.9)×10⁻⁴ | 5.8 | 0.9 |
  | κ_Si alternating contrast | — | 2.0 | 0.3 |

  At the model-error scale that the paper itself adopts, neither is significant.
- **Fix:** "the Si knots show a weak ripple near the resolution of the knot grid, possibly the gap etalon (c/2d_air = 1.19 THz)". R22 asked for "likely"; even that is generous (§6).

**3.8 [MINOR] Fig. 6(b), Table 3 and Sec. 8: the n_Si CR bands look as if the literature value is excluded, while the text says it "lies on the valley floor". The caption must reconcile the two.**
- C15: the six n_Si knots are correlated at ≥ 0.98 (0.85 for the 3-THz knot). The ±1.9×10⁻³ bands are almost entirely a *common offset*, the thickness-transfer direction.
- Knot-to-knot differences are determined about ten times better (±2×10⁻⁴).
- The knot mean is 13.5 pre-pulse bounds from 3.4175 (C12), but 2.1 residual-scaled bounds. That is consistent with the flat profile.
- Add to the Fig. 6 caption: "Si n knots are ≥ 0.98 correlated: the band is a common offset along the thickness-transfer valley; with the residual RMS the literature value lies 2 bounds from the knot mean."

**3.9 [MINOR] The hinge does not "keep κ ≥ 0" (Sec. 8, l.173), and the κ_Si values are not those of HR-Si.**
- The hinge penalises negative κ; it does not forbid it.
- The Si κ spline dips to about −2×10⁻⁴ near 2.1 THz (Fig. 6(d)). The 0.76-THz knot is −1.4×10⁻⁵, and 3.9 % of the band is negative (`k_negative_fraction_in_band`).
- The in-band κ_Si knots (up to 1.2×10⁻³; α ≈ 0.1–0.5 cm⁻¹) are 3–10× the HR-Si value [Dai2004]. They are consistent with zero only at the residual-scaled bound.
- Fix: "penalises κ < 0" in l.173, and one clause in Sec. 8: "κ_Si is consistent with zero only at the model-error scale".

**3.10 [MINOR] Eq. (1)/(2) (l.64–73): the transfer function as written is missing the air-path factor, so the counting argument does not follow from the equations.**
- `Simulate.rts_batched` returns t·exp(i n₀ D_j), i.e. the slab transmission referenced to the same thickness of ambient.
- Without that factor, the main-pulse delay would measure Σ n_j d_j, not Σ(n_j − 1)d_j as stated in Sec. 6.1.
- Fix: write t_j e^{i n₀ D_j} in Eq. (1), or add "(referenced to the same thickness of ambient)" after "stack transmission T(ω)".

**3.11 [MINOR] The Sec. 3 protocol is not sufficient to reproduce Secs. 4–6 (loop-1 item 2.23, partly).** Missing, at the cost of one sentence:
- **Sweep:** 5 starts in a box centred on the truth (n ± 0.5, κ ± 0.01, d × [0.5, 1.5]). Centring on the truth is worth stating, because it flatters the sweep. Noise is added to the sample only.
- **PI:** 8 starts; all six layer parameters plus τ and a are free.
- **Hinge weight** for the sweep and PI.
- **Quartz/Si:** the three starts (n_q 1.982 / 1.900 / 2.010, d_q 1000 / 900 / 900 µm; Si 3.40 / 0.0054 / 400 µm; κ_q 0.0071) and the 2×10⁴-iteration budget. These were in the loop-1 Table 2 caption and were cut. Fig. 3 shows the start dots, so the reader needs them.
- **Data availability (l.215):** says "on request". Brief item 4 delivers notebooks that regenerate every figure, so cite them: "notebooks regenerating every figure: [repository]". This is an author `\todo` already, but the sentence should say what exists.

**3.12 [MINOR] PI paragraph (l.94).**
- "The best starts agree on n_PI d_PI to 0.2 %" means *two* of eight starts (Table 1 caption). Say "the two best starts" in the text.
- The n_Si-fixed variant, also quoted, has 6/8 starts within 0.5 % (`b_… n_starts_within_0p5pct`). That is the more robust support and could be cited instead.
- The profilometer uncertainty and the stylus-to-THz-spot offset are still not stated (loop-1 item 2.2). A 2 % (0.5 µm) difference may lie within them. One clause, or an author `\todo`.
- "A PI-like simulation puts the 1 % limit in d at 10–15 µm, so the PI errors are systematic": correct.

**3.13 [MINOR] Fig. 3 caption (l.150): s₃ is mislabelled.**
- "s₂ = 1.54, s₃ = 0.69, optical thickness moved between layers" fits s₂ (+0.59 d_Si +0.51 n_Si −0.50 n_q −0.37 d_q).
- It does not fit s₃ (+0.81 d_q −0.45 n_q +0.31 d_Si). s₃ is mainly the quartz physical thickness against its index, i.e. the (n−1)d combination fixed by the main-pulse delay.
- Fix: "s₂ = 1.54, optical thickness moved between layers; s₃ = 0.69, physical thickness at fixed delay (mainly d_q)".

**3.14 [MINOR] Sec. 4 (l.91) and Fig. 1(a): the round-trip criterion is specific to the sweep film.**
- "Exceeds 1 % … where the round trip falls below two thirds of the reference FWHM" holds for the sweep's high-contrast film (n_f = 3.6, Fresnel product 0.135).
- The PI-like film crosses 1 % at a round trip of 0.12–0.17 ps, i.e. one-third to one-half of the FWHM (C11c). The threshold depends on contrast (D1).
- The top axis of Fig. 1(a) is computed with n_f = 3.6. The PI-like marker at 25 µm therefore sits under "0.6 ps", although its round trip is 0.29 ps.
- Fix: "for this film" in l.91, and in the Fig. 1 caption: "top axis for the simulated film; PI-like round trip at 25 µm, 0.29 ps".

**3.15 [MINOR] Sec. 8 (l.202): "the gap, timed by two interfaces".** The reason the gap is well determined is that its index is fixed (n = 1). Its echo delay gives d_air directly, whereas each plate's n–d split rests on the Fresnel amplitudes. Write "the gap, whose index is known, is timed directly and is eight times better determined…".

**What is now correct and well argued (no action):**
- The counting argument (apart from the scope caveat in 3.1).
- The holomorphy argument for κ against n.
- The non-white-noise caveat and step (iv).
- "Budget iterate on a plateau".
- The kink evidence (one-sided slope, monotonic line).
- The withdrawal of "dispersion makes Si lossless".
- Si presented as non-dispersive (C13).
- The C12 profile replacing the grid scan.
- The PI reframing around τ.
- The softplus definition (formula, β, why softplus rather than a rectifier).
- τ is defined.
- The log-cosh term is quantified.
- The combined-loss sentence is qualified.

**Conclusions (l.206).** They are very short (R32) but not wrong, except for "likely set by time zero" (3.1). The brief's main goals are stated in one line each: limits, capabilities, loss profile, and Adam against L-BFGS. The thin-film limit is not restated. Consider "a thin film is recovered down to a round trip of about half the pulse duration" if a line is free.

---

## 4. Figures against captions

All figures are set at \linewidth. The smallest printed font is ≥ 7.0 pt (coder `figure_export`) and everything is legible in `main_loop2.pdf`. All panels are described in the captions. Items:

**4.1 [MINOR] Fig. 4(e):** the y-label "Eq. (4)" should read Eq. (3) (item 2.3).

**4.2 [MINOR] Fig. 4 caption:** "offset in κ_q" (item 2.4).

**4.3 [MINOR] Fig. 1(a):** the top axis applies to the sweep film only, so the PI-like marker is mis-timed against it (item 3.14). The legend says "CR 1σ" while outside item 16 reserves σ for the noise floor; "CR bound" in the legend would be consistent. The same applies to Fig. 6. This is trivial.

**4.4 [MINOR] Fig. 6:**
- The bands need the common-mode remark (item 3.8).
- Panel (d) shows the spline going negative near 2.1 THz while the text says the hinge keeps κ ≥ 0 (item 3.9).
- The caption should name (a)–(d) individually as quartz n, Si n, quartz κ, Si κ. Currently "(a), (b) n(f), (c), (d) κ(f)" leaves the column assignment to the panel titles. Acceptable.

Nothing else in a figure contradicts the text. Fig. 2(a) shows the start-1 fit while Table 2 lists the start-3 best iterate; they are indistinguishable at plot scale (no action). Fig. 3 and Fig. 4(a)–(c) share planes and axis loadings as stated. The 3.7 % inset in Fig. 4(e) is readable.

---

## 5. Checks

**[CHECK] C15 (done by the reviewer, 3 s).**
- **Script:** `/tmp/chk/c15.py`. It reads the data and `nb5_cache.npz` only; nothing in the project was written.
- **What it does:** recomputes the 27-parameter Jacobian at the stage-2 optimum. The residual RMS is 6.388σ and reproduces C12.
- **Results:**
  - corr(n_Si knots) ≥ 0.98 (0.85 with the 3-THz knot); κ_Si knots |corr| ≤ 0.46.
  - Interior n_Si alternation (5.4 ± 0.94)×10⁻⁴; κ_Si alternating contrast 2.0 bounds.
  - Quartz at 1 THz: n_g = 1.9743, n_g d = 1928.2 µm.
- **Request:** the coder adds this to `nb5.json → loop1_checks` (or a `loop2_checks` key) if items 3.1, 3.7 or 3.8 quote it (R10).

**[CHECK] C16 (seconds, coder).** Write to `nb6.json` the products n_q d_q, n_Si d_Si, d_q + d_Si and Σ(n−1)d for the C4 variants (baseline, gap, adhesive, τ/a). From my recomputation:

| variant | n_q d_q | n_Si d_Si | d_q + d_Si | Σ(n−1)d |
|---|---|---|---|---|
| baseline | 1952.3 | 1315.3 | 1375.8 | 1891.8 |
| gap | 1931.4 | 1324.8 | 1367.7 | 1888.4 |
| adhesive | 1929.1 | 1322.9 | 1366.0 | 1886.0 |
| τ, a | 1952.2 | 1315.0 | 1303.0 | 1964.3 |

Also write to `nb5.json` the run-3 quartz group optical thickness n_g d at 1 THz. Needed for item 3.1.

**[CHECK] C17 (seconds, coder).** In nb4, the number of L-BFGS iterations (not evaluations) to complete 95 % of the loss decrease, per start (the Phase-I criterion of Sec. 6.1). Write it to nb4.json to back "about thirty iterations" (item 2.5), or drop the phrase.

---

## 6. Disagreements with rulings

- **R23 / R29 (Sec. 9 attributes the 1.2 % quartz n·d difference to time-zero reproducibility).** I disagree, with evidence that was not considered when the ruling was made (item 3.1):
  - A reference time-zero offset cannot change an echo spacing, so it cannot change any n_j d_j. That is the paper's own counting argument, and the τ-variant of C4 confirms it (n_q d_q 1952.5 → 1952.2 while d_q + d_Si moves by 73 µm).
  - The interface-gap model of C4, already run on the same data, moves n_q d_q to 1931.4 µm (adhesive: 1929.1). That is within 0.05–0.17 % of the run-3 group value of 1928.2 µm.
  - The user's answer in R29 that the plates were "clamped in contact, no adhesive" is fully consistent with an 8-µm air film. It is the expected outcome of clamping two non-optically-flat plates.
  - **Request:** amend R23/R29 so that Sec. 9 attributes the n·d accuracy to the interface model (time zero limits only d_q + d_Si, residue ≈ 7 fs). Qualify the abstract's 0.1 % as within-model reproducibility.
  - The physics change is small, but as written a THz referee would catch it immediately.
- **R22 (ripple "a likely gap-etalon mismodelling").** I proposed that wording in loop 1; C15 now shows it is not significant at the model-error scale and sits at the knot-grid resolution. Request: "possibly" instead of "likely" or the bare assertion (item 3.7). Minor.
- **R32 ("held constant outside the band" cut).** Accepted. Note for the authors only: the knots are clipped at the band edges (`interp_basis` clips f). Out-of-band spectral content is therefore fitted by the edge knots, which is part of why they are "unconstrained". If a line is free, it explains Fig. 6's flat extensions.
- **No other disagreement.** R16–R21, R24–R28, R30 and R31 are implemented correctly and I agree with them.

---

## Counts

[MAJOR] 1 · [MINOR] 23 · [CHECK] 3 (C15 done, C16, C17) · total 27.

- **MAJOR:** 3.1.
- **MINOR:** 2.1–2.8, 3.2–3.15. Items 4.1–4.4 are cross-references to 2.3, 2.4, 3.8, 3.9 and 3.14 and are not counted twice.
- **CHECK:** C15, C16, C17.
