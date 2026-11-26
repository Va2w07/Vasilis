# Domain review — loop 1

Reviewer: domain reviewer (THz-TDS parameter extraction, inverse problems, optimisation).
Material: `main.tex`, `reviews/main_loop1.pdf` (8 pp.), `figs/*.png`, `results/nb1–nb5.json`, notebook code cells (nb2–nb5, extracted and read), `code/quartzsi_model.py`, `code/Simulate.py`, the raw traces in `data/`, `coder_report_loop0.md`, `writer_report_loop1.md`, `loop1_outside.md`, `DECISIONS_LOG.md` (brief + R1–R15).
Tags: **[MAJOR]** a claim is unsupported, wrong or mis-stated; **[MINOR]** should be fixed, does not change a conclusion; **[CHECK]** a bounded computation the coder must run before the claim can stand (compute estimates are for 2 CPU cores and use the timings in the coder report).
Figure numbers follow the PDF: Fig. 1 = `fig2_thin_film`, Fig. 2 = `fig3_exp_fit`, Fig. 3 = `fig4_adam_two_phase`, Fig. 4 = `fig5_lbfgs_trajectory`, Fig. 5 = `fig6_qas_fit`, Fig. 6 = `fig7_qas_nk`. Line numbers refer to `main.tex`.

**Bottom line.** The number bookkeeping is in good shape: every number I checked traces to a JSON, apart from a handful of rounding or selection slips. The physics is not yet in good shape. Three central claims are not supported by the data as they stand:
- that the polyimide film is "recovered in agreement with profilometry";
- that the six-knot model "removes" the negative κ_Si;
- that κ is intrinsically "hundreds of times less visible" than n.

Several other statements are scaling artefacts or grid artefacts. A cleaner and more physical account of the landscape is sitting in the JSONs (item 6), and the paper would be stronger for using it. **Verdict: major revision.**

---

## 1. Number audit

**Method.** `/tmp/audit/audit.py` extracts every number in the body of `main.tex` (abstract to backmatter, comments, citations and labels stripped). It searches all numeric leaves of `nb1–nb5.json` at the printed precision, with unit scalings ×10^{±2,±3,±6,±9}. That automatic match is permissive, so every non-trivial number in tables, captions and prose was then checked by hand against the JSON field it should come from.

**Result.** All table entries in Tables 1–3 match their JSON fields at the printed precision. That includes all 24 knot values, both stage columns, all Cramér–Rao (CR) rows, the softplus row and the loss values. The figure captions also match, except where items 1.1–1.10 say otherwise. Numbers that are not in any JSON:
- the literature constants n_Si = 3.4175 and n∥(633 nm) = 1.726 (the latter is in nb2.json);
- the derived value 42.6 = 1.726 × 24.7 (Table 1);
- the instrument metadata (1000 waveforms, 15 min).

These are acceptable, but per R10 the derived 42.6 should be written to nb2.json or dropped (see 2.4).

Inconsistencies found:

**1.1 [MINOR] Sec. 7, l.156: "Adam's end point, 3.7 % lower".** nb4: L-BFGS 4.6172e-3 against Adam 4.4522e-3. Adam is 3.57 % lower; L-BFGS is 3.71 % higher. Fix: "the L-BFGS end point lies 3.7 % above Adam's".

**1.2 [MINOR] Sec. 7, l.156: "the two end points differ by 5.7 µm in d₁ and d₂".** 5.7 µm is `valley/spread_d1_um_all_six`, the spread over all six end points. For the start-1 pair the difference is 5.31 µm (d₁) and 5.30 µm (d₂). Write "by 5.3 µm (5.7 µm over all six end points)", or say "up to 5.7 µm".

**1.3 [MINOR] Sec. 6.2(iv) uses κ₂ = −1.5×10⁻³ and 77; Table 2 and the Fig. 4 caption use fit 1 (−0.00144 → −1.4×10⁻³).** The CR calculation is done at the best point of the three runs, which is start 3 (κ₂ = −1.458×10⁻³). With fit 1 the ratio is 76. Use one point for the table, the bounds and the prose, and say which.

**1.4 [MINOR] Sec. 6.2(v): "raises the loss by only 0.29 %".** The JSON value of 0.2915 % is measured against the best loss over all three runs (start 3, 4.45158e-3). The softplus fit starts from start 1, and against fit 1 (4.4522e-3, the row printed directly above it in Table 2) the increase is 0.28 %. Compare like with like.

**1.5 [MINOR] Sec. 6.2(iii): "thicknesses at 31 nm".** The bounds are 31.6 nm (d₁) and 31.1 nm (d₂), and Table 2 prints 0.032 for d₁. Write "≈ 30 nm", or "32 and 31 nm".

**1.6 [MINOR] Sec. 6.1, l.140: "κ₁ is 280 times and κ₂ 890 times less visible".** The JSON gives 276.6 and 885.1. The rounding is acceptable, but see item 2.6: these ratios are artefacts of the scaling and should probably go.

**1.7 [MINOR] Sec. 7 text says "250–300 evaluations"; the Fig. 4 caption says "250–299".** Use one. Also, from nb4 histories, the stall level is reached at evaluation 181, 167 and 253. The rest are failed restarts. Say "reaches within 170–250 evaluations (250–299 including the restarts that confirm the stall)".

**1.8 [MINOR] Fig. 4 caption, l.160: "log-cosh (178–484)".** This silently drops log-cosh start 1, which hit the 400-iteration cap after 908 evaluations at d₁ = 987.7 µm, still moving. The writer's report knows this. Add "(start 1 not converged after 908 evaluations)".

**1.9 [MINOR] Table 3, Si κ at 0.76 THz printed as 0.0000.** The JSON value is −1.4×10⁻⁵, and `k_negative_fraction_in_band/Si` = 3.9 %. With a hinge penalty in place this is cosmetic, but print "−0.0000" or "<10⁻⁴" and do not claim the band is non-negative (item 2.12).

**1.10 [MINOR] Sec. 5, l.112: "the loss falls by an order of magnitude in the first 8 % of the iterations".** Per start the loss falls by a factor of 8.2, 11.8 and 5.3 (nb3 `loss_drop_phase1_percent` 87.7 / 91.5 / 81.1). Write "by a factor of 5–12".

**1.11 [MINOR] Sec. 9: "the spread between fits … here 2 %".** That holds for d₁ (977.7–1000.0 µm, 2.3 %). For d₂ the spread is 375.7–398.1 µm (5.8 %), and for n_Si it is 3.304–3.499. Write "2 % in d_q and 6 % in d_Si". This also endorses outside item 15.

**1.12 [MINOR] Unit inconsistency for optical thickness.** The Fig. 4 caption gives n·d in mm (1.9526); Table 3 gives it in µm (1926.0). Use µm throughout. In the caption, "n₁d₁ = 1.9526" also has no unit.

**1.13 [MINOR] Table 1 row "local minimum of two starts" (l.97, l.94).** The two starts reach *different* points:
- start 1: d = 28.47 µm, n = 1.611, κ_PI = 0.0107, +24 % loss;
- start 7: d = 25.54 µm, n = 1.724, κ_PI = 0.0106, +0.9 % loss.

Only the first is printed, and the caption implies both starts landed on it. Fix the caption, or list both. Start 7 is effectively on the valley floor of the best fit, so it is not a "local minimum" in any useful sense.

**1.14 [MINOR] Sec. 6.1, l.140: ρ(n₁,d₁) = 0.97 and ρ(n₂,d₂) = 0.99 are the cosine similarities of the Jacobian columns (`sensitivity/corr_*`).** Table 1 and Sec. 4 use ρ(n_PI,d_PI) = −0.93, which is the correlation of the *estimates* from (JᵀJ)⁻¹. These are two different quantities with opposite sign conventions. The estimate correlations for quartz/Si are −0.991 and −0.998 (`cramer_rao/rho_*`). Use the estimate correlation everywhere and define it once.

---

## 2. Physics and method

### 2a. Polyimide on silicon (Sec. 4, Table 1, abstract)

**2.1 [MAJOR] The PI fit returns n_Si = 3.259, which is −4.6 % from 3.4175 and about 900 CR-σ away (bound 1.7×10⁻⁴). The text never mentions it, and it invalidates the d_PI comparison as framed.** High-resistivity silicon has n = 3.4175 to about 10⁻⁴ across 0.5–4.5 THz, so 3.26 is physically impossible.

Every PI start, including the "local minima", conserves three timing quantities (I recomputed them from nb2.json):
- n_Si·d_Si ≈ 1304.3 µm, from the Si echo spacing of 8.70 ps;
- Σ(n−1)d = 922.5 µm, from the main-pulse delay relative to the reference;
- hence d_PI + d_Si = 425.8–426.3 µm.

If n_Si is set to 3.4175, the timing gives d_Si = 381.65 µm and then d_PI = 44.4 µm with n_PI ≈ 1.0, which is absurd. The data are compatible with a 25-µm film only through the impossible silicon index. The alternative is a reference/sample time-zero offset of about 65 fs: (n−1)d mismatch of 19.7 µm, computed with n_PI = n·d/d_prof = 1.79 and d_Si = 381.65. Such an offset is quite plausible for a reference taken 23 min before the sample (file timestamps 13:13 and 13:36).

Consequences:
- The "+4.6 % vs profilometer" agreement rests on compensating systematics.
- The 18-CR-σ figure is not interpretable.
- The thin-film claim in the abstract cannot stand as written.

Required: run CHECK C1. Then reframe the PI paragraph and abstract around what the timing *does* determine. These are n_PI·d_PI (from the film-echo/pulse-shape) and d_PI + d_Si (only if time zero is trusted). This endorses and sharpens outside item 5. Ask the authors for a micrometer thickness of the wafer (question Q1).

**2.2 [MAJOR] "18 times the Cramér–Rao bound" (l.94) is misinterpreted.** The bound assumes an exact model and white noise at σ = 1.35×10⁻³. Three things break that assumption:
1. The fit residual is 2.69σ in RMS, so the model error, not the noise, sets the scatter. At minimum the bound should be scaled by √(χ²/ν) ≈ 2.7, giving about 6.8 bound-units instead of 18.
2. The noise is not white (item 2.8).
3. The 0.22-µm spread of d_PI among the six "agreeing" starts (25.755–25.976 µm) is already 3.6× the CR bound. The optimiser stops at different points of the n–d valley, so the scatter is optimiser-limited, not noise-limited.

The profilometer uncertainty and the spatial offset between the stylus trace and the THz spot (spin-coated films vary by several % across a wafer) are not stated either. Fix: drop "eighteen times the CR bound", or replace it with a statement that includes model error and the profilometer uncertainty. Report the n·d agreement among the six starts instead: 44.22–44.31 µm, 0.2 %, compared with 0.9 % in n and d separately. That number is the real evidence for "good in n·d".

**2.3 [MAJOR] The simulated sweep and the PI experiment are not comparable as claimed ("its round trip of 0.30 ps sits at the threshold", R6).** They differ in four ways:
- **Contrast.** The echo that separates n from d scales with the product of the two internal Fresnel coefficients. That is 0.565 × 0.239 = 0.135 for the sweep (film 3.6 on substrate 2.213) and 0.263 × 0.333 = 0.088 for PI/Si. At equal round trip the PI echo is therefore about 35 % weaker and the PI threshold lies at a larger thickness. Outside item 9's proposed sentence ("set by the round trip, not by the index contrast") would be false; I contradict it.
- **Noise.** The sweep adds white noise to the sample only and uses the same noise-carrying reference to generate and to fit, so reference noise cancels. In the experiment the reference noise is independent. On top of that, the sweep converges at (0.98–1.00)σ, while the PI residual is 2.69σ, i.e. it is model-error dominated.
- **Optimiser.** Both use coarse-to-fine L-BFGS with rMSE and a hinge (good). Sec. 3 does not describe this protocol at all (item 2.15).
- **Starts.** The sweep's start box is centred on the truth.

Fix: run CHECK C11(c), a PI-like sweep, and put its prediction for d ≈ 25 µm, with the CR band, into Fig. 1(a) or its caption. Otherwise downgrade the sentence to a qualitative analogy.

**2.4 [MAJOR] κ assignment in the PI fit.** κ_PI is pinned at the hinge boundary (1.5×10⁻¹² to 2×10⁻¹²) in all six "good" starts. Meanwhile κ_Si = 6×10⁻⁴, which is 5–10× the HR-Si value [Dai 2004], and corr(κ_PI, κ_Si) = −0.989 (nb2 CR correlation matrix). The fit cannot tell which layer absorbs. The two "κ-driven local minima" have κ_PI ≈ 0.0106 and κ_Si = 0, which is the physically plausible assignment for a polyimide. Also, a CR bound printed for a parameter sitting on a constraint boundary ("κ_PI < 10⁻⁶ ± 0.0018") is not a valid 1σ.

Fix:
- say that the absorption split between film and substrate is undetermined;
- mark κ_PI as "at bound" and drop its CR entry;
- do not call start 7 a local minimum;
- add the "independent" n·d of 42.6 with a note that it mixes a visible index with a THz thickness (endorses outside item 5), or drop it.

### 2b. Quartz/Si landscape (Sec. 5–6)

**2.5 [MAJOR] The quartz/Si silicon index is also a model artefact by the paper's own argument, and the paper uses that argument selectively.** Sec. 5 calls n_Si = 3.4478 "within a percent of the literature value", yet that is +0.9 % and about 100 CR-σ from 3.4175. The converged least-squares fit (nb4, rMSE) gives n_Si = 3.304 (−3.3 %), d₂ = 398.1 µm and κ₂ = −0.0073. Every step (i)–(iv) of Sec. 6.2 applies equally to n_Si.

Fix: state that both silicon constants are biased, and treat the CR argument as a test of "noise vs systematics" for *every* parameter. The honest uncertainty is the spread across the n–d valley (item 1.11), not the CR row.

**2.6 [MAJOR] "κ₁ is 280 times and κ₂ 890 times less visible than n₁" and "the total optical path is constrained three orders of magnitude [thousands of times] better than the extinction coefficients" (abstract, Sec. 6.1, Sec. 9) are artefacts of the fractional scaling.** n and κ of a layer enter T(ω) only through ñ = n − iκ, and T is holomorphic in ñ. Therefore ∂y/∂κ = Hilbert[∂y/∂n] (a −i rotation at every positive frequency), and by Parseval ‖∂y/∂κ_j‖ = ‖∂y/∂n_j‖ exactly. The JSON confirms it to four digits:
- √G_nn/n₀ = 6.0468/1.964 = 3.0788;
- √G_κκ/κ₀ = 0.021859/0.0071 = 3.0787;
- likewise 4.3017/3.4 = 1.2652 against 0.0068320/0.0054 = 1.2652.

So the "visibility ratio" 276.6 is exactly the scale ratio n₁/κ₁ = 1.964/0.0071 = 276.6. Per unit *absolute* change, n and κ are equally visible. In absolute terms κ₁ is actually bounded better than n₁ (8.7×10⁻⁶ against 6.1×10⁻⁵). The κ's are poorly determined only *relatively*, because the samples are weakly absorbing (κ ≪ n). The same applies to the "factor of 5500" span of the singular values: most of the gap to the two κ directions is the choice of units.

Fix:
- delete the 280/890 sentence;
- in the abstract and conclusions write "relative to their size, the extinction coefficients are determined two to three orders of magnitude less precisely than the optical thicknesses, because the layers are weakly absorbing";
- state once that the singular values refer to coordinates scaled by the nominal values.

Also correct the label: 8.66 … 0.0016 are the singular values of J/√L (= √eig G), not of G. The eigenvalues of G span 3×10⁷ (endorses outside item 6).

**2.7 [MAJOR] The valley structure is mis-described, and a cleaner physical description is available in the JSONs.** I recomputed four quantities over all 16 quartz/Si end points in nb3/nb4 (Adam ×3, softplus, L-BFGS on Eq. 4 ×3, MAE ×3, rMSE ×3, log-cosh ×3):

| quantity | range over the 16 end points | relative spread |
|---|---|---|
| d₁ + d₂ | 1375.49–1375.85 µm | 0.03 % (0.013 % without the MAE start-1 point) |
| Σ(n−1)d | 1891.75–1892.02 µm | 0.014 % |
| n₁d₁ | 1952.27–1952.88 µm | 0.03 % |
| n₂d₂ | 1314.58–1315.33 µm | 0.06 % |

Over the same end points, d₁ ranges over 977.7–1000.0 µm and κ₂ over −0.0073 to +0.0044.

The physics is simple: the three timings (main-pulse delay against the reference, the quartz echo spacing and the Si echo spacing) fix n₁d₁, n₂d₂ and the total physical thickness d₁ + d₂. The single soft direction is the *transfer of physical thickness from one layer to the other at constant n₁d₁ and n₂d₂*. Only the interface Fresnel amplitudes lift it. This is exactly the 5.3-µm (L-BFGS), 16.6-µm (rMSE) and 2.3-µm (softplus) displacements, each with d₁ and d₂ moving equal and opposite.

Likewise κ₁d₁ + κ₂d₂ (total attenuation) is conserved to ±2.5 % (4.67–4.90; log-cosh start 1, unconverged, gives 4.25), while κ₂ changes sign: the κ valley is the attenuation split between the layers.

The current text, "the anti-phase combination (n₂ up, d₂ down) … lifted only by the echoes, whose Fresnel coefficients depend on n alone", is wrong in two respects:
- the n–d split of a single slab is lifted by *timing* (main-pulse delay measures (n−1)d, echo spacing measures n·d);
- the residual degeneracy is the inter-layer thickness transfer.

Fix: replace the 57-word sentence with this counting argument (shorter and more transparent; it also absorbs outside items 1 and 3). It also says what the method can and cannot measure, which is what the title promises. CHECK C6 quantifies it.

**2.8 [MAJOR] The noise model behind the CR argument is not white, and the paper does not say so.** I split the first 10 ps of each trace into 2-ps segments:
- **Reference_15min_nitrogen.** Segment standard deviations range from 1.0×10⁻⁴ to 2.9×10⁻³, with non-zero segment means up to 5.8×10⁻⁴. The first-difference standard deviation ("white" part) is 4.5–7.6×10⁻⁵ in quiet segments and 4.2×10⁻⁴ in the 8–10 ps burst.
- **Quartz_Si.** The quiet white level is ≈ 4.5×10⁻⁵, against σ = 2.9×10⁻⁴.
- **PI_silicon.** The quiet white level is ≈ 2×10⁻⁵, against σ = 1.35×10⁻³.

So the canonical σ is dominated by low-frequency bursts and offsets (the coder's "bursty" remark, quantified). The white component is 6–70× smaller. The noise has strong temporal correlation, which (JᵀJ)⁻¹σ² does not describe. THz-TDS noise is also known to grow near the pulse (amplitude and timing jitter), where the information is, so the pre-pulse window cannot represent it.

Step (i) of Sec. 6.2 must state the white-noise assumption and say that it fails here. The conclusion of step (iv) is still robust, but for a different reason than given: the residual RMS of 1.2×10⁻² is 41σ (nb3 `residual/rms`; the paper only quotes MAE = 15σ) and 4–250× any defensible noise estimate. It is model error, so *any* noise-only explanation of κ₂ < 0 fails. Write step (iv) that way, quote the RMS/σ, and add CHECK C5 to give the |κ₂|/bound ratio under three noise models. Under a residual-scaled (sandwich-type) bound the ratio may drop to about 2 (my estimate: 77/41 ≈ 1.9 with σ → residual RMS). That would make "77σ" misleading if left alone.

**2.9 [MAJOR] Step (vi), "the frequency-independent model absorbing the high-frequency amplitude excess of Fig. 2(c)", is a hypothesis presented as a deduction. The mechanism as written does not work, and competing explanations are not excluded.**

Why the mechanism fails: a negative frequency-independent κ_Si adds *gain* that grows with f (α ∝ κf). That worsens a high-frequency excess rather than absorbing it. What actually changes between the Adam (κ₂ < 0) and L-BFGS (κ₂ > 0) end points is the *split* of the attenuation (κ₁ 5.4×10⁻³ against 3.0×10⁻³, item 2.7). That split acts through the echo amplitude ratios, not through the spectral shape.

A negative κ_Si means the Si-side echoes are stronger than the model can produce with κ_Si ≥ 0. Several unmodelled effects do exactly that:
- **(a) an interface layer or air gap between quartz and Si.** It raises the internal reflectivity, and "bonded" is unexplained (l.112). The same plates were later re-mounted with a 126-µm gap (Sec. 8), so they were evidently not permanently bonded.
- **(b) a reference amplitude or time-zero mismatch** (cf. item 2.1).
- **(c) frequency-dependent optical constants**, the paper's hypothesis.

Required: CHECK C4 (interface gap) and C3 (six knots on the *same* quartz/Si data with λ_κ = 0). Rewrite (vi) to report what the checks show. If only (c) survives, the claim stands. If (a) or (b) remove the sign, the paper has a different and more useful message: unmodelled interface layers show up as negative κ.

**2.10 [MAJOR] Table 2 "fit 1" is not converged, so the "fitted κ₂" is a budget-dependent snapshot.** The evidence:
- Adam's loss still falls by 1 % between 10⁴ and 2×10⁴ iterations (4.496 to 4.452×10⁻³).
- The gradient norm at the Adam end is 0.035 (nb4), the same as at the stalled L-BFGS point.
- Fig. 3(d) shows κ₂ still moving at full speed through zero.
- The starts differ by 3–5 CR-σ (n₂ 1.45×10⁻³, d 0.15 µm, κ₂ 6×10⁻⁵), so "converge to the same point to four figures" (l.112) reflects parallel drift, not convergence.
- The reference notebook was run to 5×10⁴ iterations (coder §2.2).

Required: CHECK C8. Then either base Sec. 6.2 on a converged estimate, which is the natural choice because the CR bound is the bound for the least-squares estimator (that is the rMSE L-BFGS point: κ₂ = −0.0073, about 380 CR-σ), or state explicitly that Table 2 is the 2×10⁴-iteration iterate. Replace "converge to the same point" with "follow the same valley".

**2.11 [MINOR] Eq. (4) is effectively a pure MAE.** At the optimum the log-cosh term equals 6.3×10⁻⁵ against 4.5×10⁻³ for the MAE (nb4: log-cosh loss 6.31e-5, MAE 4.54e-3). With residuals ≪ 1, ln cosh r ≈ r²/2, so the log-cosh term is half the MSE and contributes about 1.4 % of the value and about 1 % of the gradient. The statement that the loss "tolerates the few large residuals near the main pulse" is therefore a property of the MAE term alone. The log-cosh behaviour also depends on the arbitrary units of the trace (it is not scale-invariant). Say this in one clause. It also explains why L-BFGS on Eq. (4) behaves exactly like L-BFGS on MAE (nb4: identical stall).

### 2c. L-BFGS vs Adam (Sec. 7)

**2.12 [MAJOR] "Twenty times cheaper" depends on Adam's untuned η = 10⁻⁴.** Counting evaluations is fair: one Adam iteration and one L-BFGS closure each cost one forward and one backward pass (4.1 vs 4.5 ms, nb3/nb4). I endorse outside item 12 and suggest saying so. But at η = 10⁻⁴ in the multiplicative coordinates, Adam's maximum step in d is 0.1 µm per iteration, so its Phase I speed is set by η. Required: CHECK C7. Report the factor as a range over η, or state "at the η used in the rest of the paper".

**2.13 [MINOR] "Kink, not barrier" is supported by the data, but the text misdescribes the evidence.**
- The "no decrease for steps of 10⁻⁸" is a *steepest-descent probe* (nb4 `steepest_descent_probe_at_lbfgs_end`), not the line search. Along −∇L the loss *rises*: +6.4×10⁻¹² at t = 10⁻⁸, while the gradient predicts −1.3×10⁻¹¹ (|g|² = 1.33×10⁻³). That one-sided mismatch is the signature of a non-differentiable point, and it is the evidence to quote.
- The "monotonic" straight line is sampled at only 29 points (Δα = 0.05). The left and right slopes at α = 0 differ by about 10 %, which shows the kink, but the loss decreases linearly toward Adam at −2.7×10⁻⁴ per unit α. The point is therefore not a minimum of any kind.

Cite Lewis & Overton, "Nonsmooth optimization via quasi-Newton methods", Math. Program. 141 (2013) 135, for BFGS line-search failure on non-smooth functions.

Also, "with the rMSE … L-BFGS converges without a stall" is not distinguished by the JSON. The rMSE runs carry the same `stop_reason` ("line search stalled … after 3 restarts"). The evidence for real convergence is that the three starts agree to seven digits (977.7087–977.7089 µm). Quote that, and the gradient norm at the rMSE end (CHECK C9).

**2.14 [MINOR] Fig. 4(d): L-BFGS points are projected onto a plane built from Adam's Phase II path, through Adam's optimum.**
- The L-BFGS iterates sit 5–17 µm away along the thickness-transfer direction.
- Their κ₁ is 0.0031 against Adam's 0.0054. In (d) L-BFGS moves almost entirely along PC2 (κ₁), not along Adam's κ₂ valley.
- The background contours therefore do not represent the loss at the L-BFGS points.

State this in the caption, and report the out-of-plane fraction (CHECK C10). "Stops at the kink before the valley" should become "stops after redistributing the attenuation from quartz to silicon (κ₁ 0.0071 → 0.0031) without entering Adam's κ₂ valley".

### 2d. Six-knot quartz–air–Si (Sec. 8)

**2.15 [MAJOR] The removal of the κ_Si < 0 artefact is not demonstrated (abstract, l.198 "with it in the model the silicon is lossless").** Four observations:
1. A hinge penalty λ_κ⟨max(0,−κ)⟩ with λ_κ = 1 is active in *both* stages. At κ₂ = −0.0014 it would add 1.4×10⁻³ to a loss of 3.6×10⁻³, so negative κ is effectively forbidden. A constrained model cannot be credited with removing a sign it forbids.
2. Stage 1 is the frequency-independent three-layer model, and it already ends at κ_Si = +2.4×10⁻⁴ with zero penalty (Table 3, "stage 1" column). On run 3 the non-dispersive model shows no artefact, so the improvement relative to Sec. 5 comes from the different mounting, the explicitly modelled gap or the hinge, not from dispersion.
3. The stage-2 Si κ knots are 0 to 1.2×10⁻³ (plus 4.8×10⁻³ at 3 THz). That is α ≈ 0.1 cm⁻¹ at 1 THz, 3–10× the HR-Si value, so not "lossless".
4. The knots have no uncertainties.

Required: CHECKS C2 and C3, plus C12 (CR bands). Rewrite the claim according to the result. The likely honest version: "with an explicit air gap and dispersive constants the silicon extinction is small and positive; whether dispersion or the gap is responsible is tested in …".

**2.16 [MAJOR] The band-edge knots are treated inconsistently.**
- The caption acknowledges the κ_Si bump at 3 THz ("spectrum reaches the noise floor").
- It does not mention the equally large n_Si drop at the same knot (3.3916 → 3.3818, −0.3 %).
- Nor the quartz-n upturn at 3 THz (1.9780 → 1.9917, +0.7 %).
- Nor the 0.2-THz knots, also at low SNR (κ_Si = 1.17×10⁻³, the largest in-band value).

Silicon is dispersion-free to about 10⁻⁴ over 0.2–3 THz [Dai 2004]. The extracted n_Si(f) varies by 0.0145 peak-to-peak (0.4 %), with a ripple of period about 1.1–1.2 THz in both n_Si and κ_Si (Fig. 6(b, d): maxima near 0.2 and 1.3 THz, minima near 0.76 and 2.1 THz). That period equals the free spectral range of the air gap, c/(2·126 µm) = 1.19 THz. The ripple is most likely gap-etalon mismodelling (wedge or tilt) leaking into the Si constants, i.e. the dispersive freedom has moved an artefact rather than removed it.

Fix:
- call the 0.2-THz and 3-THz knots unconstrained band edges in the caption and text;
- draw the literature n_Si = 3.4175 on Fig. 6(b), which is currently off-scale;
- mention the ripple period;
- run CHECK C13 (model selection with Si knots tied). If tying the Si knots costs little loss, present Si as non-dispersive.

**2.17 [MAJOR] "A loss scan against n_Si with d_Si re-optimised puts the valley floor near 3.37 (literature value: +4 % in loss)" is a grid artefact.** The nb5 scan code (`landscape()`) does not re-optimise. It evaluates a 41×41 grid of (uniform n offset, d), with all other 26 parameters fixed, and takes the minimum over the *grid* in d. Problems:
- The d step is 0.58 µm (0.15 %). The n·d valley is far narrower than that (its CR bound is about 0.05 µm), so the grid minimum misses the floor by up to 0.29 µm.
- As a result, the scan minima (3.26×10⁻³) sit 11 % *above* the fitted loss (2.93×10⁻³). A true profile must pass through the optimum.
- Only 6 of the 41 n values are stored (spacing 0.041), so "near 3.37" cannot be told apart from the fitted 3.392.
- "+4 %" is read at 3.412, not 3.4175.

Fix: run CHECK C12 (a true profile) or delete the sentence. R15 permits one sentence on the air gap only.

**2.18 [MAJOR] The air-gap "conditioning" claim is not backed by the number given (R15).** "A single sharp minimum at 126 µm, 5.5 times lower than at the nominal 50 µm" comes from a 1-D *slice* with everything else fixed. Conditioning is a curvature or bound comparison, not a slice. From the same slice, ±4 µm (3 %) raises the loss by 26–32 %. The quartz n–d grid scan shows a comparable rise for about 1.2 % in n, so the slice does not show "far better". Fix: quote the CR bound (or profile half-width) of d_air next to those of d_q and d_Si, and of the thickness-transfer direction, from the 27-parameter Jacobian (CHECK C12). Write that number to nb5.json, which keeps the claim compliant with R15.

**2.19 [MAJOR] Sample identity between Sec. 5 and Sec. 8 is inconsistent.** The text says it is the same quartz plate and Si wafer "re-mounted with an air gap", but:
- quartz n·d is 1952.6 µm in Sec. 5 against 1926.0 µm at 1 THz in Sec. 8. Even with the group index (n_g ≈ 1.975, n_g·d ≈ 1929 µm) the difference is 1.2 %, which is 150 fs in echo spacing and trivially resolvable;
- solid thickness d_q + d_Si is 1375.7 µm against 1365.5 µm;
- Si n·d is 1315.0 µm against 1319.3 µm.

Either the samples differ, or there is a time-axis or time-zero systematic of tens of fs, the same signature as item 2.1. Ask the authors (Q2). If they are the same plates, the paper must discuss this, because it bounds the *accuracy* of n·d at about 1 %, far above the 0.01 % *precision* the paper quotes.

**2.20 [MINOR] A free reference delay or amplitude is never considered.** Across items 2.1, 2.9 and 2.19, the dominant systematic of all three experiments appears to be reference/sample timing and amplitude. One sentence in Sec. 3 or Sec. 9 is needed: absolute thicknesses require time-zero reproducibility (Δ(d₁+d₂) = −cΔt, so 10 fs ≈ 3 µm). CHECK C6 quantifies how much the d bounds grow if τ is freed.

### 2e. Simulations, protocol and conclusions

**2.21 [MAJOR] "The converged loss is (0.98–1.00)σ throughout, so the errors are set by the noise" (l.86) is a non sequitur.** A loss at σ shows the model can represent the data; it says nothing about where on a flat valley the optimiser stopped.
- At 5 µm, draw 0, four starts end at d = 4.57–4.91 µm with loss/σ between 1.00035 and 1.00049. That is a ±3.5 % spread at indistinguishable loss, comparable to the reported 4.8 % error.
- All three draws at 5 µm are biased the same way: d −3 to −5 %, n +1 to +2 %, and even n·d −2 to −3 %.

So the errors at the thin end are at least partly optimiser- or stop-criterion-limited. Fix: run CHECK C11(a) and compare the median errors with the CR bound at each thickness. Only where they agree is "noise-limited" justified. Run C11(b), more draws, because three draws make the medians and bands very noisy (e.g. n error at 100 µm, 0.05 %, is larger than at 70 µm, 0.013 %).

**2.22 [MINOR] "Fail steeply below a round trip of about 0.24 ps (10 µm)" is not what the data show.** The local log-log slope of the median d error is 3.0 (5→7 µm), 1.6 (7→10 µm), 2.4 (10→15 µm) and 5.3 (15→20 µm). The steepest change is between 15 and 20 µm, and below 10 µm the error grows smoothly. "10 µm" is simply where the median d error crosses 1 % (nb2 `thresholds`). Write that: "the median error in d exceeds 1 % below 10 µm (round trip 0.24 ps)". Relate the round trip to the pulse duration: quote the reference FWHM so that "the film echo no longer leaves the main pulse" is quantitative.

**2.23 [MAJOR] Sec. 3 misdescribes the Sec. 4 protocol.** Sec. 3 states that until Sec. 8 all parameters are multiplicative factors p = s·p⁽⁰⁾. nb1 and nb2 instead use:
- the RobustExtractor scaling (n = n₀ + 0.1x, κ = κ₀ + 0.005x, d = d₀e^{0.1x});
- a coarse-to-fine low-pass schedule (0.25 → 0.5 → 1 THz → full band);
- a hinge penalty against gain (κ ≥ 0) in the sweep and the PI fit.

None of this appears in the paper, yet it is exactly what avoids the "plain full-band" collapse to 22.3 µm that Sec. 4 reports. Fix: one sentence in Sec. 3 ("simulations and the polyimide film: coarse-to-fine L-BFGS on progressively low-pass-filtered traces, κ ≥ 0 by hinge") and the word "coarse-to-fine" in Sec. 4. Without it the 22.3-µm sentence cannot be interpreted.

**2.24 [MINOR] Abstract and conclusions against the evidence.**
- "Recovered in agreement with profilometry": no (item 2.1; endorses outside item 8).
- "Constrains the total optical path three orders of magnitude better than the extinction coefficients": unit artefact (item 2.6).
- "Cramér–Rao bounds show that the negative silicon κ … is a model artefact, removed by a six-knot dispersive model": the first half is fine if rephrased as in item 2.8; the second half is not demonstrated (item 2.15).
- "L-BFGS … with a smooth loss, is the better default": supported for speed and reproducibility, but say that its end point (κ_Si = −0.0073, n_Si = 3.30) is no closer to physical values.

**2.25 [MINOR] Sec. 7, last sentence (R9): "a combined time- and frequency-domain loss adds no information, its Cramér–Rao bounds being those of the time trace".** This is true only for a full-band complex-spectrum loss under white noise (Parseval). Add "full-band, white noise", or keep the R9 sentence but qualify it.

**2.26 [MINOR] Forward model (Sec. 2).** I checked `Simulate.rts_batched` / `RTm_batched`. Each layer is a slab in n₀, cascaded by a Redheffer-type recursion through zero-thickness n₀ gaps. That is exact, because the characteristic matrix of a zero-thickness layer is the identity. A reader seeing Eqs. (1)–(2), single-slab coefficients in n₀, will think the method ignores the neighbouring layers. Add half a sentence ("exact, since zero-thickness ambient layers inserted between slabs leave the characteristic matrix unchanged").

---

## 3. Figures

The page is 13.4 cm text width; the figures are designed at 13.3 cm and set at 0.57–0.72 \linewidth.

**3.1 [MAJOR] Every figure is below the brief/paperstyle 7-pt floor at print size.** nb3 sets axis labels at 6.5 pt before scaling. At 0.72 \linewidth (Figs. 1–4) 7 pt prints at 5.0 pt and 6.5 pt at 4.7 pt. At 0.62 and 0.57 (Figs. 5, 6) they print at 4.3 and 4.0 pt. The PCA axis loadings of Figs. 3–4 and the top axis of Fig. 1(a) are unreadable on paper. Fix: re-export at the inclusion width (e.g. `figsize(0.72)` with 7-pt fonts, so they are included at 100 %). This needs no extra space.

**3.2 [MAJOR] Fig. 6 (n(f), κ(f)) and Table 3 have no uncertainties.** This is an extraction paper whose central message is what can be measured. Add CR bands, at least, around the curves (CHECK C12), the literature n_Si line, and a mark on the unconstrained edge knots (item 2.16).

**3.3 [MINOR] Fig. 1(a): no CR curve and no PI marker.** Add the CR 1σ curve per parameter (C11a), so "noise-limited" is visible, and a vertical line at the PI round trip of 0.295 ps (or the PI-like sweep from C11c). The bands are min–max over three draws. Say "range of three draws", or better, use more draws (C11b).

**3.4 [MINOR] Fig. 1(b): the caption omits that the residual is ×10, and there is no ±σ band.** The PI film echo is not visible, which is the point; say so.

**3.5 [MINOR] Fig. 4(d), panel title "Phase-II plane of Fig. 4 (shared)".** The figure *is* Fig. 4, and the Adam figure is Fig. 3. The figure number is hard-coded in nb4. Use "of the Adam fit" without a number.

**3.6 [MINOR] Fig. 4(e): the 3.7 % gap between the L-BFGS stall and Adam's end point, the main result of Sec. 7, is invisible on the log axis.** Add an inset (linear y, 4.4–4.7×10⁻³), or horizontal reference lines. Also add "one shade per start" (endorses outside item 24).

**3.7 [MINOR] Fig. 4(d): state that the contours are those of the Adam plane and that L-BFGS points lie off-plane (item 2.14).**

**3.8 [MINOR] Fig. 2(b): the ±σ band is invisible (σ = 2.9×10⁻⁴ on a ±0.1 axis).** I endorse outside item 25. Also quote the residual RMS (41σ) in the caption.

**3.9 [MINOR] Fig. 3: the caption should say that the axes are in fractional (nominal-scaled) coordinates** and that the Phase II plane is entirely κ₁/κ₂ (u = 0.97κ₂ + 0.24κ₁, v = −0.97κ₁ + 0.24κ₂). Readers otherwise look for n–d motion in (d).

---

## 4. Code spot checks (nb3, nb4, nb5)

- **nb3.** Loss = MAE + log-cosh (`quartzsi_model.loss_from_y`, `kind='mae_logcosh'`). Adam lr = 1e-4, N = 20000: the loop runs `range(1, n_iters)`, i.e. 19 999 steps, which is negligible. Starts are `STARTS` = (1.982/1.900/2.010, 1000/900/900 µm, κ₁ = 0.0071; Si 3.40/0.0054/400 µm), as in the Table 2 caption. Parametrisation u = p/p_start (multiplicative, as in Sec. 3). The softplus map is κ = κ₀·softplus(βθ)/softplus(β), β = 5, as stated. The sensitivity and CR are evaluated at the best point of all trajectories (start 3), with scale = mean initial guess (item 1.3). The "singular values" are √eig(G) (item 2.6). **Consistent with the paper apart from items 1.3 and 2.6.**
- **nb4.** Same `QuartzSiProblem`, same loss and the same u = p/p_start; torch L-BFGS, strong Wolfe, history 20, as stated. `max_iter = 1` per step with up to 3 history-clearing restarts on a stall (the restarts are not mentioned in the paper; one clause would help). One evaluation = one closure (forward + backward), the same cost as one Adam iteration. **Consistent.**
- **nb5.** Six knots at linspace(0.2, 3.0, 6) = 0.20/0.76/1.32/1.88/2.44/3.00 THz; cubic interpolation; held at the edge values outside the band; data are `ref_3.txt` / `quartz_air_silicon_15_min_3.txt` (run 3); Adam lr = 1e-4, 15k + 30k iterations; n_scale 0.1, k_scale 0.005; hinge λ = 1 in both stages (item 2.15). The air layer has n = 1, κ = 0 frozen and d free. Stage 2 is converged: loss 2.9326–2.9329×10⁻³ over the last 1000 iterations. **Consistent with the paper.** The scans, however, are not profiles (items 2.17, 2.18).
- **nb2 (beyond the brief).** The noise is added to the simulated sample only, at σ_ref; the measured reference is used noise-as-is for both generation and fit (item 2.3). The starts are centred on the truth.

---

## 5. Extra checks for the coder (each ≤ 10 min CPU)

Priority order: C1, C3, C2, C4, C8, C12, C5, C7, then the rest.

- **C1 (PI, ~2 min).** Refit PI/Si with the same coarse-to-fine multi-start (8 starts) under three model variants: (a) n_Si fixed at 3.4175 with d_Si free; (b) as (a) plus a free reference delay τ and amplitude factor a (y = a·F⁻¹[T·X·e^{−iωτ}]); (c) everything free plus τ and a. Report n_PI, d_PI, n·d_PI, d_Si, τ, a and rMSE/σ for each. Write to nb2.json.
- **C2 (run 3, ~8 min).** nb5 stage 1 and stage 2 with λ_κ = 0, otherwise identical (if needed, cut to 10k + 20k iterations and say so). Report the sign of κ_Si (knots and the stage-1 value) and the loss.
- **C3 (Sec. 5 data, ~8 min).** The six-knot model (same `dispersive_grid64`, two layers, no gap) on the quartz/Si bilayer with λ_κ = 0: stage 1 for 10k iterations, stage 2 for 20k. This is the direct test of "dispersion removes κ_Si < 0" on the sample where the artefact occurs.
- **C4 (Sec. 5 data, ~1 min).** A three-layer frequency-independent model with a free interface gap (n = 1, κ = 0; starts at d_gap = 1, 5, 20 µm), using rMSE L-BFGS (which converges in about 200 evaluations). Also try a free-index interface layer (n_g ∈ [1.4, 1.8], adhesive). Report d_gap, κ_Si, n_Si and the loss change.
- **C5 (~3 min).** |κ₂|/bound at the nb3 best point and at the rMSE optimum under three noise models: (a) white σ from quiet-segment first differences; (b) generalised Fisher (JᵀC⁻¹J)⁻¹ with C a Toeplitz matrix from the pre-pulse autocorrelation; (c) a residual-based sandwich bound H⁻¹(JᵀΣ_rJ)H⁻¹, with a block bootstrap of residuals (block ≈ 1 ps, 100 resamples, L-BFGS rMSE warm-started, about 1 s each). Also report the residual RMS/σ (41).
- **C6 (seconds).** Add columns ∂y/∂τ (= −dy/dt) and ∂y/∂a to J for quartz/Si and PI/Si, and report the CR bounds of d₁, d₂ (d_PI, d_Si) with and without τ and a. Also report the eigen-direction for the thickness transfer (d₁↑ d₂↓ at constant n₁d₁, n₂d₂), to support the rewrite in item 2.7.
- **C7 (~4 min).** Adam with η ∈ {3e-4, 1e-3, 3e-3} from the three starts, at most 6000 iterations each (4 ms per iteration). Report the iterations needed to reach (i) the L-BFGS stall level of 4.617×10⁻³ and (ii) loss 4.5×10⁻³, and whether any run diverges.
- **C8 (~4 min).** Continue Adam from the fit-1 state (or start 1) to 6×10⁴ iterations. Report the loss, κ₂, d₁ and d₂ at 2, 3, 4, 5 and 6×10⁴, and the gradient norm.
- **C9 (seconds).** At the L-BFGS (Eq. 4) end point: one-sided finite differences (L(x+tp)−L(x))/t for t = 10⁻⁹…10⁻⁵ along −g, along the last L-BFGS direction and along x_Adam − x_LBFGS. Count the residual samples that change sign within t ≤ 10⁻⁶. Also report the gradient norm at the rMSE and log-cosh end points.
- **C10 (seconds).** For each L-BFGS path, the fraction of ‖x_end − x_start‖² lying in each Adam plane, and the out-of-plane distance of the end point (for the Fig. 4 caption).
- **C11 (≤ 9 min in total).** (a) CR bounds of n_f, κ_f and d_f per thickness for the sweep geometry, plus the reference pulse FWHM (seconds). (b) 10 noise draws at 5, 7, 10, 15 and 20 µm (5 × 10 × 5 starts, about 5 min); also re-polish the best start of each fit with tight tolerances to separate stop-criterion scatter from noise scatter. (c) A PI-like sweep: film n = 1.72, κ = 0.01, Si substrate 3.4175 / 1×10⁻⁴ / 400 µm, the PI reference and σ = 1.35×10⁻³, d = 10, 15, 20, 25, 30, 40 µm, 5 draws (about 3 min).
- **C12 (~5 min).** The 27-parameter Jacobian at the nb5 stage-2 optimum. Report CR bounds for all knots and thicknesses (bands for Fig. 6, d_air against d_q and d_Si). Replace the grid scans by true profiles: n_Si offset over 3.36…3.43 (8 points) with all other parameters re-optimised (warm-started Adam, 1500 iterations, about 15 s each), and d_air over ±10 µm (6 points).
- **C13 (~5 min).** nb5 stage 2 with the Si knots tied (quartz dispersive only). Report the loss and MAE/σ against the full model, and AIC/BIC with the effective σ.
- **C14 (seconds).** A noise table for all six traces: standard deviation over 0–5, 0–10 and 5–10 ps; quiet-segment first-difference standard deviation; segment means; pre-pulse autocorrelation length. This backs the sentence requested in item 2.8.

## 6. Questions for the authors (to the director's open-questions list)

- **Q1.** What is the mechanical (micrometer) thickness of the PI substrate wafer, and of the quartz plate and Si wafer in Sec. 5 and Sec. 8? The timing analysis predicts d_PI + d_Si = 426.0 µm and d_q + d_Si = 1375.7 µm (Sec. 5) or 1365.5 µm (Sec. 8).
- **Q2.** Are the Sec. 5 and Sec. 8 samples physically the same plates? How was the Sec. 5 stack "bonded" (adhesive? optical contact? clamped?), given that it was later re-mounted with a gap?
- **Q3.** Was the PI measurement purged? Its file names lack "nitrogen". How long passed between the reference and sample scans (23 min from the headers), and is the delay-stage zero stable at the 10-fs level?
- **Q4.** Is the PI substrate high-resistivity float-zone Si (needed for 3.4175 to apply)?
- **Q5.** Where do the sweep parameters come from (film n = 3.6, κ = 0.0203; substrate 2.213/0.0125/496 µm)? A PI-like sweep (C11c) would serve the paper better.

## 7. Outside reviewer's report: endorse or contradict

- **Endorse:** 5 (and extend: n_Si is the bigger problem, item 2.1); 6 (and add the scaling artefact, item 2.6); 8 (the abstract over-claims); 12; 15 (and add 6 % for d_Si); 16 ("77 times its CR bound"); 22; 24; 25; 28.
- **Contradict:** 9. The proposed sentence ("the threshold is set by the round trip, not by the index contrast") is physically wrong; the echo amplitude scales with the interface Fresnel product (item 2.3).
- **Refine:** 1 and 3. The "n·d valley" is in fact the *thickness-transfer valley*, with n₁d₁, n₂d₂ and d₁ + d₂ all conserved (item 2.7). Name it that way, rather than only harmonising "optical path" and "optical thickness".
- **Neutral:** 23 (merging Figs. 3 and 4). Brief item 9 asks for an L-BFGS figure like the old Fig. 5; a combined figure would satisfy it, so the director decides.

## 8. Disagreements with rulings (for the log)

- **R12 (canonical σ = standard deviation of the first 10 ps).** I do not ask to overturn it; the paper needs *a* fixed convention. But the estimator is poor. The window contains bursts with non-zero means: the quartz/Si reference has segment standard deviations from 1×10⁻⁴ to 2.9×10⁻³, and the white component is about 5×10⁻⁵ (item 2.8). All σ-relative numbers (15σ, 77, 2.69σ, 18 CR-σ, 44 dB) can therefore move by a factor of 5–70 under another reasonable convention. Request: keep R12 for the bookkeeping, but (i) state in Sec. 3 that the pre-pulse noise is bursty and non-white and that σ is a conventional scale, not a white-noise level; (ii) base the CR argument on C5, not on σ alone.
- **R6 (PI as the experimental counterpart "expected marginal in n and d, good in n·d").** The data are consistent with the n·d half: n·d agrees to 0.2 % across six starts. But the fit only reaches a 25-µm film through an impossible n_Si (item 2.1). Unless C1 resolves this, the PI paragraph must report that the extraction exposes a substrate/time-zero systematic, rather than present the result as agreement with profilometry. This changes the framing in R6, not its intent (thin-film experimental data stay in the paper).
- **R15 (air gap: one sentence with the number from nb5.json).** Fine in form, but the number now used (the 5.5× slice ratio) does not support "far better conditioned". The sentence must quote a bound or profile width (C12) written to nb5.json; otherwise drop "far better conditioned".
- **R11 (L-BFGS reported as found).** Agree. Two additions: the η-dependence of "twenty times" (C7), and the steepest-descent probe as the kink evidence (item 2.13).
- **No conflict with the brief.** Nothing above asks for Bayesian seeding or a polynomial model, or moves Sec. 7 away from L-BFGS vs Adam. The requested text changes are roughly length-neutral: the counting argument of item 2.7 replaces longer prose, items 2.6 and 2.17 delete sentences, and the figure re-export (3.1) frees no space but costs none.

---

## Counts

[MAJOR] 20 · [MINOR] 29 · [CHECK] 14 (C1–C14) · total 63.

- MAJOR: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 2.10, 2.12, 2.15, 2.16, 2.17, 2.18, 2.19, 2.21, 2.23, 3.1, 3.2.
- MINOR: 1.1–1.14, 2.11, 2.13, 2.14, 2.20, 2.22, 2.24, 2.25, 2.26, 3.3–3.9.
