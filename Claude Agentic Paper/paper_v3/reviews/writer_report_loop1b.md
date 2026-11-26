# Writer report — loop 1b (revision after the loop-1 reviews and checks C1–C14)

**`main.tex`: 8 pages** including the references. `briefs/compile.sh` reports no errors and no undefined references, and the log has 0 overfull boxes.

Every number was checked against `results/nb1–nb6.json` with a python script (listing in the session; the keys used are named below). The PI/Si numbers are those of `nb2 loop1_checks/fig2b_variant`; the variant with `n_Si` fixed is quoted as the cross-check.

Numbering in the PDF:
- **Figures:** Fig. 1 = `fig2_thin_film`, Fig. 2 = `fig3_exp_fit`, Fig. 3 = `fig4_adam_two_phase`, Fig. 4 = `fig5_lbfgs_trajectory`, Fig. 5 = `fig6_qas_fit`, Fig. 6 = `fig7_qas_nk`.
- **Tables:** Table 1 = PI/Si; Table 2 = quartz/Si end points by estimator and model, with the CR row; Table 3 = six-knot knots, thicknesses and CR bounds.
- **Equations:** Eqs. (1)–(3).

## Layout decisions needed to reach 8 pages (please check)
The new full-width figures (Figs. 2–7 at `\linewidth`, as instructed) plus the new content first gave 11 pages. These steps brought it to 8:
1. **`fig1_five_layer` dropped** (R4/R14). The five-layer result is now one sentence: 0.002 % in n and d, 0.3 % in κ.
2. **Captions set at full text width** with `\captionsetup{width=\linewidth}`. The class default is 4.25 in. This is a one-line template override, commented in the preamble.
3. **Eqs. (1)–(2) merged into one display.** I introduced g±_j = ½(n0/ñ ± ñ/n0); the content is identical. The Adam update is now inline.
4. **`\raggedbottom\allowdisplaybreaks`** replace `\flushbottom`.
5. **Singular-value mini-table moved into the Fig. 3 caption** (the R26 fallback).
6. **The C3/C4 six-knot rows of Table 2 moved into its caption.** The gap and adhesive rows stay in the table.
7. **Prose cut from about 3100 to about 1900 words.** Abstract: 195 words, two numbers (0.1 %, six-knot).

## Domain review — MAJOR items
- **2.1 (impossible n_Si in PI):** stated. The timing-invariant argument is given. τ and a are fitted as nuisance parameters: τ −70.1 fs, n_Si 3.430, d_PI 25.24 µm, rMSE 1.70σ against 2.69σ. The no-τ fit is kept as a table row, and the n_Si-fixed variant is in the caption.
- **2.2 (18 CR-σ):** deleted. Reported instead: n·d agreement (0.22 % over the 2 of 8 starts within 0.5 %; all 8: 3.0 %), and CR bounds with τ and a in J.
- **2.3 (sweep vs PI):** PI-like sweep (C11c) added. Fig. 1(a) hollow markers and caption: 25 µm recovered to 0.16 % in d. Text: the 1 % limit lies at 10–15 µm, so the PI errors are systematic. Outside item 9's sentence is not added (D1).
- **2.4 (κ assignment):** with τ the attenuation moves from the substrate to the film; the text says the split is undetermined. κ_Si is marked "at bound" with no CR entry. The "local minimum" row is removed. The independent n·d of 42.63 µm carries a visible-index × mechanical-thickness footnote.
- **2.5 (n_Si also biased):** Sec. 5 and step (iv) of Sec. 6.2 say n_Si is +0.9 % off, about 100 bounds; the CR argument is applied to n_Si too. Realistic uncertainty = spread between estimators (2 % in d_q, 6 % in d_Si; Sec. 9).
- **2.6 (κ visibility artefact):** the 280/890 sentence is deleted. The holomorphy argument is given, with relative-precision wording in the abstract and Sec. 6.1, and the nominal-scaling caveat. The values are now called "singular values s_k of J/√L".
- **2.7 (valley structure):** replaced by the counting argument (R17), naming the thickness-transfer valley (s₄) and the κ valley (s₆). The 16-end-point invariants are quoted (< 0.06 %; d_q spans 2 %; κ_Si changes sign).
- **2.8 (non-white noise):** Sec. 3 says σ is a bookkeeping scale and that the white part is 6–72× smaller (C14). Step (iv) rests on residual RMS = 41σ and the C5 bounds (below 2 at the Adam point, 4–10 at the least-squares point), with 77 kept alongside.
- **2.9 (step vi a hypothesis):** rewritten per R28 using C3, C4 and the τ variant. Six knots and an interface gap or adhesive do not remove the sign; free τ gives κ_Si > 0 only with n_Si 3.64. The remedy is the softplus constraint. "Dispersion makes Si lossless" is withdrawn everywhere.
- **2.10 (Adam not converged):** described as a "budget iterate on a plateau" (C8: < 0.02 % change to 6×10⁴ iterations). "Follow the same valley" replaces "converge to the same point". Table 2 uses the best (start-3) point throughout.
- **2.12 (η dependence):** evaluation counting stated. Ranges from C7: Adam 5000–5900 iterations at 1e-4 and 360–460 at 3e-3 (none diverge), against L-BFGS 123–154 evaluations (250–299 with restarts). The advantage is stated as "no step size to tune", with Phase I in about 30 iterations.
- **2.15 (six knots credited with removing the sign):** withdrawn. Sec. 8 says stage 1 already gives κ_Si > 0 without the hinge (C2), so the sign owes nothing to dispersion.
- **2.16 (band-edge knots, ripple):** edge knots are called unconstrained band edges in the Fig. 6 caption (n_Si drop, κ_Si rise, quartz n upturn; hollow markers). The ripple at the gap FSR (1.19 THz) is mentioned. Si is presented as non-dispersive (C13: +0.7 % loss).
- **2.17 (3.37 / +4 % grid artefact):** deleted. Replaced by the C12 profile: flat to 0.05 % over n_Si 3.38–3.43, with the literature value on the floor.
- **2.18 (air-gap conditioning):** the "5.5×" sentence is deleted. Table 3 and the text give CR 27 nm (gap) against 221/220 nm (plates), with ρ = −0.98 and d_q + d_Si bound 39 nm ("eight times better").
- **2.19 (sample identity):** per R29, same plates clamped without adhesive. Sec. 9: quartz n·d differs by more than a percent between the two measurements, so accuracy is about 1 % against a precision of 0.01 %, likely set by time zero.
- **2.21 ("noise-limited" non sequitur):** now justified by C11b (scatter = CR to within 25 %, no bias; CR curves in Fig. 1(a)).
- **2.23 (protocol misdescribed):** Sec. 3 describes the RobustExtractor scaling, the coarse-to-fine low-pass schedule, the hinge and multi-start.
- **3.1 (font size):** all figures are at `\linewidth` (exported for that width).
- **3.2 (no uncertainties in Sec. 8):** CR bands are in Fig. 6 (caption gives the 6.4× residual factor), and a CR column is added to Table 3.

## Outside review — MAJOR items
- **1 (two valleys):** named once and used consistently: thickness-transfer ("n d") valley and κ valley.
- **2 (echo):** defined in the Introduction.
- **3 (optical thickness):** defined in Sec. 3; "optical path" is no longer used.
- **4 (CR, rMSE, MAE, ρ, softplus defined before use):**
  - CR is spelled out at first use (Fig. 1 caption), then defined in Sec. 6.2.
  - rMSE and MAE are defined in Sec. 3.
  - ρ is defined as a correlation coefficient in the Table 3 caption. The old ρ(n,d) of PI is no longer used.
  - softplus is defined with its formula.
- **5 (PI reads as a failure):** reframed around τ and n·d, with a footnote on the visible-index product.
- **6 (Sec. 6.1 hardest passage):** the counting argument replaces it, and the singular values are in the Fig. 3 caption.
- **7 (PCA plane):** a one-sentence plain definition is given in Sec. 6.1.
- **8 (abstract over-claims, long sentence):** rewritten. No sentence exceeds 45 words; 195 words.
- **16 (σ meanings):** σ is the noise floor only. CR multiples are written "77 bounds"; table rows are labelled "CR bound"; singular values are s_k.
- **22 (PCA captions):** the shading, axes (variance fraction and composition, nominal-scaled, 1 = q, 2 = Si), curvature scaling and markers are explained in the Fig. 3 caption, and Fig. 4 refers to it.

## Rulings applied
- **R16–R31:** all applied.
- **Lewis & Overton (2013):** added to `refs.bib` as `Lewis2013` and verified on Crossref: Math. Program. 141(1–2), 135–163, doi 10.1007/s10107-012-0514-2.
- **Author note:** the planned PI re-measurement is recorded as a LaTeX comment below the PI paragraph.
- **Author `\todo` items kept:** author block, PI sample details and substrate type, quartz/Si nominal and mechanical thicknesses, sample material, mechanical thicknesses in Sec. 9, funding, acknowledgments, disclosures, data link.

## Not addressed, or reduced for space
- Five-layer figure: dropped (see above).
- Domain 1.8: "log-cosh start 1 not converged after 908 evaluations" was cut from the Fig. 4 caption for space.
- Outside 15: the conclusions no longer summarise the dispersive model beyond the accuracy statement.
- Outside 14 / R22: "held constant outside the band" (six-knot model) was cut. Fig. 6 shows it.
- Outside 26: figure-internal labels were not changed (coder's domain).
- Captions are dense. The Fig. 3 caption carries the singular-value list (about 80 words in one sentence of semicolon-separated items).
- Outside 17: q/Si subscripts are used in the text; the figure axis labels still use 1/2, and this is stated in the Fig. 3 caption.
