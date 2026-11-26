# Writer report — loop 0

**File:** `main.tex` (body assembled by `briefs/build_body2.py`, a writer-only helper). Compiles with `briefs/compile.sh` without errors or overfull boxes. **Page count: 8** (pdfinfo), with placeholder PNGs for figs 2–5 (2:1 aspect) and the coder's real `fig6_qas_fit`, `fig7_qas_nk`. Abstract 197 words, no numerals.

## Structure
1. Introduction (2 paragraphs) · 2. Forward model, eqs (1)–(3) · 3. Inverse problem: scaling, rMSE/MAE+log-cosh, Adam update eq, L-BFGS paragraph, one Bayesian sentence · 4. Simulated extractions and a thin film at the limit (sweep + polyimide/Si; **Fig. 1** = `fig2_thin_film`, **Table 1** PI/Si) · 5. Experimental quartz/Si (**Fig. 2** = `fig3_exp_fit`, **Table 2** fits + softplus row + CR row) · 6.1 Sensitivity matrix and the two-phase descent (**Fig. 3** = `fig4_adam_two_phase`, sensitivity panels merged in as (a,b)) · 6.2 Cramér–Rao step by step (i)–(vi), softplus in three lines, no figure · 7. L-BFGS vs Adam (**Fig. 4** = `fig5_lbfgs_trajectory`; landing points in the caption) · 8. Quartz–air–silicon, six knots (**Fig. 5** = `fig6_qas_fit`, **Fig. 6** = `fig7_qas_nk`, **Table 3** knots/thicknesses, numbers from `results/nb5.json`) · 9. Discussion and conclusions (≈0.4 page) · 16 references.

## What was cut
Old Figs 1, 3–7, 10–13 and Tables 1, 3, 4; the Bayesian subsection; the polynomial model and its coefficients table; the pair-slice and constrained-κ figures; the combined-loss section (one sentence remains, R9); all loss-function/Hessian digressions; ~70 % of the prose (20 pp → 8 pp). Citations reduced from 39 to 16 (unused keys left in `refs.bib`). `Greenall2017` dropped: flagged as unverifiable by the refs reviewer and its note broke the build.

## Deviations from the brief (all forced by the page budget)
- **Sec 4.1 and `fig1_five_layer` removed** (R4 foresees this); the five-layer result survives as one sentence in Sec 4 with `\nb{1}`.
- `fig_sensitivity` merged into Fig. 3 as panels (a,b) — flagged in the caption for the coder.
- Old Table 4 removed; the Adam/L-BFGS landing points are in the Fig. 4 caption (brief allowed removal).
- Sec 6.1 and 6.2 merged into one subsection; Sec 6.4 (softplus) folded into the Cramér–Rao subsection as three lines (brief: "a few lines, no figure").
- Table 2 shows fit 1, the softplus fit and the CR row; fits 2 and 3 are stated as identical to four figures. The three starts are in the caption.
- Figures are set at 0.64–0.82 `\linewidth`; with the coder's full-width designs the fonts will scale to ≈6.5–8 pt. **With real figs 2–5 (fig4 up to 9 cm tall) the paper will very probably exceed 8 pages again**; the next cut candidates are Table 1's secondary-minimum row, the Fig. 3/4 captions, and 2–3 references.
- Numbers from the old manuscript that nb3 will regenerate (Table 2, SNR, 90σ, 0.25 %) are kept but marked `\nb{3}` so the reviewer sees them.

## `\todo` numbers awaiting the coder
- **nb1**: five-layer errors (one `\nb{1}` sentence; no numbers printed).
- **nb2**: Fig. 1 caption (σ, SNR, error level, thickness threshold, round-trip threshold); Table 1 entire (best fit, CR bounds, secondary minimum, n·d, loss/σ).
- **nb3**: Fig. 2 caption (σ, residual/σ); Table 2 (all rows, old values shown); Fig. 3 caption (σ1/σ6, correlations, phase-split iteration, agreement %, variance captured, valley axis, loss gain, budget fraction); Sec 6.2 (κ2 bound, fitted κ2, 90σ, 0.25 %); Sec 5 (SNR, "first few percent").
- **nb4**: Fig. 4 caption (evaluations Adam/L-BFGS, d1, d2, n·d agreement); Sec 7 ("three orders of magnitude" spread/bound).
- **nb5**: none outstanding (Table 3, Fig. 5/6 captions and Sec 8 scans filled from `results/nb5.json`).

Author `\todo`s kept: author block, sample details (PI/Si, quartz/Si bonding), sample material, mechanical thickness, funding, acknowledgments, disclosures, data link.
