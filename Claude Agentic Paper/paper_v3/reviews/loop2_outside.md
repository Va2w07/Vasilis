# Outside reviewer report — loop 2

Reviewer persona: optics/photonics researcher who fits models to data, not a THz-TDS specialist.
Material read: `DECISIONS_LOG.md` (R1–R32), `reviews/loop1_outside.md`, `reviews/writer_report_loop1b.md`, `reviews/main_loop2.pdf` (8 pp., read page by page), `main.tex`.
Tags: [MAJOR] blocks comprehension, mis-states a result or breaks internal consistency; [MINOR] should be fixed, does not block; [LANGUAGE] wording only. Rewrites are drop-in replacements.

**Readability verdict.** Much improved. The counting argument in Sec. 6.1 is the best paragraph in the paper: I could follow why the timings fix n_j d_j and d_q + d_Si and what the two soft directions are. I can now follow the story from the abstract to the end of Sec. 8. The cost of the cut from ~3100 to ~1700 prose words is a telegraphic style. Sec. 3, step (iv) of Sec. 6.2, Sec. 7 para 2 and the opening of Sec. 8 are now readable only by someone who already knows the result. The conclusions (≈110 words) no longer summarise the paper. Three new internal inconsistencies have appeared:
- the "sixteen end points agree to 0.06 %" claim is contradicted by two rows of Table 2;
- Fig. 4(e) cites Eq. (4), which no longer exists;
- Fig. 4 says "Adam optimum", which contradicts R21.

Table 3 is printed after the reference list.

---

## Part 1 — Status of loop-1 items

| # | Topic | Status | Where / remark |
|---|---|---|---|
| 1 | Two valleys distinguished | Resolved | Sec. 6.1 para 1 names the thickness-transfer and κ valleys. "valley" is used in Sec. 5 before this definition (new item 14). |
| 2 | "Echo" defined | Resolved | Sec. 1 para 1 |
| 3 | Optical thickness one term | Resolved | Sec. 3 line 1. The gloss still varies ("n d valley", "quartz n d"; new item 27). |
| 4 | CR, rMSE, MAE, ρ, softplus before use | Resolved | CR: Fig. 1 caption, Sec. 6.2. rMSE: Sec. 3. MAE: below Eq. (3). ρ: Table 3. softplus: Sec. 6.2. |
| 5 | PI reads as failure | Resolved | Sec. 4 para 2 reframed around τ and n·d. Footnote d in Table 1. |
| 6 | Sec. 6.1 hardest passage | Resolved | Counting argument; s_k listed in the Fig. 3 caption. The holomorphy sentence is a new hard spot (item 16). |
| 7 | PCA plane in plain words | Partly | Sec. 6.1 para 3 gives the plain sentence. The acronym "PCA" is never expanded. The iteration labels on the paths are not explained (item 22). |
| 8 | Abstract over-claim / long sentence | Partly | The PI sentence and the long sentence are fixed. There is a new over-claim ("an order of magnitude" against "eight times" in the body), and the Adam two-phase result is missing (item 7). |
| 9 | Threshold vs index contrast | Withdrawn | D1 ruling accepted |
| 10 | Forward reference to Sec. 5 | Resolved | Sec. 4 "with the reference of Section 5" is explicit, which is acceptable. |
| 11 | Why L-BFGS in Sec. 4, Adam in Sec. 5 | **Not** | Sec. 3 now describes two different pipelines (variables, optimiser, loss) without saying why (item 8). |
| 12 | One Adam iteration = one evaluation | Resolved | Sec. 7 para 1 |
| 13 | Kink argument | Partly | Sec. 7 para 2 says "kink, not a barrier". The plain-word cause (the corner of \|r\| at zero residual) is missing, and the citation is mis-characterised (item 9). |
| 14 | Sec. 8 96-word sentence | Partly | The sentence is split, but it opens with a symbol, has no motivation sentence, and "held constant outside the band" has been cut (item 19). |
| 15 | Sec. 9 "2 %": which fits; dispersive summary | Partly | "Estimator spread" is named, but the dispersive model and the thin-film limit are absent from the conclusions (item 1). |
| 16 | σ with three meanings | Partly | Fixed in the prose. The legends of Fig. 1(a) and Fig. 6 still read "CR 1σ" (item 5). |
| 17 | Notation (subscripts, θ/x/s, J_ti, n_f, V-shape, SNR) | Partly | J_ki, n_f and the q/Si subscripts in the text are fixed. Figures still use 1/2 (stated in the caption, acceptable). There are new clashes of x, n₀, s and θ (item 5). SNR (Fig. 1) and peak/σ (Fig. 2) remain two conventions. |
| 18 | Sentences > 45 words | Mostly | Two remain, in Sec. 5 para 1 and Sec. 8 para 2 (item 28). |
| 19 | Colloquialisms | Resolved | — |
| 20 | Inverted / ambiguous constructions | Resolved | There are new ambiguous "it"/"this" (items 29–31). |
| 21 | Affiliations 4, 5 unassigned | Not | Author to-do; acceptable for now |
| 22 | PCA captions: shading, axes, markers | Resolved | Fig. 3 caption. Only the numeric path labels are missing (item 22). |
| 23 | Merge Figs. 3 and 4 | Withdrawn | D2 ruling accepted |
| 24 | Shades per start | Resolved | Fig. 4 caption |
| 25 | Fig. 2(b) invisible band | Resolved | Caption and in-panel text |
| 26 | Inset axis, "(arb. u.)" | Partly | Fig. 2 is fixed. Fig. 5(b) "residual" has no unit (item 24). |
| 27 | Fig. 6 panel letters | Resolved | — |
| 28 | Print size | Mostly | All figures are full width at ≥ 7 pt. There are label overlaps in Fig. 3(b) and in the Fig. 4(d) legend (items 22–23). |
| 29 | Fig. 4 caption data dump | Resolved | A new cryptic parenthesis has appeared (item 23). |
| 30 | Fig. 1 caption | Partly | n_f and the draws are stated, but the caption is now garbled and the link between the top axis and the threshold is not made (item 21). |
| 31 | Units of Fig. 3(d) vs (a)–(c) | Resolved | Caption |
| 32 | Table 1 footnotes, rMSE | Resolved | — |
| 33 | Table 2 header / softplus | Resolved | — |
| 34 | σ = 9.0e-4 vs 2.9e-4 unexplained | **Not** | Fig. 5 caption (item 20) |
| 35 | Number formats | Resolved | ×10⁻ⁿ throughout |
| 36 | Five-layer numbers | Resolved | One sentence |
| 37 | Instrument details to caption | Resolved | Fig. 2 caption |
| 38 | Sec. 6.1 numbers | Resolved | — |
| 39 | Sec. 6.2 numbers | Resolved | — |
| 40 | Sec. 7 numbers | Partly | 5.7 µm is gone. New argument-free numbers have appeared (item 33). |
| 41 | Sec. 8 numbers | Resolved | — |
| 42 | Merge candidates | Resolved / moot | — |
| 43 | THz-TDS sentence; geometry; contributions | Partly | (i) is done. (ii) is an author to-do. (iii) Sec. 1 para 2 lists topics, not contributions (item 12). |
| 44 | Duplicate labels | Resolved | — |

Summary: 28 resolved or mostly resolved, 11 partly, 3 not (11, 21, 34), 2 withdrawn by ruling (9, 23).

---

## Part 2 — New items

### A. Story, consistency and over/under-claiming

**1. [MAJOR] The conclusions no longer conclude the paper.** Sec. 9 has ≈110 words. It does not mention three of the four goals of the brief:
- the thin-film limit;
- the typical loss profile (the two-phase Adam descent along the κ valley);
- the six-knot model.

It also does not mirror the abstract. A reader who jumps from the abstract to Sec. 9 (as most will) finds a different paper. R32 deferred this "unless loop 2 finds room". Room can be found:
- Table 2 caption: drop ~40 words (item 4);
- Fig. 1 caption: drop ~25 words (item 21);
- Sec. 7: drop ~25 words (item 33);
- Table 3 currently floats after the references (item 2).

*Rewrite (≈190 words, replaces Sec. 9 except the to-do):*
"A differentiable TMM fitted to one transmitted pulse measures what the pulse timings measure: the optical thickness n_j d_j of each layer and, if time zero is trusted, the total thickness. For a thin film this holds down to a round trip of about two thirds of the pulse width, where the simulated errors are still at the Cramér–Rao bound. In the experiments the error budget is set by the model, not the noise. The residuals exceed the noise floor by more than an order of magnitude. A time offset between reference and sample shifts the total thickness by about 3 µm per 10 fs, and must be fitted. The same plates measured twice differ in quartz n·d by about 1 %, against a precision of 0.01 %. The loss landscape has two soft directions: a transfer of thickness between the layers and a transfer of attenuation between them. Adam completes the well-conditioned descent quickly and then spends most of its budget drifting along the second. L-BFGS completes the descent without a step size to tune, but stalls at a kink of the mean-absolute loss, so it is best used with a smooth loss. A realistic uncertainty is the spread between estimators (2 % in d_q, 6 % in d_Si). The sign of a weak absorption is set by a constraint, not by a richer model. A six-knot dispersive model leaves silicon non-dispersive and fixes an air gap eight times better than either plate."

**2. [MAJOR] Table 3 is printed after the reference list (p. 8).** An Optica copy-editor will return this. Change `\begin{table}[!b]` to `[!htb]` and move the environment above the "First, the gap…" paragraph, so that it lands on p. 7 or the top of p. 8. If it then pushes the paper to 9 pages, use the savings listed in item 1.

**3. [MAJOR] "All sixteen quartz/Si end points … agree on d_q + d_Si and both n_j d_j to 0.06 %" (Sec. 6.1 para 1) is contradicted by Table 2.** The products in Table 2 give n_q d_q ≈ 1952.3–1952.8 µm for the four two-layer rows. The "+ gap" and "+ adhesive" rows give 1931.4 and 1929.1 µm, which is 1.1 % lower; their n_Si d_Si is 0.7 % higher. Their d_q + d_Si is also 8–10 µm short, unless the gap is added in.

The reader cannot tell which sixteen end points are meant: Table 2 has six rows, and its caption mentions others. There is a second consequence. Sec. 9 attributes the 1.4 % difference in quartz n·d between Sec. 5 (1952.5) and Sec. 8 (1926.0) to time zero. But the Sec. 5 data fitted *with an interface gap* give 1931.4, only 0.3 % from Sec. 8. The difference may therefore be the interface model, not time zero.

*Request:* the director or domain reviewer should confirm before Sec. 9 keeps "likely set by time zero".

*Rewrite (Sec. 6.1):* "All sixteen two-layer end points of this paper (three Adam starts, softplus, L-BFGS on both losses, and the variants in the Table 2 caption) agree on d_q + d_Si and on both n_j d_j to 0.06 %, while d_q spans 2 % and κ_Si changes sign. Adding an interface layer (last two rows of Table 2) moves n_q d_q by 1 %."

*Rewrite (Sec. 9), if confirmed:* "…differ in quartz n·d by 1.4 %, of which an interface layer in the model accounts for 1.1 %; time-zero reproducibility and the interface model together limit the accuracy of n·d to about 1 %, against a precision of 0.01 %."

**4. [MAJOR] The κ_Si-sign thread is hard to follow across Secs. 5, 6.2 and 8.**
- *Why κ < 0 is a problem is never said in plain words.* Sec. 3's "hinge penalty against gain" is the only hint. Add to Sec. 5, after "κ_Si is negative": "a negative extinction coefficient means gain, which a passive silicon wafer cannot have".
- *Step (vi) and Sec. 8 appear to contradict each other.* Step (vi) says richer models keep the sign. Sec. 8 then says stage 1 "already gives κ_Si > 0 without the hinge". Yet Sec. 8 para 1 has just said that the hinge keeps κ ≥ 0. The reader does not know whether stage 1 ran with or without the hinge, or why the sign flips on the re-mounted stack.
  *Rewrite (Sec. 8 para 2, second sentence):* "Second, silicon is non-dispersive. Tying its knots costs only 0.7 % in loss, and the knot ripple follows the free spectral range of the gap, c/2d_air ≈ 1.19 THz, which points to an imperfectly modelled gap etalon. A control run without the hinge (λ_κ = 0) gives κ_Si > 0 already at stage 1. The positive sign on this stack therefore comes from the mounting, not from dispersion, as on the bilayer of Section 5 (step vi)."
- *Table 2 caption.* The κ_Si values of the six-knot runs (the evidence for step (vi)) sit only in this 117-word caption. Move them into two table rows, "Adam, rMSE | six knots" and "Adam, rMSE | six knots + gap 6.4 µm", with the knot ranges in the κ_Si and n_Si cells. Then cut the caption to the gap/adhesive definitions. This makes step (vi) checkable at a glance and saves ~40 caption words.

**5. [MAJOR] Symbol clashes introduced by the compression.**
- **x** has three meanings:
  - the reference pulse x(t) (Eq. 2);
  - the simulation optimiser variable (Sec. 3, "n = n₀ + 0.1x");
  - the nominal-scaled coordinate x_i = (p_i − p_i^opt)/p_i^(0) (Sec. 6.1).

  Rename the Sec. 3 variable to u: "the optimiser variables u enter as n = n_start + 0.1u, κ = κ_start + 0.005u, d = d_start e^{0.1u}".
- **n₀** is the ambient index in Sec. 2 and the start value in Sec. 3. Use n_start, κ_start, d_start as above.
- **s** is used both for the singular values s_k (Sec. 6.1) and for the knot scale s in "θ^init + sϑ" (Sec. 8). Rename the knot scale to Δ: "each knot value is θ^init + Δ ϑ with a trainable ϑ and a fixed scale Δ".
- **θ** is the ratio p_i/p_i^(0) in Sec. 3 and a knot value θ^init in Sec. 8. Write "n^init" in Sec. 8, or say "the same scaling as in Section 3".
- **σ**: the legends of Fig. 1(a) and Fig. 6 read "CR 1σ". Change them to "CR bound" (loop-1 item 16).
- **Eq. number**: the y-axis of Fig. 4(e) reads "loss, Eq. (4)". The loss is Eq. (3) since the merge.
- **Adam "optimum"**: the Fig. 4(d) legend ("+ Adam optimum") and caption ("plane through Adam's optimum") contradict Sec. 5 and R21 ("a budget iterate, not a converged optimum"). Change both to "Adam end point".

**6. [MAJOR] Step (iv) of Sec. 6.2 is the crux of brief item 7, and it is unreadable to an outsider.** "With the residual in place of σ, or a sandwich covariance, the ratio is below 2 at the Adam point and 4–10 at the least-squares point." Four things are unclear:
- "sandwich covariance" is undefined and has no reference;
- "the ratio" has no antecedent;
- "the least-squares point" is not identified (presumably the "L-BFGS, rMSE" row);
- "77 bounds" is an odd unit.

*Rewrite (iv):* "(iv) Measured in CR bounds, κ_Si lies 77 bounds below zero and n_Si about a hundred bounds above the literature value. But the bound assumes that the residual is noise, and the residual RMS (41σ) exceeds any defensible noise estimate. Noise can therefore neither explain the sign nor make it significant. If the bound is recomputed from the actual residuals, either by scaling it with the residual RMS or with the robust 'sandwich' covariance [ref.], |κ_Si| is less than two bounds from zero at the Adam point and 4–10 bounds at the rMSE (least-squares) point of Table 2."

Add a reference for the sandwich estimator, e.g. H. White, Econometrica 48, 817 (1980), or Aster et al. [8], if it covers it. Drop steps (ii) and (iii) as numbered steps. They are bookkeeping, not argument: fold them into (i) ("…evaluated at the Adam point, Table 2").

**7. [MAJOR] Abstract vs body.**
- (a) "determines the air gap **an order of magnitude** better than either plate". The body says "eight times" (27 vs 221 nm). Write "eight times".
- (b) The abstract never mentions Adam's two-phase behaviour. Adam completes the descent and then drifts along the κ valley for most of its budget. This is the "typical loss profile" the brief names as a main goal, and it is Fig. 3, the paper's central figure. This is an under-claim.
- (c) "A polyimide film on silicon is recovered only once a reference time offset is fitted" slightly over-claims: d_PI is still 2 % above profilometry and the text calls the errors systematic.
- (d) "fix each layer's optical thickness and the total thickness" omits the time-zero caveat that the PI sentence and Sec. 9 make.

*Rewrite (192 words, three numbers):*
"Terahertz time-domain spectroscopy can return the thickness and complex refractive index of each layer of a stratified sample from one transmitted pulse, but the inverse problem is unstable for thin or multiple layers. We fit a differentiable transfer-matrix model to the measured pulse with Adam and L-BFGS and ask what the data determine. Simulations with a measured reference pulse show noise-limited recovery of a thin film down to a round trip comparable with the pulse duration. For a polyimide film on silicon, a time offset between reference and sample must be fitted, so absolute thicknesses rest on time-zero reproducibility. On a quartz-on-silicon bilayer the pulse timings fix each layer's optical thickness, and the total thickness given time zero, to better than 0.1 %. Two soft directions remain: a transfer of thickness between the layers and a split of the total attenuation between them. Adam completes the well-conditioned descent quickly, then drifts along the attenuation direction for most of its budget. The sign of the silicon extinction coefficient therefore depends on the estimator; a softplus constraint, not a richer model, fixes it. L-BFGS needs no step-size tuning but stalls at a kink of the mean-absolute loss. A six-knot dispersive model of a quartz–air–silicon stack determines the air gap eight times better than either plate."

**8. [MAJOR] Sec. 3 describes two pipelines without saying why (loop-1 item 11, not done).** The simulations and the PI film use:
- x-variables;
- L-BFGS run coarse to fine;
- the rMSE with a hinge.

Quartz/Si uses ratio variables θ, Adam and the MAE + log-cosh loss. An outsider asks whether the results are comparable and why. One sentence suffices.

*Insert after the first sentence of Sec. 3:* "We use two pipelines: a robust one for the simulations and the polyimide film, and a plain one for the quartz/Si fit, whose slow phase (Section 6) we want to see and compare with L-BFGS (Section 7)."

Also replace "with a hinge penalty against gain" by "with a hinge penalty that forbids κ < 0 (gain)".

**9. [MAJOR] Sec. 7 kink paragraph: the plain-word cause is missing, and the citation may be mis-characterised.**
- **Missing cause.** "as two residuals change sign" is the only hint. The reader needs the mechanism: |r| has a corner at r = 0.
- **Citation.** "The known failure of BFGS on non-smooth functions [16]" is questionable. Lewis & Overton (2013) is best known for showing that full BFGS with an inexact line search is *surprisingly effective* on non-smooth functions. The failures they report concern the line search at kinks and the limited-memory variant. As written, the sentence attributes to [16] the opposite of its headline message. Please have the domain reviewer check.

*Rewrite:* "L-BFGS then stops 3.7 % above Adam's end point. The cause is a kink, not a barrier. The MAE term |ŷ_k − y_k| has a corner wherever a residual crosses zero, and at the stopping point a step along −∇L makes two residuals change sign. The loss then rises at once (one-sided slope +0.023, where the gradient predicts −0.036), and the line search fails. Yet the loss falls monotonically along the straight line to Adam's end point. Such line-search stalls at kinks are documented for (limited-memory) BFGS on non-smooth functions [16]."

### B. Concepts introduced before or without plain words

**10. [MINOR] Time-offset nuisance parameter (Sec. 4 para 2).**
- The physics is there, but the reader is not told what τ and a represent physically, or the sign convention of τ (is −70 fs a sample arriving early?).
- The step "10 fs is 3 µm of total thickness" depends on the counting argument of Sec. 6.1, which has not yet been given.

*Rewrite:* "We therefore let the reference be shifted in time by τ and scaled by a, the drift of the delay-stage zero and of the source power between the two scans, y = a F⁻¹[T F[x] e^{−iωτ}] (τ < 0: the sample trace arrives early). Because the timings fix every n_j d_j (Section 6.1), a 10 fs shift of time zero moves d_PI + d_Si by 3 µm."

**11. [MINOR] The PI paragraph's two closing inferences are compressed past the point of logic.**
- "the attenuation moves from substrate to film, so its split is undetermined": a move between two fits is not by itself a proof.
- "A PI-like simulation puts the 1 % limit at 10–15 µm, so the PI errors are systematic": this skips the middle step.

*Rewrite:* "The two fits put the attenuation in different layers at almost the same loss, so its split between film and substrate is not determined. A simulation with PI-like constants puts the 1 % limit in d at 10–15 µm; the 25 µm film lies above it, where the simulation recovers d to 0.16 %, so the 2 % found here is systematic, not noise."

**12. [MINOR] Sec. 1 para 2 lists topics, not contributions (loop-1 43 iii).**
*Rewrite of its last sentence:* "We show (i) how thin a film can be before thickness and index separate only at the noise limit (Section 4), (ii) what the loss landscape of a real bilayer fit looks like and what it does and does not determine (Section 6), (iii) how Adam and L-BFGS travel through it (Section 7), and (iv) what a dispersive model adds (Section 8)."
Also "fitting the whole trace [7] does not remove this": "this" is ambiguous. Use: "fitting the whole time trace [7] avoids phase unwrapping but not the non-convexity."

**13. [MINOR] Kink, band-edge knot and residual-aware bound are named in captions or the abstract but explained nowhere in the text.**
- Kink: item 9.
- Residual-aware bound: item 6.
- Band-edge knots: add to Sec. 8 para 1: "The outermost knots, at 0.2 and 3.0 THz, are pinned only by the weak spectrum at the band edges and are poorly constrained; Fig. 6 marks them hollow."

**14. [MINOR] Sec. 5 uses "follow the same valley" and "budget iterate on a plateau" before either concept is defined.** The Adam budget (2 × 10⁴ iterations) is never stated in the prose; it appears only in the title of Fig. 3(d).
*Rewrite (Sec. 5 para 1, replacing the clause after "from three starts"):* "…and fitted with Adam (η = 10⁻⁴) from three starts for 2 × 10⁴ iterations. The three end points differ only along a soft direction of the loss (Section 6.1; Table 2). The loss is then on a plateau: three times as many iterations change it by less than 0.02 %. We call this point the Adam end point; it is not a converged optimum."

**15. [MINOR] Phase II: the link to Sec. 3 is not made.** Sec. 3 sets it up well: Adam moves every parameter at rate η "however weakly the loss depends on it". Sec. 6.1 then says Adam drifts along the κ valley "at Adam's full step size" without connecting the two. "Descends an n–d bowl" introduces a new word ("bowl").
*Rewrite (Sec. 6.1, last two sentences):* "In Phase I, until every start has completed 95 % of its loss decrease, each descends the well-conditioned n–d directions. In Phase II, the remaining 92 % of the budget, all three drift along the κ valley. Because Adam normalises each step by the gradient size (Section 3), it keeps moving κ at the full rate η although the loss hardly changes."

**16. [MINOR] The holomorphy sentence (Sec. 6.1 para 2, 44 words) is the new hardest passage.**
*Rewrite:* "Because T depends on n_j and κ_j only through ñ_j = n_j − iκ_j, the pulse is exactly as sensitive to κ_j as to n_j: ‖∂y/∂κ_j‖ = ‖∂y/∂n_j‖. The extinction coefficients are two to three orders of magnitude less precise only relative to their own small values. Much of the four-decade span of s_k comes from scaling each parameter by its nominal value."

**17. [MINOR] Undefined acronyms and symbols.**
- PCA (Sec. 6.1): write "principal-component analysis (PCA)".
- FWHM (Sec. 4): write "full width at half maximum (FWHM)".
- SNR (Fig. 1 caption): say whether this is the power or the peak-field ratio, and use the same convention as "peak/σ" in Fig. 2.
- L-BFGS: write "limited-memory BFGS (L-BFGS)" at first use in Sec. 1.
- ⟨·⟩ in the Sec. 8 hinge: write "mean over the band".
- G (Fig. 3 caption only): acceptable.
- HR-Si in the Fig. 6 legend: spell it out in the caption as already done ("high-resistivity").
- Checked and fine: σ, s_k, ρ, τ, a, rMSE, MAE, CR and FSR (written out, no acronym used).

**18. [MINOR] "Quartz" vs fused silica.** n_q = 1.96 is the index of fused silica, not of crystalline quartz (n_o ≈ 2.1). Sec. 5 says "The quartz index is that of fused silica" with a to-do. An optics reader will stop here. Until the authors confirm the material, write "a fused-silica ('quartz') plate". If the plate is fused silica, rename "quartz" throughout, including the title of Sec. 8.

**19. [MINOR] Sec. 8 has no motivation and changes two things at once.** It introduces dispersion *and* a different sample (re-mounted, with a 126 µm air gap, and σ three times larger) without saying why.
*New opening sentence:* "To test whether the misfit of Section 5 comes from the missing frequency dependence, we let the optical constants disperse, using a later measurement of the same plate and wafer re-mounted with a deliberate air gap."
Then: "In this model, n_j(f) and κ_j(f) are natural cubic splines through six equally spaced knots over 0.2–3.0 THz, the band above the noise floor, and are held constant outside it…" This restores R32's cut clause in five words and avoids opening the sentence with a symbol.

**20. [MINOR] σ = 9.0 × 10⁻⁴ (Sec. 8) vs 2.9 × 10⁻⁴ (Sec. 5) for the same plates is unexplained (loop-1 34).** Add to the Fig. 5 caption: "(shorter acquisition than Fig. 2)", or whatever the cause is. Also say in Sec. 8 para 2 what "the loss profile" is: "the loss, minimised over all other parameters at fixed n_Si, varies by less than 0.05 % from 3.38 to 3.43".

### C. Figures and captions

**21. [MINOR] Fig. 1 caption is garbled.** "median, range (bands) of film errors, 13 noise draws at 5–20 µm, 3 above" does not parse. "Hollow: PI-like film (1.72/0.01 on Si…): n 0.043, d 0.16, κ 17 %" mixes units. The median d errors (2.9, 0.88, 0.12 %) duplicate what the plot shows.
*Rewrite:* "Fig. 1. Thin film on a substrate. (a) Simulated thickness sweep (film n_f = 3.6, κ = 0.0203; substrate n = 2.213, κ = 0.0125, d = 496 µm; σ = 1.58 × 10⁻³). Symbols: median relative error of the film parameters, over 13 noise draws from 5 to 20 µm and 3 draws above; bands: range; dashed: CR bound. The top axis gives the film round trip; the 1 % error in d is crossed between 0.17 and 0.24 ps. Hollow symbols, offset for clarity: a PI-like film (n = 1.72, κ = 0.01 on Si, five draws) at 25 µm. (b) Polyimide on Si fitted with a free time offset τ (σ = 1.35 × 10⁻³); residual and ±σ band magnified ×10."
In panel (b), mark where the PI film echo falls (it is the point of the paragraph); only "Si echo" is labelled.

**22. [MINOR] Fig. 3.**
- The numbers along the paths (10, 100, 1,000, 10,000) are iteration counts; the caption does not say so. In (b), "10" overlaps the start marker.
- (e): the dashed line "I | II it. 1,615" is clear, but the y-axis should read "loss, Eq. (3)" for consistency with Fig. 4(e).
- In (a)–(c), say what "+" marks: the Adam end point of that start, projected onto the plane.
*Add to the caption:* "Numbers on the paths: iteration."

**23. [MINOR] Fig. 4.**
- (d): the legend ("● start ★ end + Adam optimum") is clipped at the top-left and overlaps the axes. Move it below the plot or into the empty lower-right, and change "Adam optimum" to "Adam end point" (item 5).
- (e): the y-label "loss, Eq. (4)" should be Eq. (3). The inset label "linear, 10⁻³" should be "×10⁻³, linear scale". The main panel wastes its upper half (10⁰–10¹) behind the inset; set y_max ≈ 0.1 and place the inset there.
- The caption parenthesis "(end point 0.05 out of, 1.13 in plane)" is cryptic. Rewrite: "the L-BFGS end point lies 0.05 off this plane and 1.13 along it (fractional units), so the contours, computed in the plane, are those of Adam's fit."
- Name the colours: "L-BFGS orange, Adam blue (faint)".

**24. [MINOR] Fig. 5(b).** The y-axis "residual" needs "(arb. u.)". The grey ±σ band is invisible at this scale; say "(±σ band thinner than the line)" as in Fig. 2.

**25. [MINOR] Fig. 6.**
- Quartz is blue and Si orange. Elsewhere orange means "fit" (Figs. 1, 2, 5) and "L-BFGS" (Fig. 4). This is tolerable, but a neutral second colour for Si would avoid a false association.
- The caption says the bands are 6.4 times wider with the residual RMS, while Table 2 says 41-fold. One clause would stop the reader asking: "(the residual is closer to the noise here than in Section 5)".

**26. [MINOR] Tables.**
- *Table 1 caption (77 words) is hard to parse.*
  *Rewrite:* "Polyimide on silicon (d in µm, τ in fs). Rows: best of eight starts with free time offset τ and amplitude a; its CR bound (none for κ_Si, which sits at the bound κ ≥ 0); the fit without τ and a; independent values (ᵃwhite-light interferometry at 633 nm; ᵇprofilometer; ᶜRef. [14]; ᵈproduct of a and b, a visible index times a mechanical thickness, a guide only). With n_Si fixed at 3.4175: d_PI = 25.25, τ = −65.4, rMSE 1.70σ. The two starts within 0.5 % of the best rMSE agree on n_PI d_PI to 0.22 % (all eight: 3.0 %)."
- *Table 2*: see item 4 (two six-knot rows in, caption cut). Significant figures: d to 0.01 µm with a CR bound of 0.03 µm is justified only for comparing estimators. Keep it, but the caption could say "differences between rows, not the absolute values, are significant at this precision".
- *Table 3*: move it above the references (item 2). The row label "d_q, d_air, d_Si (µm)" with three numbers per cell is fine. The caption should say why n·d is quoted at 1 THz ("n·d varies with f in this model").

### D. Language (Optica copy-editing)

**27. [LANGUAGE] One gloss for optical thickness.** The text uses "n d valley" (Sec. 6.1, in quotation marks), "quartz n d" (Sec. 9), n_j d_j, "n·d" (Table 3 caption "n d"), and "n_PI d_PI". Use n_j d_j (or n_q d_q) in all math contexts. Drop the quoted gloss "(“n d”)" from Sec. 6.1: the name "thickness-transfer valley" is clear on its own. Also unify the sample names to "quartz/Si" (body) and "quartz–air–Si" (Sec. 8); the paper currently also uses "quartz-on-silicon" and "Quartz–air–silicon".

**28. [LANGUAGE] Sentences of 45 words or more.**
- Sec. 5 para 1, first sentence (~50 words): replace with the item-14 rewrite.
- Sec. 8 para 2, "Second, silicon is non-dispersive: …" (~45 words, three colons): replace with the item-4 rewrite.
- Sec. 6.2 (i)–(iii) (~45 words in one sentence): replace with item 6.
- The abstract has none.

**29. [LANGUAGE] Ambiguous pronouns and dangling phrases.**
- Sec. 7: "Yet it falls monotonically" → "Yet the loss falls monotonically along the straight line to Adam's end point".
- Sec. 6.2: "The sign, like n_Si, is a model artefact" → "The sign of κ_Si, like the 0.9 % excess of n_Si, is a model artefact".
- Sec. 4: "…lies 2 % above profilometry, also with n_Si fixed" → "…lies 2 % above the profilometer value, whether n_Si is free or fixed at the literature value".
- Sec. 1: "does not remove this" → item 12.
- Sec. 6.1: "The soft n–d direction left, s₄," → "The singular vector of s₄".

**30. [LANGUAGE] Telegraphic phrasing the copy-editor will query.**
- Sec. 3: "For quartz/Si it sees factors θ_i = p_i/p_i^(0) of the start" → "For quartz/Si the optimiser variables are the ratios θ_i = p_i/p_i^(0) of each parameter to its starting value".
- Sec. 3: "the log-cosh term … adds about a percent" → "adds about 1 % to the loss".
- Sec. 3: "(reference, for the sweep)" → "(of the reference trace for the simulated sweep, which has no sample trace)".
- Sec. 3: "Fits use several starts." → "Every fit is run from several starts and the lowest loss is kept."
- Sec. 4: "10 fs is 3 µm of total thickness" → item 10.
- Sec. 6.1: "The geometry is set by the sensitivity" → "The shape of the loss near its minimum is set by the sensitivity".
- Sec. 6.1: "lifted only by the interface Fresnel amplitudes" → "constrained only by the echo amplitudes, through the interface Fresnel coefficients".
- Sec. 6.2: "at 0.3 % in loss" → "at a cost of 0.3 % in loss".
- Sec. 7: "Adam needs 5000–5900 iterations" → "Adam needs 5000–5900 iterations to reach the same loss".
- Sec. 7: "crosses the Phase-I bowl in about thirty iterations" → "completes the Phase-I descent in about thirty iterations (each with a few line-search evaluations)".
- Sec. 9: "a weak absorption's sign" → "the sign of a weak absorption".
- Sec. 8: do not open a sentence with a symbol (item 19).
- Title: "multi-layer" → "multilayer" (Optica house style).

**31. [LANGUAGE] Tense and spelling.** Present tense is used for analysis and past for the measurements ("were measured"), which is consistent. British spelling is consistent (optimiser, artefact, parametrised, modelling, normalised); a grep finds no -ize/-yze/color/center forms in the prose.

**32. [MINOR] Data availability ("available from the authors on request").** Optica discourages this wording, and brief item 4 delivers notebooks that regenerate every figure. Use a repository DOI (e.g. Zenodo). The to-do is present; I flag it so that it is not forgotten.

### E. Numbers in prose that carry no argument (brief item 3)

**33. [MINOR] Sec. 7 has acquired new decorative numbers.**
- "memory 20": move to the Fig. 4 caption.
- "123–154 evaluations (250–299 with the confirming restarts)": keep one range in the prose; I suggest "about 150 evaluations (300 with confirming restarts)", with exact ranges in the caption.
- "one-sided slope +0.023, gradient prediction −0.036": these are evidence and can stay, but they need a unit or "(scaled units)".
- "gradient down six decades, starts agreeing to 0.2 nm": move to the caption.

Elsewhere the following can go to captions: the low-pass schedule "0.25, 0.5, 1 THz" in Sec. 3 (method detail), and "(η = 10⁻⁴)" plus "6 × 10⁴" in Sec. 5 (handled by the item-14 rewrite). The remaining prose numbers (0.06 %, 2 %, 41σ, 77, 0.9 %, 0.7 %, 0.8 %, eight times, 3.7 %, 5000 vs 150) all carry an argument and should stay.

---

## Disagreements with rulings

- **R32 (dispersive summary / fuller conclusions deferred "unless loop 2 finds room").** I disagree that this is optional. The conclusions are the part most readers read, and at present they omit the thin-film limit, the Adam two-phase profile and the six-knot result (item 1). Room can be made from captions (items 4, 21, 33); the paper does not need to grow.
- **R23/R29 (1.2–1.4 % quartz n·d difference attributed to time zero).** Table 2's interface-layer rows give n_q d_q within 0.3 % of the Sec. 8 value (item 3). Before the paper says "likely set by time zero", the domain reviewer should confirm that the interface model is not the larger cause.
- **R21 ("budget iterate, not a converged optimum").** I agree. Fig. 4 has not yet been brought into line with it (item 5).
- D1 and D2: I accept the rulings.

## Counts

[MAJOR] 9 (items 1–9) · [MINOR] 19 (items 10–26, 32, 33) · [LANGUAGE] 5 (items 27–31) · **total 33**.
Loop-1 items: 28 resolved or mostly resolved, 11 partly, 3 not resolved (11, 21, 34), 2 withdrawn by ruling (9, 23).
