# Outside reviewer report — loop 3 (final pass, corrections only)

Reviewer persona: optics researcher who fits models to data, not a THz specialist.
Material read: `DECISIONS_LOG.md` (R1–R37), `reviews/loop2_outside.md`, `reviews/writer_report_loop2.md`, `reviews/main_loop3.pdf` (8 pp., read page by page), `main.tex`. Spot checks were made against `results/nb5.json` (C16, sigma_note) and `reviews/loop2_domain.md`.
Tags: **[MUST]** is an error, an inconsistency, an undefined symbol or a binding ruling not carried out. **[SHOULD]** is a wording or readability fix. All rewrites are drop-in replacements.

**Verdict in one line.** All nine loop-2 MAJOR items are resolved. The paper reads as one story. Six small [MUST] corrections remain. Two of them are new inconsistencies between the conclusions and the body.

---

## Part 1 — Status of loop-2 items

| # | Topic | Status | Where / remark |
|---|---|---|---|
| 1 | Conclusions summarise the paper | Resolved | Sec. 9 covers the thin-film limit, the model-set error budget, interface accuracy, Adam drift, the L-BFGS kink, the estimator spread, the constraint and the six knots (≈130 words). Two of its statements disagree with the body (C2, C6 below). |
| 2 | Table 3 before the references | Resolved | p. 7 |
| 3 | "Sixteen end points" vs the gap rows | Resolved | Sec. 6.1 now reads "within the two-layer model … an interface layer moves n_q d_q by 1 %". Sec. 9 follows R33. |
| 4 | κ_Si thread | Resolved | "means gain, impossible for a passive wafer" is in Sec. 5. "A separate run without the hinge" is in Sec. 8. The six-knot rows are in Table 2. One leap remains: "the sign comes from the mounting" (C17). |
| 5 | Symbol clashes | Resolved | u, ξ, n_a, n_st, Δ, ϑ are in place. The legends now read "noise-limited 1σ", "loss" and "Adam end point". There are two new small clashes: footnote *a* vs amplitude *a* (C3), and index *k* (C10). Δ and ϑ are undefined (C5). |
| 6 | Step (iv) readability | Resolved | Sec. 6.2 (ii): the sandwich covariance is defined and cited [16], and the least-squares row is named. |
| 7 | Abstract vs body | Resolved | It says "eight times", mentions the Adam drift, and states the time-zero and interface caveats. |
| 8 | Two pipelines, why | Resolved | Sec. 3 para 1 |
| 9 | Kink cause, Lewis & Overton | Resolved | Sec. 7 para 2. The plain-word corner and the weak/strong-Wolfe point are both given. |
| 10 | τ, a meaning and sign | Resolved | Sec. 4 para 2 |
| 11 | PI closing inferences | Partly | The 25 µm / 0.16 % step is in. "The two fits put the attenuation in different layers" still does not say which two fits (C19). |
| 12 | Contributions (i)–(iv); "this" | Resolved | Sec. 1 para 2 |
| 13 | Kink, residual bound, band-edge knots explained | Partly | The kink and the residual bound are explained. Band-edge knots appear only in the Fig. 6 caption ("Hollow: unconstrained"). That is acceptable. |
| 14 | Sec. 5 budget, "end point" | Resolved | Sec. 5 para 1. "0.16 µm in d" is ambiguous (C9). |
| 15 | Phase II linked to Adam's normalised step | Resolved | Sec. 6.1 para 3 |
| 16 | Holomorphy sentence | Resolved | Sec. 6.1 para 2 |
| 17 | Acronyms | Mostly | PCA, FWHM, L-BFGS and peak/σ are fixed, and ⟨·⟩ is now "band mean". "HR-Si" in the Fig. 6 legend is still not spelled out (C16). |
| 18 | "Quartz" vs fused silica | Not (author) | The `\todo` remains in Sec. 5. It must be settled before submission. |
| 19 | Sec. 8 motivation, "constant outside" | Resolved | Sec. 8 para 1 |
| 20 | σ 9.0×10⁻⁴ vs 2.9×10⁻⁴; "loss profile" | **Not** | R37 asks for one clause in Sec. 8. It is absent from the text and from the Fig. 5 caption (C1). "Loss profile" is still not defined (C18). |
| 21 | Fig. 1 caption | Partly | The PI echo is now marked in (b) and the unit mixing has gone. "13 noise draws at 5–20 µm, 3 above" is still telegraphic (C14). |
| 22 | Fig. 3 labels | Resolved | It now has "path numbers: iteration", a "+" explanation and the y-label "loss". |
| 23 | Fig. 4 legend, inset, caption | Mostly | The labels, the inset label, the caption and the colour names are fixed. The main-panel y-range is unchanged (cosmetic). |
| 24 | Fig. 5(b) unit | Resolved | "(arb. u.)" is present. |
| 25 | Fig. 6 colour; 6.4× vs 41× | Not | Optional. No action is needed for submission. |
| 26 | Table captions | Partly | The Table 2 rows have been added. The Table 1 caption is still dense, with the a×b clash (C3). The Table 3 caption does not explain the "stage 1" and "CR bound" columns (C13). |
| 27 | One gloss for optical thickness; sample names | Partly | The quoted "n d" is gone. "Quartz n_j d_j" (Sec. 9) and "n d" (Table 3) remain. The sample appears under four names (C20). |
| 28 | Sentences > 45 words | Partly | Four remain: Sec. 3, 6.2 (iv), Sec. 7 kink, Sec. 8 para 1 (C11, C12, C15). |
| 29 | Pronouns | Mostly | "The soft n–d direction left" is unchanged (C10). A new "imposes it" in the abstract (C7). |
| 30 | Telegraphic phrasing | Mostly | "at 0.3 % in loss" is unchanged (C11). |
| 31 | Tense, British spelling | Resolved | A grep finds no -ize/-yze/color/center/behavior/modeling in the prose. |
| 32 | Data availability | Partly (author) | The notebooks are named. The DOI is still a `\todo`. |
| 33 | Decorative numbers in Sec. 7 | Resolved | These were moved to the Fig. 4 caption. |

**Summary.**
- MAJOR (9): all 9 resolved.
- MINOR (19): 8 resolved, 3 mostly, 5 partly, 3 not. Of the three not resolved, item 18 is an author decision, item 20 is a ruling not carried out, and item 25 is optional.
- LANGUAGE (5): 1 resolved, 2 mostly, 2 partly.

---

## Part 2 — Corrections

### [MUST]

**C1. [MUST] The R37 σ clause is missing (loop-2 item 20).** Sec. 8 and the Fig. 5 caption still give σ = 9.0×10⁻⁴ without saying why it is three times the Sec. 5 value for the same plates. The numbers below are from `nb5.json loop2_checks/sigma_note`.
*Replace the Fig. 5 caption opening* "Quartz–air–Si, six knots ($\sigma=\num{9.0e-4}$)." *with:*
"Quartz–air–Si, six knots ($\sigma=\num{9.0e-4}$, three times that of Fig.~2 because two noise bursts fall in the pre-pulse window of this 3000-waveform, \SI{115}{\pico\second} scan; its white-noise level, \num{2.7e-5}, is lower than in Fig.~2, \num{4.5e-5})."
(Alternatively, put the same clause at the end of Sec. 8 para 1, after "(Figs. 5, 6)".)

**C2. [MUST] Conclusions: the 1.2 % / 0.2 % comparison cannot be checked from the paper.** Sec. 9 says that quartz n·d "differs by 1.2 % between two mountings, but by 0.2 % once an interface gap is modelled (time zero then accounts for only ~7 fs)". The only run-3 value printed is the *phase* n d at 1 THz, 1926.0 µm (Table 3). From that value a reader computes 1952.5/1926.0 = **1.4 %** and 1931.4/1926.0 = **0.3 %**. The 1.2 %/0.2 % figures use the *group* value 1928.2 µm (`nb5.json C16_run3`, nd_group_1THz_um), which appears nowhere. The "7 fs" and the 2.2 µm it stands for also appear only in the conclusions. A result should not first appear in the conclusions.
*Fix, two edits:*
- (a) Table 3 caption: replace "$n\,d$ at \SI{1}{THz} 1926.0 (q), \SI{1319.3}{\micro\metre} (Si)" with "$n\,d$ at \SI{1}{THz} 1926.0 (q; group value $n_{\mathrm g}d$ = 1928.2, the quantity the echo timing measures) and \SI{1319.3}{\micro\metre} (Si)".
- (b) Add to the end of Sec. 8: "Fourth, the quartz group optical thickness is \SI{1.2}{\percent} below that of the clamped bilayer (Section~\ref{sec:exp}) but within \SI{0.2}{\percent} of its fit with an interface gap (Table~\ref{tab:exp}); the remaining \SI{2.2}{\micro\metre} in $\dq+\dS$ corresponds to only about \SI{7}{\femto\second} of time zero."
- Then shorten Sec. 9 to: "Quartz $\nq\dq$ differs by \SI{1.2}{\percent} between two mountings but by \SI{0.2}{\percent} once an interface gap is modelled, so the interface model, not time zero, sets the accuracy of $n_jd_j$: \SI{1}{\percent} against a precision of \SI{0.01}{\percent}."

**C3. [MUST] Table 1: footnote "$^d$ $a\times b$" clashes with the amplitude *a* in the same table** (column *a* = 1.024). A reader will multiply by 1.024.
*Replace* "$^d$$a\times b$, a guide only" *with* "$^d$product of the values marked a and b (visible index × mechanical thickness), a guide only".

**C4. [MUST] Reference 6 is typeset as "K. ichi Sawanaka".** BibTeX reads the lower-case "ichi" of "Ken-ichi" as a "von" particle.
*In `refs.bib` (Yasui2005) replace* `Ken-ichi Sawanaka` *with* `Sawanaka, Ken-ichi`. This prints "K.-i. Sawanaka".

**C5. [MUST] Sec. 8: Δ and ϑ are undefined** ("each knot is its start value plus $\Delta\,\vartheta$"). They are used nowhere else.
*Replace with:* "each knot is its start value plus a trainable variable $\vartheta$ times a fixed scale $\Delta$". Simpler alternative, with the symbols dropped: "each knot is fitted as an offset from its start value".

**C6. [MUST] The conclusions mis-state the thin-film result.** Sec. 9 says "Thin films reach the CR bound down to a round trip of two thirds of the pulse width". Sec. 4 says something different. The sweep is noise-limited (at the CR bound) at *every* thickness. Two thirds of the FWHM is where the median error in d crosses 1 %.
*Replace with:* "In simulation, thin-film errors stay at the CR bound and remain below \SI{1}{\percent} in $d$ down to a round trip of about two thirds of the pulse width (FWHM)."

### [SHOULD]

**C7. Abstract, "a softplus constraint imposes it".** The antecedent of "it" is two clauses back. Replace with "a softplus constraint imposes the physical sign, κ ≥ 0". This adds 4 words, and the abstract stays under 200.

**C8. Abstract and conclusions: the PI time-offset finding** is in the abstract but not in Sec. 9. Optional, if space allows, add to Sec. 9 after the first sentence: "Absolute thicknesses need the reference time offset to be fitted (\SI{10}{\femto\second} ≈ \SI{3}{\micro\metre})."

**C9. Sec. 5, "\SI{0.16}{\micro\metre} in $d$".** Say which thickness: "in $\dq$" (or "in $\dq$ and $\dS$", whichever the nb3 spread gives).

**C10. Sec. 6.1.**
- (a) "whose singular values $s_k$" reuses *k*, the time-sample index of $J_{ki}$ and $\hat y_k$. Write "whose singular values $s_1\ge\dots\ge s_6$".
- (b) "The soft $n$--$d$ direction left ($s_4$)" → "The remaining soft $n$--$d$ direction ($s_4$)".

**C11. Sec. 6.2.**
- (a) Item (iv) has 57 words, and "do not remove it" has no clear antecedent. Replace with: "(iv) Richer models do not remove the negative sign: six knots leave $\kS<0$ across the band, a gap makes it swing from $-0.013$ to $+0.005$, and a free $\tau$ gives $\kS>0$ only with $\nS=3.64$. We fit $\tau$ only where, as for PI, it cuts the misfit strongly and brings a fixed-index substrate into agreement."
- (b) "It drives $\kS$ to zero at \SI{0.3}{\percent} in loss" → "It drives $\kS$ to zero at a cost of \SI{0.3}{\percent} in loss".

**C12. Sec. 3, the robust-pipeline sentence (~60 words).** Split it: "In the robust one the optimiser variables $u$ enter as … about the start values. L-BFGS runs coarse to fine, on traces low-pass filtered at 0.25, 0.5 and \SI{1}{THz} and then on the full band; the loss is the root-mean-square error (rMSE) plus a hinge penalty $\sum_j\max(0,-\kap_j)$ (weight 1) against gain."

**C13. Table captions.**
- Table 2:
  - Replace "Six knots: Section~\ref{sec:disp}, no hinge, ranges over knots." with "Six knots: the spline model of Section~\ref{sec:disp} applied to these data, no hinge; ranges over the knots (quartz knots not listed)." As written, a reader thinks these rows are the Sec. 8 sample.
  - Replace "rMSE change: … (rows 5–8)" with "rMSE change relative to row 4: …".
- Table 3: add "Stage 1: knots tied (non-dispersive); CR bound: range over the knots."

**C14. Fig. 1 caption.** Replace "Median error (symbols), range (bands), 13 noise draws at \SIrange{5}{20}{\micro\metre}, 3 above; dashed, …" with "Symbols: median error; bands: range, over 13 noise draws (5–\SI{20}{\micro\metre}) or 3 (thicker films); dashed: …". Replace "Hollow, offset:" with "Hollow symbols (offset for clarity):".
Note: the caption's PI-like round trip (0.29 ps, n = 1.72) and the panel (b) label "+0.30 ps" (fitted n = 1.784) are both correct for different films, but they look like a slip. Optionally add "(b) measured film, round trip \SI{0.30}{\pico\second}".

**C15. Sec. 7.**
- (a) Kink sentence (55 words). Replace with: "L-BFGS stops \SI{3.7}{\percent} above Adam's end point, at a kink rather than a barrier. The MAE term $|\hat y_k-y_k|$ has a corner where a residual crosses zero; along $-\nabla\Lt$ two residuals change sign, so the loss rises (one-sided slope $+0.023$ in scaled units, against a predicted $-0.036$), yet it falls monotonically along the straight line to Adam's end point."
- (b) "A combined time--frequency loss adds no information (white noise, full band)." is opaque to an outsider. Replace with "Adding a frequency-domain term to the loss adds no information: for white noise over the full band its CR bounds are identical."

**C16. Fig. 6 caption.** It has no title sentence, and "HR-Si" in the legend is undefined. Replace the opening with: "Optical constants from the six-knot quartz–air–Si fit: (a), (b) $n(f)$ and (c), (d) $\kap(f)$ of quartz and Si; stage~1 dotted." Replace "the literature value (dashed)" with "the high-resistivity Si (HR-Si) literature value \cite{Dai2004} (dashed)".

**C17. Sec. 8, "A separate run without the hinge gives $\kS>0$ in stage 1: the sign comes from the mounting."** "The mounting" is a leap for an outsider. Replace with: "…in stage 1, so on this stack the positive sign comes from the modelled air gap, not from dispersion."

**C18. Sec. 8, "the loss profile is flat to \SI{0.05}{\percent}"** (loop-2 item 20, second half). Replace with "the loss, minimised over all other parameters at fixed $\nS$, varies by less than \SI{0.05}{\percent} over $\nS=3.38$--3.43, which includes the literature value."
Also in Sec. 8 para 1:
- "whether the misfit of Section 5 is missing dispersion" → "whether the misfit of Section 5 comes from missing dispersion";
- "with knots tied (stage 1)" → "with each parameter's knots tied to one value (stage 1, non-dispersive)".

**C19. Sec. 4, "The two fits put the attenuation in different layers"** → "The fits with and without $\tau$ put the attenuation in different layers (Table~\ref{tab:pi}), so its split between film and substrate is not determined."

**C20. Naming consistency.**
- (a) Sec. 9 "Quartz $n_jd_j$" → "Quartz $\nq\dq$".
- (b) Table 3 caption "$n\,d$" → "$n_jd_j$".
- (c) The sample is called "quartz/silicon bilayer" (abstract), "quartz on silicon" (Sec. 5 title), "quartz/Si" (body) and "Quartz–air–Si" (Sec. 8). Use "quartz/Si" and "quartz–air–Si" in the body, and keep the spelled-out "quartz/silicon" only in the abstract and the Sec. 5 title.
- (d) "THz-TDS" is defined in Sec. 1 and never used again. Drop the acronym.

**C21. Fig. 3 caption.**
- "loss $-$\SIrange{81}{92}{\percent}" and "loss $-$\SIrange{22}{29}{\percent}" read as negative ranges. Write "loss falls by 81–\SI{92}{\percent}" and "…by 22–\SI{29}{\percent}".
- The iteration count is "1615" in the caption but "1,615" in panels (d) and (e). Use one form (Optica: 1615).

**C22. References.**
- (a) Refs. 10 and 15 print "Systems 32, (2019)" with a stray comma, and Ref. 15 has no pages. Add `publisher = {Curran Associates}` to both entries, so that the bst prints "(Curran Associates, 2019)". Add pages to Li2018 only if the authors confirm them (refs_verification.md: dblp gives 6391–6401).
- (b) Ref. 11 prints "(2015). Int. Conf. on Learning Representations (ICLR), San Diego, 2015; arXiv:1412.6980", with the year twice. Make it `@inproceedings{Kingma2017, booktitle = {3rd International Conference on Learning Representations (ICLR)}, address = {San Diego}, year = {2015}, note = {arXiv:1412.6980}}`.

---

## Part 3 — Abstract and conclusions

**Abstract (≈190 words).**
- It is readable, and every claim has a home in the body: noise-limited thin film (Sec. 4), PI time offset (Sec. 4), 0.1 % within one interface model with a 1 % shift (Sec. 6.1, 0.06 %), the two soft directions (Sec. 6.1), Adam drift (Sec. 6.1), the estimator-dependent sign imposed by softplus (Sec. 6.2), the L-BFGS kink (Sec. 7), and "eight times" (Table 3: 221/27 = 8.2).
- There is no over-claim.
- There is one under-statement: the finding that the residual is 41σ (model-limited) is not in the abstract. It is optional.
- Only C7 is needed.

**Conclusions.**
- They now mirror the abstract and add the honest-uncertainty statement: 2 % in d_q (Table 2: 977.7–999.6, 2.2 %) and 6 % in d_Si (376.1–398.1, 5.7 %, two-layer rows).
- Two statements disagree with the body: the thin-film sentence (C6) and the 1.2/0.2 % comparison, which rests on an unprinted group value (C2).
- The style is telegraphic but acceptable.

## Part 4 — Can a non-specialist follow the story?

Yes, with effort in two places. The story runs as follows:
- the counting argument of Sec. 6.1 tells the reader what the pulse timings fix;
- the thin-film sweep shows where that stops working;
- the polyimide film shows the first practical trap, time zero;
- the quartz/Si fit shows the second, the model-limited residual and the two soft valleys;
- Adam's two phases and the L-BFGS kink show how optimisers behave in such a landscape;
- the six-knot section closes the loop on whether dispersion explains the misfit.

Each section now begins with its purpose, and every acronym but one is defined. The two harder places are:
- Sec. 6.2, where four steps and a sandwich covariance are packed into one paragraph;
- Sec. 8, where tied knots, stage 1/2, the hinge control run and the 1.19 THz ripple arrive in about 120 words.

Both are readable after C11 and C17/C18. The main remaining risk for an outside reader is not comprehension but trust: the conclusions quote one number not visible in the body (C2), and two notation slips would be noticed by a careful referee (C3, C5). Once the six [MUST] items are fixed, the paper tells a coherent story at a level an optics generalist can follow.

---

## Counts

- Loop-2 items (33):
  - MAJOR 9/9 resolved.
  - MINOR: 8 resolved, 3 mostly, 5 partly, 3 not (18 author, 20 ruling not carried out, 25 optional).
  - LANGUAGE: 1 resolved, 2 mostly, 2 partly.
- New corrections: **[MUST] 6** (C1–C6) · **[SHOULD] 16** (C7–C22) · **total 22**.
