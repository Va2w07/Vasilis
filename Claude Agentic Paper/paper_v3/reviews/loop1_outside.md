# Outside reviewer report — loop 1

Reviewer persona: optics/photonics researcher who fits models to data, no THz-TDS experience.
Material read: `reviews/main_loop1.pdf` (8 pp.), `main.tex`, `DECISIONS_LOG.md`, figure PNGs at full resolution.
Tags: [MAJOR] blocks comprehension or mis-states a result; [MINOR] should be fixed, does not block; [LANGUAGE] wording only.
Where a rewrite is offered it is a drop-in replacement for the quoted sentence.

Readability verdict in one line: the paper is honest, tightly written and the story (method → limit → landscape → optimiser behaviour → dispersive fix) is there, but it is written for someone who already knows what an echo, an n·d valley and a sensitivity matrix are. An outsider loses the thread in Sec. 6.1 and never fully recovers it until Sec. 9.

---

## A. Clarity of the story (terms used before definition, jumps)

**1. [MAJOR] Two different "valleys" are never distinguished.** Abstract and Sec. 6.1/Fig. 3 speak of a "flat κ valley"; Sec. 7 ends with "yet another point of the n d valley"; Sec. 9 uses "the per-layer n–d split". A reader who is not a THz specialist cannot tell whether these are the same soft direction or two. They are two (anti-phase n–d combination that preserves optical thickness, and the κ direction), and this is the central result of the paper.
*Location:* Sec. 6.1, first paragraph, after the singular-value sentence.
*Rewrite:* add one sentence: "The landscape therefore has two soft directions: the anti-phase n–d combination of each layer, along which the optical thickness n·d and hence the pulse delay are unchanged (the n·d valley), and the two extinction coefficients (the κ valley); the rest of the paper refers to them by these names."
Then use "n·d valley" and "κ valley" consistently (Sec. 7 "n\,d valley" → "n·d valley"; Sec. 9 "the per-layer n–d split only by the echoes" → "the n·d valley is lifted only by the echoes").

**2. [MAJOR] "Echo" / "Fabry–Pérot echoes" is used from the abstract onward and never explained.** For an outsider "echo" is the key physical idea (the delayed replica of the pulse from internal reflection is what separates n from d), and it is first explained implicitly only in Sec. 4 ("once the film echo no longer leaves the main pulse").
*Location:* Sec. 1, first paragraph, "for thin samples or several layers internal reflections overlap".
*Rewrite:* "for thin samples or several layers the delayed replicas of the pulse produced by internal reflections (the Fabry–Pérot echoes) overlap with the main pulse, the transfer function admits several solutions and the objective becomes non-convex [5], so thickness and index can no longer be separated, precisely in the thickness-gauging regime where THz-TDS is most in demand [6]."

**3. [MAJOR] "Optical path", "optical thickness", "n·d", "total optical path", "n₁d₁" are used as synonyms without ever saying so.** Abstract: "total optical path"; Sec. 6.1: "optical path", "optical-thickness degeneracy", "n₂d₂"; Table 3: "n d quartz, Si (µm)". Define once and pick one term.
*Location:* Sec. 3, after Eq. (4) or at the first use in Sec. 4.
*Rewrite (Sec. 4, para 2, last sentence):* "…so the errors are set by the noise, and the echoes [5, 6] are what separates n from d: the main pulse fixes only the optical thickness n·d of each layer (its delay), and the echo amplitudes, through the Fresnel coefficients, fix n alone." Then use "optical thickness n·d" everywhere (drop "optical path").

**4. [MAJOR] Cramér–Rao bound, rMSE, MAE, ρ(n,d) and softplus appear before they are defined.**
- "Cramér–Rao" first appears in Sec. 4 (Table 1 row, "eighteen times the Cramér–Rao bound") and is defined in Sec. 6.2. Add a forward reference in Sec. 4: "…and a discrepancy eighteen times the noise-limited (Cramér–Rao) uncertainty of Section 6.2."
- "rMSE/σ" (Table 1 header) and "rMSE" (Sec. 7) — the acronym is never introduced; Sec. 3 says "root-mean-square error". Write "root-mean-square error (rMSE)" in Sec. 3.
- "MAE" (Fig. 5 caption, Table 3) is never expanded. Sec. 5 says "mean absolute residual". Write "mean absolute error (MAE, the first term of Eq. (4))" in Sec. 5 and use MAE thereafter.
- "ρ(n_PI, d_PI) = −0.93" (Sec. 4) — ρ is never defined. Rewrite: "with a correlation coefficient ρ(n_PI, d_PI) = −0.93 between the two estimates".
- "softplus" appears in the Table 2 caption and row (Sec. 5) before Sec. 6.2 defines it. Caption rewrite: "from start 1 with κ constrained to be non-negative (softplus map, Section 6.2)".

**5. [MAJOR] The polyimide result as written reads as a failure in n·d, not a success.** The prose says the extraction "should be … good in n·d, and so it is", but Table 1 shows n·d fit = 44.26 µm against "independent" 42.6 µm (3.9 % off, hardly better than the 4.6 % in d). The reader is not told that the independent n is a visible-wavelength value, so the independent n·d is not a THz quantity and the comparison is not meaningful. Either remove the "independent" n·d entry and say so, or add a sentence.
*Location:* Sec. 4, para 3.
*Rewrite:* "…so the extraction should be marginal in n and d separately and reliable only in their product n·d, and so it is: six of eight starts agree on n·d to within 1 % at a loss of 2.7σ, but they split it so that d_PI lies 4.6 % above the profilometer value and n_PI 0.8 % below the visible-band index n_∥(633 nm) (no independent THz value of n_PI or n·d exists, so the last column of Table 1 is a guide only), with ρ(n_PI, d_PI) = −0.93 and a discrepancy in d eighteen times the noise-limited uncertainty."
Also: Table 1 shows n_Si = 3.259 against 3.4175 (−4.6 %) without a word of comment, while Sec. 5 finds 3.45 for the same material. An outsider will trip over this. One clause is needed ("the substrate index is likewise pulled along the n·d valley, by −4.6 %").

**6. [MAJOR] Section 6.1, first paragraph, is the hardest passage in the paper and is where a non-specialist stops following.** Six singular values in prose, "eigenvectors … singular values" of the same matrix, "after two directions that redistribute optical path comes a gap of an order of magnitude to the anti-phase combination (n₂ up, d₂ down) that leaves n₂d₂ unchanged, the optical-thickness degeneracy (ρ(n₁,d₁)=0.97, ρ(n₂,d₂)=0.99), lifted only by the echoes, whose Fresnel coefficients depend on n alone" (57 words, four nested ideas). The brief says numbers go in tables/captions; R15 allows the sensitivity panel "if it fits".
*Suggestion:* put the six eigenvalues and their dominant parameter combinations in a six-row mini-table (two columns: eigenvalue, direction), or as a sixth panel of Fig. 3, and reduce the prose to:
"G has six eigen-directions whose eigenvalues span a factor of 5500 (Table X). The stiffest is the total optical thickness, fixed by the pulse delay; the two next stiffest redistribute optical thickness between the layers. An order of magnitude below them lies the anti-phase combination of n and d within a layer, which leaves n·d unchanged and is lifted only by the echo amplitudes; another order of magnitude below lie the two extinction coefficients, which the pulse sees 280 and 890 times less clearly than n₁ for the same fractional change."
Also decide whether 8.66 … 0.0016 are eigenvalues of G or singular values of J/√L (they differ by a square; the "three orders of magnitude" of the abstract depends on which).

**7. [MAJOR] The PCA-plane construction needs one plain sentence before the jargon.** "the plane the trajectory itself chooses [13]: the two leading principal components of the centred iterates span it" is opaque to anyone who has not read Li et al.
*Location:* Sec. 6.1, para 2, first sentence.
*Rewrite:* "To see how an optimiser moves through this landscape we follow Li et al. [13]: the sequence of parameter vectors visited by the optimiser is treated as a point cloud, its two directions of largest spread (principal components, PC1 and PC2) define a plane, the loss is evaluated on a grid in that plane and the path is projected onto it. The axis labels of Figs. 3 and 4 give, for each PC, the fraction of the path variance it carries and its composition in the scaled parameters."

**8. [MAJOR] The abstract over-claims the polyimide result and is at the 200-word limit.** "a polyimide film on silicon at that threshold is recovered in agreement with profilometry" — the body says d is 4.6 % high and eighteen bounds off, i.e. marginal, as predicted.
*Rewrite:* "…and a polyimide film on silicon at that threshold behaves as predicted: its optical thickness is reproducible, its thickness alone is recovered to within 5 % of profilometry."
The 74-word fourth sentence should be split: "On a quartz-on-silicon bilayer the sensitivity matrix of the forward model constrains the total optical thickness three orders of magnitude better than the extinction coefficients and separates n from d only through the Fabry–Pérot echoes. The optimiser trajectory, projected onto the principal planes of its own path, shows Adam creeping along a flat κ valley for most of its budget, whereas L-BFGS completes the well-conditioned part of the descent in a few hundred evaluations." (Word count stays ≈200; trim "in principle" and "of every layer" in sentence 2 to compensate.)

**9. [MINOR] Why does the simulated sweep transfer to the polyimide film?** The sweep uses film n = 3.6 on a substrate n = 2.213 (high-index film on low-index substrate); the experiment is n = 1.7 on n = 3.4 (the reverse contrast). The paper transfers the threshold purely through the round-trip time. One clause: "the threshold is set by the round trip 2n_f d/c, not by the index contrast, so it applies to the polyimide film" (if that is what the authors believe; otherwise say the contrast differs).

**10. [MINOR] Sec. 4 uses "the measured reference of Section 5" before Sec. 5 exists.** Rewrite: "using the measured reference pulse of the quartz/Si experiment (Section 5) and white noise at its σ".

**11. [MINOR] Why L-BFGS in Sec. 4 but Adam in Sec. 5?** The simulations and the polyimide film use multi-start L-BFGS; the quartz/Si fit uses Adam; Sec. 7 then compares them. A reader asks why. Add to Sec. 5: "fitted … with Adam (so that the slow phase of the descent, Section 6, is visible; L-BFGS follows in Section 7)".

**12. [MINOR] Sec. 7: say explicitly that one Adam iteration is one function (and gradient) evaluation**, otherwise the "twenty times cheaper" and Fig. 4(e) axis "function evaluations" vs Fig. 3(e) "iteration" look like different currencies. Rewrite: "reaches, in 250–300 evaluations, the loss that Adam reaches after 5000–5900 iterations (one evaluation each); then it stops."

**13. [MINOR] Sec. 7, the kink argument.** "yet the loss decreases monotonically along the straight line to Adam's end point" — tell the reader why that proves a kink and not a barrier: "…so no barrier separates the two points; the line search fails because the mean-absolute term of Eq. (4) has a non-differentiable kink wherever a residual sample changes sign, and at the stop the search direction crosses one." Also Fig. 4(d) shows L-BFGS leaving Adam's path at a steep angle; one clause saying the stop is on a different parameter combination (mostly PC2) than Adam's valley would reconcile text and figure.

**14. [MINOR] Sec. 8 opens with a 96-word sentence containing three new constructs** (natural cubic interpolants, θ^init + sϑ, hinge penalty). Split into three sentences:
"Every diagnostic so far points to the missing frequency dependence of the optical constants. We add it in an interpretable, well-scaled form: n_j(f) and κ_j(f) of each layer are natural cubic splines through their values at six equally spaced knots across the band in which the sample spectrum stays above the noise floor, 0.2–3.0 THz, and are held constant outside it. Each knot value is parametrised as θ^init + s ϑ with a trainable ϑ and a fixed scale s (the scaling of Section 3), and κ ≥ 0 is enforced by a hinge penalty λ_κ⟨max(0, −κ_j(f))⟩ with λ_κ = 1."

**15. [MINOR] Sec. 9 "here 2 %"** — say which fits: "…the spread between fits that weight the model error differently (Adam, L-BFGS on Eq. (4), L-BFGS on the rMSE: Fig. 4 caption), here 2 % in d."
Sec. 9 also never mentions the dispersive model or the air gap, which the abstract presents as the resolution. Add: "Letting the index disperse removes the artefact and leaves a residual carrying only the echo periodicity; the remaining 0.7 % offset of n_Si is the n·d valley again, not noise."

---

## B. Notation and language

**16. [MAJOR] σ has three meanings.** (i) Noise floor (Sec. 3, "2.7σ", "15σ", "±σ" bands, "(0.98–1.00)σ"); (ii) Cramér–Rao standard deviation ("Cramér–Rao 1σ" rows, "77σ from zero"); (iii) the text says "singular values" so σ_i is avoided there, but a reader of Eq. (6) will expect σ_i. Fix: write "77 times its Cramér–Rao bound" in Sec. 6.2(iv) and label the table rows "Cramér–Rao bound (1 s.d.)"; keep σ for the noise floor only.

**17. [MINOR] κ₂ vs κ_Si, n₁/d₁ vs n_q/d_q, θ vs s vs p vs x.**
- Table 2 and Sec. 6 use subscripts 1, 2; the same quantities are κ_Si in Sec. 6.2, 7 and 9; Table 3 uses d_q, d_air, d_Si. Pick layer subscripts (q, Si, air) throughout; they cost nothing and help the outsider.
- Sec. 3 introduces p_i (physical), s_i (scale factor), x_i (centred scaled), then Eqs. (4)–(5) use θ with no link. Add "θ = (s_i) is the vector the optimiser sees". Sec. 8 then uses ϑ and θ^init for yet another parametrisation; say it is the same scaling.
- Eq. (6): J_{ti} = ∂ŷ(t_t)/∂x_i — the double t is confusing. Use J_{ki} = ∂ŷ(t_k)/∂x_i.
- n_f in the Fig. 1 top-axis label and caption is undefined (film index).
- "N layers" (Sec. 3) next to "n" indices — fine, but say "N layers (N = 2 here)".
- "V-shaped" (Sec. 3) vs "cone-shaped" (Sec. 6.1) for the same loss: pick one.
- Power SNR "44 dB" (Fig. 1) vs "peak-to-σ ratio 7.9×10³" (Sec. 5): one convention.

**18. [LANGUAGE] Sentences over 45 words.** Twenty prose sentences exceed 45 words (abstract 74; Sec. 3 62; Sec. 4 68, 64; Sec. 6.1 57, 59, 70; Sec. 6.2 68; Sec. 7 66; Sec. 8 96, 54, 51). The worst are addressed in items 6, 8, 14. Three more:
- Sec. 3: "Both are local and the basin of a THz fit is narrow, so each fit is run from several starts around the nominal guess and the lowest loss kept; Bayesian optimisation of the start, used in earlier work, was abandoned because in fifteen dimensions it moved the start into the wrong echo-assignment basin at a cost of more evaluations than it saved." → split at the semicolon into two sentences, and define "echo-assignment basin": "…moved the start into a basin in which the model assigns the echoes to the wrong interfaces, at a cost of more evaluations than it saved."
- Sec. 6.1 last sentence: split at "; in Phase II".
- Sec. 7: "The two end points differ by 5.7 µm … how the loss weights the model error." → split at ": what L-BFGS reaches".

**19. [LANGUAGE] Colloquialisms that an Optica copy-editor will query.**
- "the simulated failure mode in the flesh" → "the failure mode of the simulation, observed experimentally".
- "Letting the model be right: …" (Sec. 8 title) → "A dispersive model: quartz–air–silicon with six frequency knots".
- "the negative κ₂ buys almost nothing" → "the negative κ₂ lowers the loss by almost nothing".
- "the honest uncertainty" → "the realistic uncertainty".
- "creep/creeping" ×4 (abstract, Fig. 3 caption, Sec. 7, Sec. 9) — fine once or twice; vary ("drifts").
- "at Adam's full speed" → "at the full step size η of Adam".

**20. [LANGUAGE] Inverted or ambiguous constructions.**
- Sec. 6.2 opening: "Whether the negative silicon extinction coefficient is a statistical fluctuation or a systematic effect the sensitivity matrix decides, step by step." → "The sensitivity matrix decides, step by step, whether the negative silicon extinction coefficient is a statistical fluctuation or a systematic effect."
- Sec. 3: "which tolerates the few large residuals near the main pulse … and is V-shaped around the optimum" — "which" refers to the loss, two lines above. → "This loss tolerates …".
- Sec. 3: "(of the measured reference for the simulated sweep)" → "; for the simulated sweep, which has no sample trace, the reference trace is used."
- Sec. 5: "converge to the same point to four figures" → "to four significant figures".
- Sec. 5: "2425 points at 33 fs" → "2425 samples at 33 fs spacing".
- Sec. 7: "history 20" → "memory of 20 gradient pairs"; "steps of 10⁻⁸" → "steps of 10⁻⁸ in the scaled parameters".
- Sec. 6.1: "whose Fresnel coefficients depend on n alone" — the echoes' Fresnel coefficients; say "whose amplitudes, through the Fresnel coefficients, depend on n alone".
- Sec. 6.1: "Another order of magnitude down are the two extinction coefficients" — 0.073 → 0.019 is a factor 4, 0.073 → 0.0016 a factor 45; say "one to two orders of magnitude".
- Tense: predominantly present, with "was measured", "We therefore repeated" in past. Acceptable; keep present for the analysis and past for the measurement, but make Sec. 7 "We therefore repeat the fit".
- British spelling is consistent (optimiser, analysed, artefact, normalised). No US forms found.

**21. [LANGUAGE] Affiliations 4 and 5 are listed but no author carries them** (author block is marked to-be-confirmed; noting it so it is not forgotten).

---

## C. Figures and captions

**22. [MAJOR] Figs. 3 and 4 (PCA planes): the caption does not say what the grey shading or the axis labels mean.** A reader must guess that lighter = lower loss, that "PC1 (78 %)" is the variance fraction, that "−0.69n₁ + 0.60d₂" are loadings in scaled coordinates, and that "+" is the final Adam point of Table 2. Add to the Fig. 3 caption: "Grey scale: loss on a log scale, lighter is lower. Axis labels: fraction of the path variance carried by each principal component and its composition in the scaled parameters x_i; (a)–(c) axes are additionally scaled by √G_ii so that a circle means equal curvature. + is the Adam end point of Table 2." Fig. 4 can then say "as in Fig. 3".

**23. [MINOR] Figs. 3 and 4 largely duplicate each other.** Fig. 4 already contains Adam's Phase-I and Phase-II paths (faint blue) and the same five panels. One figure with both optimisers (Adam solid, L-BFGS with markers) would convey the two-phase story and the comparison together and would free roughly a third of a page — enough to restore the five-layer figure (R14) or to enlarge Figs. 5–6. If the director keeps both, make Fig. 3's role distinct (e.g. drop (e) from Fig. 3, since Fig. 4(e) supersedes it).

**24. [MINOR] Fig. 4(e): three shades of blue and of orange are not explained** (one per start). Add "(one shade per start)" to the caption.

**25. [MINOR] Fig. 2(b) (quartz/Si residual): the ±σ band in the legend is invisible on the plot** (σ = 2.9×10⁻⁴ against an axis of ±0.1). Either drop the band from the legend, or state "±σ band is thinner than the line" in the caption, or add an inset at ×100.

**26. [MINOR] Fig. 2(a) inset "echoes" has no axis labels** (time in ps, same units as the main axis); add "(ps)" or state in caption. Fig. 2(c) and Fig. 5(c) y-axis "|spectrum|" has no units: "(arb. u.)".

**27. [MINOR] Fig. 6 caption refers to "top/bottom/left/right" while the panels carry (a)–(d).** Rewrite: "(a), (b) n(f) and (c), (d) κ(f) of quartz (left) and silicon (right)".

**28. [MINOR] Print size.** Figs. 5 and 6 are set at 0.62 and 0.57 of the line width; axis tick labels come out at ≈6 pt and the legend text in Fig. 5(b) is smaller still. The loadings in the axis labels of Figs. 3 and 4 (0.72 line width) are at the limit. Re-export Figs. 5–6 at narrower design widths (the writer already notes this).

**29. [MINOR] Fig. 4 caption end-point list is a positional data dump** ("999.6, 376.1 µm, +0.0044, 1.9528, 1.3146 mm" must be matched to d₁, d₂, κ₂, n₁d₁, n₂d₂ by counting). If the mini-table cannot be afforded, label each: "L-BFGS, same loss (250–299 evaluations): d₁ = 999.6, d₂ = 376.1 µm, κ_Si = +0.0044, n_q d_q = 1.9528 mm, n_Si d_Si = 1.3146 mm".

**30. [MINOR] Fig. 1 caption: "median (symbols) and range (bands)"** — say "over three noise draws" directly after "range", and define n_f. The top axis "film round trip 2n_f d/c" is the single most useful quantity in the figure; say in the caption that the 0.24 ps threshold of the text corresponds to the gap between the 10 and 7 µm points.

**31. [MINOR] Fig. 3(d) and Fig. 4(d): x-axis "[fraction of nominal]" appears only on panel (d)**, whereas (a)–(c) are curvature-normalised with no unit shown. Put "(curvature-normalised)" on the (a)–(c) axes or in the caption so the different units are visible.

---

## D. Tables

**32. [MINOR] Table 1.** Column "rMSE/σ" undefined (item 4). Row "independent" mixes a visible-band index (n_∥ at 633 nm), a profilometer thickness and a literature silicon index; add a footnote marker per entry ("a: white-light interferometry at 633 nm; b: profilometer; c: Ref. [14]") so the row is self-explanatory. κ_PI "<10⁻⁶" is fine. Significant figures: d_PI 25.84 (4 s.f.) vs independent 24.7 (3 s.f.) — acceptable since they come from different instruments, but say so in the footnote.

**33. [MINOR] Table 2.** "loss" column: say "loss, Eq. (4)" in the header. The row label "fit 1, κ ≥ 0 (softplus)" uses an undefined term (item 4). Bounds row: 6×10⁻⁵ next to 0.032 (µm) — fine, but the column header already carries the unit, so write 0.032 not "0.032 µm" in prose.

**34. [MINOR] Table 3.** Three values in one cell ("976.6, 125.9, 388.9") is readable but the header must say which is which — it does. σ = 9.0×10⁻⁴ here against 2.9×10⁻⁴ in Sec. 5 for the same sample re-mounted; one clause in the caption ("shorter acquisition") would stop the reader asking.

**35. [MINOR] Tables 1–3 use mixed formats for the same quantity type**: "1.35e-3" style (siunitx \num) in captions, "4.4522×10⁻³" in Table 2, "2.93e-3" in Table 3 body. Unify to ×10⁻³.

---

## E. Numbers in prose that carry no argument (brief item 3)

**36. [MINOR] Sec. 4, para 1 (five-layer stack):** fifteen parameter values and three error figures in a single sentence with no figure or table. Either restore the half-width five-layer figure and move the parameters into its caption (R14), or cut the parameters: "A noise-free five-layer stack with fifteen free parameters (values in notebook nb1), fitted by multi-start L-BFGS from a hand guess within n ± 0.3, κ ± 0.005, d ± 30 %, is recovered to 10⁻³ % in n and d and 0.3 % in κ."

**37. [MINOR] Sec. 5 instrument details** (1000 waveforms, 15 min, 2425 points, 33 fs, peak-to-σ 7.9×10³) belong in the Fig. 2 caption or a Methods-style sentence; the argument needs only σ.

**38. [MINOR] Sec. 6.1:** the six singular values, ρ = 0.97/0.99, 280/890, "95 %", 1615, 81–92 %, 92 %, 22–29 % — the phase boundary (1615) is already in the Fig. 3 caption; the rest can move there or to the proposed mini-table (item 6). Keep 5500 and "an order of magnitude" in prose.

**39. [MINOR] Sec. 6.2:** "31 nm", "1.9×10⁻⁵", "β = 5", "0.6 %", "0.015 %" are in or belong in Table 2 / its caption. Keep 77× and 0.29 % (they are the argument).

**40. [MINOR] Sec. 7:** "5.7 µm", "0.012 %", "0.034 %" duplicate the Fig. 4 caption. Keep 250–300 vs 5000–5900 and 3.7 %.

**41. [MINOR] Sec. 8:** "126 µm, 5.5 times, 50 µm", "3.4175", "0.01 %", "3.37", "+4 %" could go to the Table 3 caption; keep "0.7 %" in prose because it carries the point. λ_κ = 1 and the knot positions are in Table 3 already.

---

## F. Structure and length

**42. [MINOR] Merge candidates for space:** (i) Figs. 3 and 4 (item 23); (ii) Sec. 6.2's last three lines (softplus map) are an implementation note that could be the Table 2 caption; (iii) the single sentence on the combined time+frequency loss (R9) could be the last clause of the preceding sentence. Space gained would be best spent on the sensitivity mini-table (item 6), the five-layer figure (R14) and larger Figs. 5–6.

**43. [MINOR] Missing for a reader outside the field:** (i) one sentence on what a THz-TDS measurement is (a reference pulse through air and a sample pulse through the stack, both time-resolved fields) — the abstract says "one transmission trace" and Sec. 2 assumes x(t) is known; (ii) the sample geometry in Sec. 5 is a to-do; (iii) the Introduction has no explicit statement of the three contributions (limit, landscape, optimiser comparison) — one sentence would help the abstract and conclusions line up.

**44. [LANGUAGE] Heading "6.2 Noise-limited uncertainties, and why κ_Si < 0 is a model artefact"** is fine; but the subsection carries two labels (sec:cr, sec:constrained) and Sec. 6.1 two labels (sec:sens, sec:pca) — harmless, but the writer should check that no reference intends a subsection that no longer exists.

---

## Counts
[MAJOR] 10 (items 1–8, 16, 22) · [MINOR] 28 · [LANGUAGE] 6 · total 44.

## Disagreements with rulings
- R9/brief item 9 asks for a *new* L-BFGS trajectory figure. I do not disagree with showing the L-BFGS path, but Fig. 4 already contains Fig. 3; two near-identical five-panel figures cost a third of a page that the paper needs for the sensitivity mini-table and the five-layer figure (R14). A single combined figure still satisfies the brief.
- Brief item 3 (numbers in tables/captions): Sec. 6.1 currently violates it most. Item 6 proposes the fix; R15 already allows the sensitivity panel/table "if it fits".
- No other disagreements.
