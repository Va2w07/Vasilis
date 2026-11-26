# Domain review — loop 3 (final pass)

Reviewer: domain reviewer (THz-TDS extraction, inverse problems, optimisation).

**Material read:** DECISIONS_LOG.md (R1–R37), reviews/loop2_domain.md, writer_report_loop2.md, coder_report_loop1.md (tail), nb4/nb5/nb6 `loop2_checks` (C15, C16, C17, sigma_note), main.tex, reviews/main_loop3.pdf (8 pp., read page by page as images), figs/*.png as placed in the PDF, results/nb1–nb6.json.

**Conventions:**
- Line numbers refer to main.tex.
- Figure numbers follow the PDF (Fig. 1 = fig2_thin_film … Fig. 6 = fig7_qas_nk).
- **[MUST]** marks an error: a wrong number, a wrong or unsupported statement, or a ruling not implemented.
- **[SHOULD]** marks an improvement. No new analyses are requested.

**Verdict: accept with the corrections listed below.**
- The seven [MUST] items are each a word-to-sentence edit.
- None changes a conclusion.
- No new runs are needed: every corrected value is already in the JSONs.

---

## 1. Status of the loop-2 items

| loop-2 item | topic | status | where / remark |
|---|---|---|---|
| **3.1 MAJOR** | n·d accuracy: interface model, not time zero | **resolved** | Abstract l.47 ("reproducibly … within one interface model; an interface gap shifts it by 1 %"); Sec 6.1 l.141 (scope "interfaces without intermediate layers"; "within the two-layer model"; "interface layer moves n_q d_q by 1 %"); Sec 9 l.205. **Residue:** two numbers in the new Sec 9 sentence are wrong or unreproducible (M3, M4). |
| 2.1 | end-point spread wording | resolved | l.112; Table 2 caption "lowest-loss iterate (start 3)" |
| 2.2 | "4–10" → 3.6–9.6 | resolved | l.154 |
| 2.3 | Fig 4(e) y-label | resolved | PDF: "loss" (also Fig 3(e)) |
| 2.4 | "offset in κ_q" | resolved | Fig 4 caption: "0.05 off it (in n, d), 1.13 along it" |
| 2.5 | "about thirty iterations" | resolved | l.160: 8–39 iterations (22–96 evaluations) = nb4 C17 |
| 2.6 | SNR 44 dB → peak/σ | resolved | Fig 1 caption: 2.1×10³ |
| 2.7 | "order of magnitude" → "eight times" | resolved | abstract, Sec 8, Sec 9 |
| 2.8 | "mean Si knot 0.8 % low" | resolved | l.201 |
| 3.2 | "imposes"; n_Si +1.5 %; "weakly determined" | resolved | abstract; l.154 (iii); l.156 |
| 3.3 | "richer models do not remove it" | resolved | l.154 (iv) |
| 3.4 | criterion for free τ | resolved, with a wording gap | l.154. The "cuts the misfit strongly" half does not discriminate (S2). |
| 3.5 | Adam at η = 3×10⁻³ below the L-BFGS stall | **partly** | l.160 is garbled (M2). |
| 3.6 | Lewis & Overton: strong vs weak Wolfe | resolved | l.162 |
| 3.7 | ripple demoted | resolved | l.201 ("possibly", "not significant at residual-scaled bounds"); C15 ratios 0.90 / 0.58 |
| 3.8 | n_Si band is a common offset | resolved | Fig 6 caption |
| 3.9 | hinge "penalises"; κ_Si vs HR-Si | resolved (first part) | l.172. The κ_Si-vs-HR-Si clause was not added. No action: the text makes no claim about κ_Si values. |
| 3.10 | air-path factor in Eq. (1) | resolved | l.67–69 |
| 3.11 | protocol / data availability | resolved | l.78 (box, starts, hinge weight 1); Table 2 caption (starts); l.112 (budget); l.214 |
| 3.12 | PI: "two best starts"; profilometer uncertainty | resolved | l.93; the uncertainty is an author `\todo` |
| 3.13 | s₃ label | resolved | Fig 3 caption |
| 3.14 | round-trip criterion specific to the film | resolved | l.91 "for this film"; Fig 1 caption "top axis: simulated film only", 0.29 ps |
| 3.15 | gap "of known index" | resolved | l.201 |
| C15 / C16 / C17 | checks | done | nb5, nb6, nb4 `loop2_checks` |

Summary: 22 of the 23 minors are resolved; 3.5 is partly resolved (M2). The MAJOR is resolved in substance, with two number slips.

**R37 σ clause in Sec 8: not implemented** (M1). Neither main.tex nor main_loop3.pdf contains it.

---

## 2. Number audit

**Method.**
- **Script** (`/tmp/a3/audit.py`, reused from loop 2): it flattens nb1–nb6.json (7836 leaves) and matches every number of the abstract-to-backmatter body at its printed precision. It reports 255 tokens, two of them unmatched: 5900 (a rounded range end, 5929) and 620 (see M2).
- **Hand check:** every number in the prose, all captions and Tables 1–3, against the named JSON key. I also diffed the loop-2 and loop-3 PDF texts to isolate new numbers.

**Verified correct** (unchanged from loop 2 unless marked new):
- **Fig 1 caption:** 1.58e-3, 2.1e3 (2107), 13 / 3 draws. PI-like film: n 1.72, κ 0.01, 5 draws, 0.29 ps.
- **Sec 4:** 7–10 µm, 0.36 ps, two thirds. PI: 3.26, 23 min, 10 fs ↔ 3 µm, misfit −37 %, n_Si +0.36 %, d_PI +2.2 % (free and fixed n_Si), 0.2 % (0.217), 0.16 %, 10–15 µm.
- **Table 1:** all entries, including the CR row (`nb2 loop1_checks/fig2b_variant`).
- **Sec 5:** 2.3e-4 in n_q; 0.16 µm in d (d_Si 0.157, d_q 0.153); <0.02 %; +0.9 %; 41σ.
- **Fig 2 caption:** 2.9e-4, 7.9e3, 15σ, 41σ, 2425 × 33 fs (80 ps, matching sigma_note).
- **Table 2:**
  - all rows and the CR row;
  - caption rMSE changes −5.0 / −4.3 / −57 / −59 % (−58.6);
  - starts;
  - "41-fold".
  - The knots + gap row is run 1 of C3_C4 (see S4).
- **Sec 6.1:** 0.06 % (invariants 0.026 / 0.031 / 0.057 %), 2 % (2.25), 1 % (1952.3 → 1931.4 = −1.07 %).
- **Fig 3 caption:** 1–1615, −81–92 %, 99.7 %, −22–29 %, s₁–s₆.
- **Sec 6.2:** 77, ≈100 (108), <2, 3.6–9.6, −0.007…+0.004, −0.013 / +0.005, 3.64, 0.3 %, +1.5 % (+1.48).
- **Sec 7:** 8–39 / 22–96 (C17); 300 (299); 5000–5900 (4998–5929); 360–460 (357–460); 3.7 %; +0.023 / −0.036.
- **Fig 4 caption:** 12–55 %, 97–99.8 %, 0.05, 1.13, six decades, 0.2 nm.
- **Table 3:** all 24 knots (`nb5 stage2/knots`), the stage-1 column, d 976.6 / 125.9 / 388.9 and 974.5 / 126.8 / 389.1, CR ranges, 1926.0 / 1319.3, −0.98.
- **Fig 5 caption:** 9.0e-4, 3.62 → 2.93e-3, −19 %, 4.0σ → 3.2σ.
- **Fig 6 caption:** 6.4, "two bounds".
- **Sec 8:** eight times, 0.7 %, 1.19 THz, 0.8 %, 0.05 %, 3.38–3.43.
- **Sec 9:**
  - 2 % / 6 % (2.25 / 5.83);
  - 0.2 % (1931.4 vs group 1928.2 = 0.165 %);
  - ~7 fs (d_q + d_Si: 1367.70 − 1365.55 = 2.15 µm = 7.2 fs; derived from C16 values).

**Discrepancies** (all listed as corrections in §5):

| where | printed | correct (JSON) |
|---|---|---|
| l.83 log-cosh share | "about 1 %" | 1.6 % of the loss at the Adam row (loss 4.4516e-3 − MAE 4.3798e-3 = 7.2e-5; `nb3 best`) → M5 |
| l.91 scatter vs CR | "within 25 %" | std/CR 0.77–1.26 (`nb2 C11_sweep/b_ten_draws ratio_std_err_to_cr`, n at index 3 = 1.259) → "within 26 %" (M6) |
| l.160 "passing below it in about 620" | — | Adam *reaches* the L-BFGS level in 357–460 (`iterations_to_lbfgs_stall`); 621–662 is the count to 4.5×10⁻³ (`iterations_to_4p5em3`), 2.5 % below the stall → M2 |
| l.160 "about 150 evaluations" | about 150 | 123–154 (`lbfgs_evaluations_to_stall`) → S1 |
| l.205 "0.01 % precision" | 0.01 % | Nothing in the JSONs gives 0.01 %. The 16 end points spread 0.031 % in n_q d_q and 0.057 % in n_Si d_Si (`nb4 end_points_timing_invariants`), and Sec 6.1 prints 0.06 %; the CR bound of n_q d_q is 8 nm = 0.0004 % → M3 |
| l.205 "1.2 % … 0.2 %" | consistent only with the *group* n d (1928.2, `nb5 loop2_checks/C16_run3`) | Table 3 prints only the *phase* n d (1926.0). With it a reader gets 1.4 % and 0.3 % → M4 |
| Table 2, knots + gap row | 971.79 / 398.34 | Run 1. C16 and the lowest rMSE (by 4×10⁻⁷ relative) use run 3: 971.83 / 398.31 → S4 (cosmetic) |

---

## 3. Referee reading: abstract, Sec 6.2, Sec 7, Sec 8, Sec 9

- **Abstract.** Correct and within the data. Every claim has a body counterpart: noise-limited, τ for PI, 0.1 % within one interface model, 1 % gap shift, two soft directions, Adam drift, estimator-dependent sign, "imposes", the L-BFGS kink, eight times. No correction.
- **Sec 6.2 (κ_Si argument).** Steps (i)–(iv) are now correct and well ordered. One logical gap, in the τ criterion (S2):
  - The text adopts τ "only where, as for PI, it cuts the misfit strongly and restores a fixed-index substrate".
  - On quartz/Si a free τ also cuts the misfit strongly: rMSE −30.4 % (`nb6 C4 summary/tau_amp`), against −37 % for PI.
  - So only the second condition discriminates. The statement is not false, but a referee will see that the first condition is idle.
  - **Log note for the director:** R35's premise "on quartz/Si it changes the misfit by < 1 %" is contradicted by nb6 (−30 %). The text does not use it.
- **Sec 7.** One error (M2). Otherwise correct: evaluation counting, kink evidence, weak/strong Wolfe, rMSE convergence.
- **Sec 8.**
  - The σ clause is missing (M1).
  - The opening question ("whether the misfit of Section 5 is missing dispersion") is not answered in Sec 8. The answer is in Table 2: six knots on the Sec 5 data cut the rMSE by 57 % but leave κ_Si < 0 (S5).
  - Without the σ clause a reader compares "MAE 15σ" (Sec 5) with "3.2σ" (Sec 8) and infers a five-fold better fit. In absolute terms the MAE improves only 1.5-fold (4.38e-3 → 2.92e-3), because σ is three times larger in run 3. This is why M1 matters.
- **Sec 9 / conclusions.**
  - M3 and M4 (numbers).
  - M7: "Thin films reach the CR bound down to a round trip of two thirds of the pulse width" mis-states Sec 4. The sweep sits at the CR bound at *every* thickness, including 5 µm (std/CR 0.77–0.80). Two thirds of the FWHM is where the error in d passes 1 %, and only for the simulated high-contrast film (the PI-like film crosses at one third to one half; D1).
  - The rest is correct.

---

## 4. Figures against captions (final)

All six figures match their captions: panels, legends, markers, "noise-limited 1σ", "Adam end point", "loss", the Phase-I/II split at 1615, the 3.7 % inset, hollow edge knots, the dashed HR-Si line, stage 1 dotted. Fonts are legible at print size. One cosmetic point:
- **S7:** in Fig 4(d) the legend entries "start" and "end" (grey markers) sit on the darkest contour region at top left and are barely readable. Move the legend or give it an opaque white box (coder, nb4 re-export only).

---

## 5. Corrections

### [MUST]

**M1. Sec 8: the σ clause required by R37 is missing.**
- **Check against nb5.json → loop2_checks/sigma_note:**
  - run 3: 3000 waveforms over 115 ps, against 1000 over 80 ps in Sec 5;
  - pre-pulse std 9.0e-4 against 2.9e-4;
  - two 2-ps bursts: 1.1e-3 at 4–6 ps and 1.7e-3 at 8–10 ps;
  - other segments 6.9e-5 to 1.7e-4;
  - white (first-difference) level 2.7e-5 against 4.5e-5.
- All of these agree with the `segment_std` and `white_min` arrays.
- **Proposed text** (Fig 5 caption or after "(stage 2)" in l.172): "Its σ, 9.0×10⁻⁴ (3000 waveforms, 115 ps), is three times that of Section 5 because two bursts fall in the 10-ps pre-pulse window; its white level is lower (2.7 against 4.5×10⁻⁵), so multiples of σ are not comparable between Sections 5 and 8."

**M2. Sec 7 l.160: "360–460 at 3×10⁻³ (none diverging), passing below it in about 620" is wrong.**
- Adam is already *at* the L-BFGS level after 357–460 iterations.
- 621–662 is the count to 4.5×10⁻³, which is 2.5 % below the stall level.
- Write: "…360–460 at 3×10⁻³ (none diverging), and reaches 4.5×10⁻³, below anything L-BFGS attains on this loss, in 620–660".

**M3. Sec 9 l.205: "1 % against a 0.01 % precision".**
- 0.01 % is in no JSON.
- The paper's own precision figure is 0.06 % (Sec 6.1: 16 end points; quartz alone 0.03 %).
- Write "against a 0.06 % reproducibility". Alternatively "0.03 %" if only quartz n·d is meant.

**M4. Sec 9 l.205 and Table 3 caption: the 1.2 % / 0.2 % comparison uses the run-3 *group* optical thickness (1928.2 µm), but only the *phase* value (1926.0 µm) is printed.**
- From the printed values a reader obtains 1.4 % and 0.3 %.
- The group value is the right comparison, because a frequency-independent fit times the echoes.
- Add it to the Table 3 caption: "$n\,d$ at 1 THz 1926.0 (q; group 1928.2), 1319.3 µm (Si)". In Sec 9 write "Quartz group $n_jd_j$ differs by 1.2 %…".

**M5. Sec 3 l.83: "the log-cosh term … adds about 1 % to the loss".** At the Adam row it is 1.6 % (7.2e-5 of 4.452e-3; `nb3 best`). Write "about 2 %" or "1.6 %".

**M6. Sec 4 l.91: "equals the CR bound within 25 %".** The largest std/CR is 1.26 (`C11b`, n). Write "within 26 %", or "(std/CR 0.77–1.26)".

**M7. Sec 9 l.205: "Thin films reach the CR bound down to a round trip of two thirds of the pulse width".** This mis-states Sec 4: the sweep is at the CR bound at all thicknesses, and two thirds is the 1 % threshold for the simulated film only. Write: "A simulated thin film is recovered at the CR bound; its error in d passes 1 % where the round trip falls below two thirds of the pulse width".

### [SHOULD]

- **S1.** l.160 "about 150 evaluations" → "120–150" (123–154).
- **S2.** l.154 (iv), τ criterion: "(a free τ also lowers the misfit here, but drives n_Si to 3.64; we adopt τ only where, as for PI, it restores a fixed-index substrate)". Uses only numbers already printed.
- **S3.** l.87: "in fifteen dimensions" refers to the five-layer problem, which is no longer in the paper. Write "in higher dimensions (fifteen parameters)", or "for five-layer stacks".
- **S4.** Table 2, knots + gap row: use the lowest-rMSE run (971.83 / 398.31, as in nb6 C16), or note "start 1". The difference is 0.04 µm.
- **S5.** l.172: the opening "To test whether the misfit of Section 5 is missing dispersion" is not answered in Sec 8. Either rephrase ("To see what a dispersive model adds, we fit…") or add "(on the Section 5 data six knots cut the rMSE by 57 % but keep κ_Si < 0, Table 2)".
- **S6.** l.78: "Until Section 8 the indices are frequency-independent" → add "(except two rows of Table 2)".
- **S7.** Fig 4(d) legend legibility (§4).
- **S8.** The index k is used for both the time sample (Eq. 3, J_ki) and the singular values s_k (l.141). Write s₁…s₆ or s_m.
- **S9.** Ref. [6] prints "K. ichi Sawanaka". Brace the given name in refs.bib: `{Ken-ichi}`.
- **S10.** Table 2 caption "Six knots: Section 8" → "six-knot model of Section 8 on these data" (it is the Sec 5 trace, not run 3).
- **S11.** Sec 9 "~7 fs" is derived arithmetic from C16 (1367.70 − 1365.55 µm). For strict R10 traceability the coder can write it to nb6 `loop2_checks`. Optional.

---

## 6. Disagreements with rulings

None on the text.
- **Factual note on R35** ("on quartz/Si [τ] changes the misfit by < 1 %"): nb6 C4 gives −30.4 %. The manuscript does not repeat the error, but S2 makes the criterion honest.

---

## Counts

**[MUST] 7** (M1–M7) · **[SHOULD] 11** (S1–S11) · total 18.

Loop-2 status: MAJOR 3.1 resolved; 22 of 23 minors resolved; 1 partly resolved (3.5 → M2); R37 σ clause not implemented (M1).
