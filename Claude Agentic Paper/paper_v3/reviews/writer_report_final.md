# Writer report — final pass (loop 3)

**`main.tex`: 8 pages** including the references. `bash briefs/compile.sh` gives 0 errors and 0 undefined references, and `main.log` has 0 overfull boxes. Table 3 is set on p. 7, before the reference list on p. 8. The abstract is 200 words.

Backups of the pre-pass files are in `scratch/main_pre_final.tex` and `scratch/refs_pre_final.bib`.

## Numbers typed in this pass
Each was checked against the JSONs with a python one-liner:
- **log-cosh share 1.6 %:** (4.4516e-3 − 4.3798e-3)/4.4516e-3 = 1.61 % (`nb3 best`).
- **26 %:** max std/CR = 1.259 (`nb2 C11_sweep/b_ten_draws`).
- **Adam and L-BFGS counts:**
  - 120–150 L-BFGS evaluations (123/152/154);
  - 360–460 iterations at η = 3e-3 (357/437/460);
  - 620–660 iterations to reach 4.5e-3 (621/662/622);
  - 4.5e-3 is below the stall level of 4.617e-3 (`nb4 loop1_checks/C7_adam_learning_rate`).
- **0.06 %:** the 16-end-point invariants spread 0.026 / 0.031 / 0.057 % (`nb4 end_points_timing_invariants`).
- **Group n·d 1928.2:** `nb5 loop2_checks/C16_run3_optical_thickness nd_group_1THz_um` = 1928.18. From it:
  - 1952.54 vs 1928.18 gives 1.25 %, printed as 1.2 %;
  - 1931.36 vs 1928.18 gives 0.17 %, printed as 0.2 %.
- **2.2 µm / ~7 fs:** 1367.70 − 1365.55 = 2.15 µm (`nb6 C16` gap row and `nb5 C16_run3`).
- **σ clause:**
  - 3000 waveforms over 115 ps;
  - σ ratio 9.04e-4 / 2.90e-4 = 3.1, printed as "three times";
  - white levels 2.69e-5 vs 4.51e-5 (`nb5 loop2_checks/sigma_note`).
- **Table 2, knots + gap row:** switched to run 3, the lowest rMSE (971.83 / 398.31; κ_Si −0.0126 to 0.0047; `nb6 C3_C4_six_knots_plus_gap`). The text's "+0.005" still holds.

## MUST items

| item | status | location |
|---|---|---|
| Domain M1 / Outside C1: σ clause | done | Sec 8 para 1, last sentence: 3000 waveforms, 115 ps; three times Sec 5 because two noise bursts fall before the pulse; white level 2.7e-5 vs 4.5e-5; multiples of σ do not compare. Fig. 5 caption: "raised by pre-pulse noise bursts; white level below that of Fig. 2". |
| M2: η = 3e-3 sentence | done | Sec 7 para 1: "reaches it in 5000–5900 iterations at η=1e-4 and 360–460 at 3e-3 (none diverging), and 4.5e-3, below L-BFGS's best, in 620–660" |
| M3: 0.01 % precision | done | Sec 9: "1 %, against a 0.06 % reproducibility" |
| M4 / C2: group n·d | done | Table 3 caption: "1926.0 (q; group value 1928.2, which the echo timing measures)". Sec 9: "Quartz group n_q d_q …"; the 7 fs is traced to "2.2 µm in d_q + d_Si", which follows from Tables 2 and 3. The extra Sec 8 "Fourth" sentence (C2b) was not added, for space. |
| M5: log-cosh share | done | Sec 3: "adds 1.6 % to the loss at the Adam end point" |
| M6: 25 % → 26 % | done | Sec 4 para 1 |
| M7 / C6: thin-film conclusion | done | Sec 9: "A simulated thin film stays at the CR bound; its error in d passes 1 % where the round trip falls below two thirds of the pulse width." |
| C3: Table 1 footnote | done | Footnote now reads "$^d$product of the values marked $^a$ and $^b$, a guide only" |
| C4: Sawanaka | done | refs.bib `Sawanaka, Ken-ichi`; prints "K.-i. Sawanaka" |
| C5: Δ, ϑ undefined | done | Sec 8: "a trainable variable ϑ times a fixed scale Δ" |
| R38: τ criterion | done | Sec 6.2 (iv): "we fit τ only where, as for PI, it brings a fixed-index substrate into agreement and the acquisition motivates it". The manuscript nowhere states the "< 1 %" misfit change and no longer says "cuts the misfit strongly". |

## SHOULD items

**Done:**
- S1: "120–150 evaluations".
- S2: in part. The criterion follows R38, but the clause "a free τ also lowers the misfit here" was cut for space.
- S3: "for five-layer stacks".
- S4: knots + gap row uses run 3.
- S5: Sec 8 opens with "To see what a dispersive model adds".
- S8 / C10a: s₁ ≥ … ≥ s₆.
- S9: same as C4.
- S10 / C13 (Table 2): "the model of Section 8 on these data … ranges over the Si knots" and "rMSE change relative to row 4".
- C7: abstract now reads "imposes the physical sign, κ ≥ 0".
- C9: "0.16 µm in d_q and d_Si".
- C10b: "remaining soft n–d direction".
- C11a: "do not remove the negative sign".
- C11b: "at a cost of 0.3 % in loss".
- C12: the Sec 3 sentence is split.
- C14: Fig. 1 caption ("Symbols: median error; bands: range"; "Hollow symbols (offset for clarity)").
- C15b: "A frequency-domain loss term adds no information (same CR bounds, white noise, full band)".
- C16: Fig. 6 caption has a title, and HR-Si is spelled out and cited.
- C17: "on this stack the positive sign comes from the modelled air gap, not dispersion" (shortened to "here …").
- C18:
  - the loss profile is defined as the loss minimised at fixed n_Si;
  - the stage-1 wording is "each parameter's knots tied (stage 1, non-dispersive)".
- C19: "the fits with and without τ".
- C20a: Sec 9 reads "n_q d_q". C20b: the Table 3 caption reads "n_j d_j". C20d: the THz-TDS acronym is dropped.
- C21: the Fig. 3 caption reads "loss falls by 81–92 % / 22–29 %".
- C22a: `publisher = {Curran Associates}` added to Paszke2019 and Li2018. The Li2018 pages were not added because they are unconfirmed.
- C22b: Kingma is now `@inproceedings` with a single year.
- Bayesian-seeding sentence: reworded.

**Not done:**
- S6: "(except two rows of Table 2)", cut for space.
- S7: Fig. 4(d) legend. This is the coder's job (nb4 re-export).
- S11: "7 fs" into nb6. This is the coder's job and optional.
- C8: PI time offset in the conclusions, for space.
- C11a: the full 57-word rewrite was not done; only the antecedent was fixed.
- C13 (Table 3 column explanations): partly covered by the text ("stage 1 … tied"); the caption addition was cut for space.
- C15a: the kink-sentence split, for space.
- C20c: sample-name harmonisation beyond the above.

## Layout changes to stay at 8 pages (please check)
- Fig. 3 (fig4_adam_two_phase) and Fig. 4 (fig5_lbfgs_trajectory) are now 0.95\linewidth, and Fig. 6 (fig7_qas_nk) is 0.9\linewidth.
- Short trailing lines were removed by light rewording:
  - Sec 1 (i)–(iv);
  - Sec 2: "field continuity gives"; "exactly" dropped;
  - softplus paragraph: "unlike a rectifier, never freezes";
  - Sec 7: "L-BFGS's advantage is needing no step size".
- p. 7 is now a float page (Table 3, Figs 5 and 6).
