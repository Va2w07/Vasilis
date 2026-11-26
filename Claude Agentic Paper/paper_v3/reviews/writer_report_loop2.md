# Writer report — loop 2

**`main.tex`: 8 pages** including the references. `briefs/compile.sh` reports no errors and no undefined references, and the log has 0 overfull boxes. Table 3 is now set before the reference list, on p. 7.

New numbers come from these JSON keys:
- `nb6.json loop2_checks/C16`: gap 1931.4, adhesive 1929.1, τ-free 1952.2.
- `nb5.json loop2_checks/C16_run3`: group n·d 1928.2.
- `nb4.json loop2_checks/C17`: 8–39 iterations, 22–96 evaluations.
- `nb4 C7`: Adam at η = 3e-3 passes below the stall in 621–662 iterations.
- `nb2 partA/peak_over_sigma`: 2.1×10³.

`White1980` (sandwich covariance) is added to `refs.bib`, checked on Crossref (Econometrica 48(4), 817).

## How the space was found (please check)
The new conclusions, the abstract draft and all the requested definitions took the paper to 9 pages. These steps brought it back to 8:
- Tables 1–3 are set in `\scriptsize`.
- The float separations are 6 pt instead of 8 pt.
- The five-layer sentence is removed (R4/R14 allow it).
- Captions are cut hard: the Fig. 3 singular-value list no longer carries loadings, and the Table 2 caption is compressed.
- Prose is condensed throughout, to about 2000 words.
- The two six-knot rows are back in Table 2 (outside item 4).

## Domain MAJOR
- **3.1 (time zero vs interface model, R33):**
  - Sec. 9 now says that quartz n·d differs by 1.2 % between the two mountings, but by only 0.2 % once an interface gap is modelled. Time zero accounts for just ~7 fs. So the interface model sets the accuracy (~1 %) against a precision of 0.01 %.
  - The counting argument is scoped to "interfaces without intermediate layers". The 16-end-point statement now reads "within the two-layer model", and the text adds that an interface layer moves n_q d_q by 1 %.
  - In the abstract, "0.1 %" now reads "reproducibly … within one interface model; an interface gap shifts it by 1 %".

## Outside MAJOR
1. **Conclusions summarise the paper:** they now cover the thin-film limit, the model-set error budget, interface accuracy, Adam drift, L-BFGS kink, estimator spread, the constraint and the six-knot result. They are about 130 words rather than 190, for space.
2. **Table 3 before the references:** done.
3. **16 end points vs gap rows:** handled by the R33 wording above.
4. **κ_Si thread:**
   - negative κ is stated to mean gain;
   - the hinge sentence now reads "in a separate run without the hinge";
   - the six-knot rows are in Table 2.
5. **Symbol clashes:**
   - x is now the reference pulse only; the optimiser variable is u and the scaled coordinate is ξ_i;
   - the ambient index is n_a, start values are n_st, and the knot scale is Δ;
   - captions match the new legends ("noise-limited 1σ", "Adam end point", "loss").
6. **Step (iv) readability:**
   - the sandwich covariance is defined with its formula and cited [White 1980];
   - "least-squares point" is named as the L-BFGS rMSE row;
   - the ratio is given as "below 2 … and 3.6–9.6".
7. **Abstract:** based on the outside reviewer's draft, with Adam's drift, "eight times", and the time-zero and interface caveats.
8. **Two pipelines:** one sentence says why.
9. **Kink:** the corner of |r| is explained in plain words; Lewis & Overton is now cited for strong-Wolfe versus weak-Wolfe line searches.

## R35 minors done
- "imposes", with n_Si at +1.5 %
- "do not remove it", with the −0.013…+0.005 swing
- a criterion for when τ is fitted
- the η = 3e-3 result
- the ripple demoted to "knot-grid resolution, not significant at residual-scaled bounds"
- the n_Si band reconciled as a common offset, two residual-scaled bounds from the knot mean
- "penalises κ < 0"
- the e^{i n_a D} factor in Eq. (1)
- s₃ relabelled
- the Fig. 1 top axis described as "simulated film only", with the PI round trip of 0.29 ps
- protocol additions: starts and box, hinge weight 1, quartz/Si starts in the Table 2 caption, the 2×10⁴ budget in Sec. 5
- the number slips: 3.6–9.6, the end-point spread, peak/σ, C17 iterations, and "mean Si knot 0.8 %"
- PCA, FWHM and L-BFGS expanded at first use
- the title now reads "multilayer"
- the data-availability sentence now names the notebooks, with the repository DOI left as a `\todo`

## Not done, or reduced, and why
- **C15 numbers (knot correlation ≥ 0.98; alternation significance):** these are not in any JSON. The text states the conclusions without numbers (R10); the coder should add them if numbers are wanted.
- **Outside item 20 (σ 9.0e-4 vs 2.9e-4):** the cause is unknown to me. It needs an author note.
- **Outside items 22–25 (figure-internal fixes such as label overlap and the Fig. 5(b) unit):** these belong to the coder. "±σ band thinner than the line" was cut from the Fig. 5 caption for space.
- **"Quartz" vs fused silica naming (outside item 18):** left as an author `\todo`.
- **Conclusions length:** about 130 words rather than the outside reviewer's 190, for space.
- **"Held constant outside the band":** restored.
