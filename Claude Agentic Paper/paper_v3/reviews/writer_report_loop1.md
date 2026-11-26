# Writer report — loop 1

**`main.tex`: 8 pages** with the coder's real figures (`briefs/compile.sh`: no errors, no undefined references, no overfull boxes). Every `\XX`/`\nb{}` placeholder replaced with values read directly from `results/nb1–nb5.json` (macros kept in the preamble, unused). Author `\todo{}` items kept (author block, sample details ×2, sample material, mechanical thickness, funding, acknowledgments, disclosures, data link).

## Numbers now in the paper (all from the JSONs)
- nb1 (Sec 4, one sentence, no figure): max errors n 0.001 %, d 0.002 %, κ 0.3 %.
- nb2 A (Fig. 1 caption, Sec 4): σ = 1.58e-3, SNR 44 dB; medians at 5/10/20 µm: n 1.9/0.29/0.03 %, d 4.8/0.99/0.08 %, κ 25/8.3/6.3 %; threshold 0.24 ps ≈ 10 µm; converged loss (0.98–1.00) σ.
- nb2 B (Table 1, Sec 4): n_PI 1.712, d_PI 25.84 µm (+4.6 % vs 24.7, 18 CR-σ), n·d 44.26 µm, n_Si 3.259, d_Si 400.2 µm, rMSE 2.69 σ; CR bounds (d_PI 62 nm), ρ(n,d) = −0.93; local minimum of 2/8 starts (28.47 µm); collapsed full-band fit from 16 µm (22.3 µm, +48 %). Note: n_PI(THz) is 0.8 % *below* the visible n∥, so the text no longer claims the THz index exceeds the visible one; it presents the n–d split as marginal, as the sweep predicts. "Six of eight starts agree to within 1 %" (JSON spread 0.9 %, not the 0.1 % of the coder report).
- nb3 (Sec 5–6, Table 2): σ = 2.9e-4, peak/σ 7.9e3, MAE/σ 15, fit 1 and softplus rows, CR bounds (31 nm, κ₂ 1.9e-5), |κ₂|/bound 77, loss +0.29 % under the constraint, singular values 8.66…0.0016 (×5500), ρ 0.97/0.99, κ₁ 280× and κ₂ 890× less visible, phase boundary 1615, Phase I −81…92 %, Phase II 92 % of budget for −22…29 %, plane variance 99.7 %, valley axis 0.97κ₂+0.24κ₁.
- nb4 (Sec 7, Fig. 4 caption): L-BFGS 250–299 evaluations to the level Adam reaches at 5000–5900 iterations, stall 3.7 % above Adam at a kink (line between end points monotonic), end points 999.6/376.1 µm with κ_Si +0.0044 vs Adam 994.3/381.4, n·d within 0.012 %/0.034 %; rMSE and log-cosh converge to 977.7/398.1 µm (κ_Si −0.0073; start 1 of log-cosh still creeping at 987.7 µm after 908 evaluations).
- nb5 (Sec 8, Table 4, Figs 5–6): unchanged from loop 0 plus the d_air scan sentence (minimum 126 µm, 5.5× lower than at 50 µm) and the n_Si valley scan (floor near 3.37, literature value +4 % in loss).

## Rulings applied
R11 Sec 7 rewritten (kink, not barrier; smooth loss as companion; numbers in the Fig. 4 caption — the separate end-point table was tried and removed for space). R12 σ defined once in Sec 3. R14 five-layer figure tried at 0.5\linewidth: the paper then ran to 9 pages, so it went (R14's own condition); the five-layer result is one sentence with its parameters. R15 `fig_sensitivity` not used (numbers in Sec 6.1); air-gap scan in one sentence.

## Cuts made to return to 8 pages
Figure widths 0.57–0.72\linewidth (Figs 2–4 at 0.72, Fig. 5 at 0.62, Fig. 6 at 0.57: fonts ≈ 6–7 pt, below the brief's 7 pt for Figs 5–6 — the coder could re-export those two at narrower design widths); Discussion halved; Table 1 collapsed-fit row removed (numbers in prose); "fits 2, 3" row of Table 2 moved to its caption; captions and prose trimmed throughout; the L-BFGS table dropped in favour of the caption. 15 references.

## Remaining issues
- Figs 5 and 6 are set below the 7 pt font floor; everything else above it.
- Sec 7's "twenty times cheaper" compares 250–300 evaluations with 5000–5900 Adam iterations (one evaluation each); the 20 000-iteration Adam end point is 3.7 % lower in loss and is never reached by L-BFGS on this loss.
- The old κ_Si-sign number "90σ" is now 77σ and MAE/σ 15 (R12).
