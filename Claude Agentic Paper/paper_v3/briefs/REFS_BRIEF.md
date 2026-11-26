# Reference verification brief

File: `/home/claude/work/paper_v3/refs_original.bib` (39 entries). Produce `/home/claude/work/paper_v3/refs.bib` (corrected) and `/home/claude/work/paper_v3/reviews/refs_verification.md`.

For every entry: verify authors, title, journal, volume, issue, pages (or article number), year and DOI with web search (publisher page, Crossref, arXiv, Google Scholar). Use `WebSearch` and `WebFetch`; for DOIs the Crossref API `https://api.crossref.org/works/<doi>` is authoritative. Fix what is wrong in `refs.bib` (keep the bibtex keys unchanged). In the report give one line per entry: key — status (verified / corrected: what / not found / flagged) — source URL. Entries marked `[check]` in the original need special care. Flag entries that are preprints or that you could not find at all so the authors can decide. Do not invent any bibliographic detail; if a field cannot be verified, leave it and flag it.

Also note which keys are *cited* in `/home/claude/work/paper_v3/main_v2_accepted.tex` (grep `\cite{`), so that uncited entries can be dropped later (do not drop them yourself).
