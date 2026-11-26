# Reference verification report

Date: 2026-11-23. Input: `refs_original.bib` (39 entries). Output: `refs.bib` (same 39 keys).
Method: Crossref API (`https://api.crossref.org/works/<doi>`) for every journal article; publisher / proceedings pages (JMLR, NeurIPS, ICML/dblp, Elsevier, Wiley) and web search for the rest. DOIs were added to every article entry that has one. Only fields that were verified were changed; nothing was invented.

Summary: 24 verified unchanged (DOI/url added where available), 9 corrected, 1 not found / flagged (Greenall2017), 5 flagged for author decision (Kingma2017, Mockus1978, Li2018, Snoek2012, Greenall2017).

## Per-entry status

| Key | Status | Source |
|---|---|---|
| Yeh1988 | verified (Wiley, 1988) | https://en.wikipedia.org/wiki/Optical_Waves_in_Layered_Media ; https://ci.nii.ac.jp/ncid/BA04528813?l=en |
| Duvillaret1996 | verified; DOI added | https://api.crossref.org/works/10.1109/2944.571775 |
| Dorney2001 | verified; DOI added | https://api.crossref.org/works/10.1364/JOSAA.18.001562 |
| Pupeza2007 | verified; DOI added | https://api.crossref.org/works/10.1364/OE.15.004335 |
| Peretti2019 | verified (vol 9, no 2, pp 136-149, 2019 confirmed); `[check volume/pages]` note removed; DOI added | https://api.crossref.org/works/10.1109/TTHZ.2018.2889227 |
| Zhang2024 | corrected: full author list (Zhang, Ren, Li, D.-D. Zhang, Gu, Xue, Chen) replaces "and others"; DOI added; vol 147, art. 103216, 2024 confirmed | https://nufind.nu.edu.sa/EdsRecord/edselp,S0963869524001816 (ScienceDirect record, DOI 10.1016/j.ndteint.2024.103216) |
| Kingma2017 | corrected: year 2017 -> 2015; note now "ICLR 2015; arXiv:1412.6980" (2017 is only the date of arXiv v9). Key left unchanged -- authors may wish to rename to Kingma2015 | https://arxiv.org/abs/1412.6980 |
| Paszke2019 | corrected: full 21-author list replaces "and others"; pages 8024-8035 added | https://researchr.org/publication/PaszkeGMLBCKLGA19/bibtex ; https://papers.nips.cc/paper/9015-pytorch |
| Li2018 | verified (NeurIPS 31, 2018). Pages not added: dblp/researchr give 6391-6401, the printed proceedings use a different pagination; authors may add if their style requires pages | https://researchr.org/publication/Li0TSG18 ; https://papers.nips.cc/paper/2018/hash/a41b3bb3e6b050b6c9067c67f663b915-Abstract.html |
| Rasmussen2006 | verified (MIT Press, 2006) | well-known monograph; not changed |
| Mockus1978 | corrected: entry type @article -> @incollection; booktitle "Towards Global Optimization 2", editors L. C. W. Dixon and G. P. Szegö, North-Holland, Amsterdam, 1978 (book record verified). Pages 117-129 are as universally cited but were not checked against the book itself -- flagged | https://vufind.lboro.ac.uk/Record/174716/Similar?print=1 (book record) |
| Srinivas2010 | corrected: pages 1015-1022 and publisher Omnipress added | https://bibtex.github.io/ICML-2010-SrinivasKKS.html ; https://icml.cc/imls/conferences/2010/papers/422.pdf |
| Snoek2012 | verified (NIPS 25, 2012). Pages not added: dblp gives 2960-2968, NIPS site gives 2951-2959 (two paginations) -- authors decide | https://researchr.org/publication/SnoekLA12 ; https://papers.nips.cc/paper/4522-practical-bayesian- |
| Grischkowsky1990 | verified; DOI added | https://api.crossref.org/works/10.1364/JOSAB.7.002006 |
| Duvillaret1999 | verified; DOI added | https://api.crossref.org/works/10.1364/AO.38.000409 |
| Dai2004 | verified; DOI added | https://api.crossref.org/works/10.1364/JOSAB.21.001379 |
| Jepsen2005 | verified; DOI added | https://api.crossref.org/works/10.1364/OL.30.000029 |
| Naftaly2007 | verified; DOI added | https://api.crossref.org/works/10.1109/JPROC.2007.898835 |
| Withayachumnankul2008 | verified; DOI added | https://api.crossref.org/works/10.1364/JOSAB.25.001059 |
| Naftaly2009 | verified; DOI added | https://api.crossref.org/works/10.1364/OL.34.001213 |
| Scheller2009 | verified; DOI added | https://api.crossref.org/works/10.1016/j.optcom.2008.12.061 |
| Scheller2014 | corrected: issue number 8 added; DOI added | https://api.crossref.org/works/10.1007/s10762-014-0053-4 |
| Withayachumnankul2014 | verified (print issue vol 35(8), 2014; online Dec 2013); DOI added | https://api.crossref.org/works/10.1007/s10762-013-0042-z |
| Yasui2005 | corrected: second author is Takashi Yasuda (not Tsutomu); DOI added | https://api.crossref.org/works/10.1364/AO.44.006849 |
| Krimi2016 | verified; DOI added | https://api.crossref.org/works/10.1063/1.4955407 |
| Naftaly2019 | verified; DOI added | https://api.crossref.org/works/10.3390/s19194203 |
| Neu2018 | verified; DOI added | https://api.crossref.org/works/10.1063/1.5047659 |
| Greenall2017 | NOT FOUND / FLAGGED. No journal article with this title or this 10-author list exists in Crossref (searched by title and by authors Greenall/Valavanis/Burnett), on White Rose Research Online, or in JIMTW vol 38. The real Greenall items are: PhD thesis "Parameter Extraction and Uncertainty in Terahertz Time-Domain Spectroscopic Measurements" (Univ. Leeds, 2017, https://etheses.whiterose.ac.uk/id/eprint/19045/); IRMMW-THz 2015 paper "Accurate material parameter extraction from broadband terahertz spectroscopy" (10.1109/IRMMW-THz.2015.7327467); IRMMW-THz 2016 "Multilayer extraction of complex refractive index ..." (10.1109/IRMMW-THz.2016.7758468); IRMMW-THz 2017 "Time-frequency extraction of sample thickness ..." (10.1109/IRMMW-THz.2017.8066976); and JIMTW 38(3):325-338 (2017) "The development of a Semtex-H simulant ..." (10.1007/s10762-016-0336-z), which is on a different topic. Entry left as in the original with a FLAGGED note; authors must replace or drop it | https://api.crossref.org/works?query.author=Greenall+Valavanis+Burnett&query.bibliographic=terahertz+parameter+extraction&rows=5 |
| Baydin2018 | verified; url added | https://jmlr.org/papers/v18/17-468.bib |
| Liu1989 | verified; issue 1-3 and DOI added | https://api.crossref.org/works/10.1007/BF01589116 |
| Nocedal2006 | verified (Springer, 2nd ed., 2006) | well-known monograph; not changed |
| Shahriari2016 | verified; DOI added | https://api.crossref.org/works/10.1109/JPROC.2015.2494218 |
| Frazier2018 | verified (arXiv:1807.02811, 2018); preprint, never journal-published -- fine to cite as arXiv | https://arxiv.org/abs/1807.02811v1 |
| Kay1993 | verified (Prentice Hall, 1993) | well-known monograph; not changed |
| Tarantola2005 | verified (SIAM, 2005) | well-known monograph; not changed |
| Aster2018 | verified (3rd ed., Elsevier, 2018, ISBN 978-0-12-804651-7) | https://shop.elsevier.com/books/parameter-estimation-and-inverse-problems/aster/978-0-12-804651-7 |
| Goodfellow2015 | corrected: entry type @article (with venue in `journal`) -> @inproceedings with booktitle "Int. Conf. on Learning Representations (ICLR)", 2015; arXiv note kept | https://research.google/pubs/qualitatively-characterizing-neural-network-optimization-problems/ ; https://arxiv.org/abs/1412.6544v4 |
| Jolliffe2002 | verified (Springer, 2nd ed., 2002) | well-known monograph; not changed |
| BornWolf1999 | verified (CUP, 7th ed., 1999) | well-known monograph; not changed |

Notes
- Optica (OSA) Crossref records give only the first page; the end pages in the bib (Dorney 1571, Pupeza 4350, Grischkowsky 2015, Duvillaret 415, Dai 1386, Jepsen 31, Withayachumnankul 1072, Naftaly 1215, Yasui 6856) agree with the standard citations and were left as is.
- Preprint-only entries: Kingma2017 (ICLR 2015 / arXiv), Frazier2018 (arXiv), Goodfellow2015 (ICLR 2015 / arXiv). Kingma2017 and Goodfellow2015 were peer-reviewed at ICLR.
- Eight monographs (Yeh, Rasmussen, Nocedal, Kay, Tarantola, Aster, Jolliffe, Born & Wolf) were checked against publisher/catalogue data or standard knowledge; no changes.

## Citation usage in `main_v2_accepted.tex`

All 39 keys are cited (grep `\cite{`): Aster2018, Baydin2018, BornWolf1999, Dai2004, Dorney2001, Duvillaret1996, Duvillaret1999, Frazier2018, Goodfellow2015, Greenall2017, Grischkowsky1990, Jepsen2005, Jolliffe2002, Kay1993, Kingma2017, Krimi2016, Li2018, Liu1989, Mockus1978, Naftaly2007, Naftaly2009, Naftaly2019, Neu2018, Nocedal2006, Paszke2019, Peretti2019, Pupeza2007, Rasmussen2006, Scheller2009, Scheller2014, Shahriari2016, Snoek2012, Srinivas2010, Tarantola2005, Withayachumnankul2008, Withayachumnankul2014, Yasui2005, Yeh1988, Zhang2024.

Keys not cited: none.
