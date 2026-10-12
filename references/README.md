# References

Local library of papers and datasets cited in the code and docs.

**PDFs and third-party data files are git-ignored** (copyright). Only this index is committed.
Files keep the names they had in the author's library, so look items up by the File column.

In-code citation format: `[KEY]p{page}e{eq}`, e.g. `[PAP91]p54e39`. Page numbers are the
journal/book page numbers printed on the page, not PDF page indices.
Legacy code uses numeric keys from the 2017 refs page: `[1]` = PAP91, `[2]` = Bastin90, `[3]` = Bastin86.

Status: ✅ in folder · ⬜ not available · 🔓 free copy online.
DOIs on ⬜ items were written from memory and are unverified.

## Core: PAP / XPP and thin-film methods

| Key | Citation | File / source | |
|---|---|---|---|
| PAP91 | Pouchou J-L & Pichoir F (1991). Quantitative analysis of homogeneous or stratified microvolumes applying the model "PAP". In *Electron Probe Quantitation*, Heinrich & Newbury (Eds.), pp. 31–75. Springer. | `Pouchou_Pichoir_1991_PAP.pdf` | ✅ |
| PAP84a | Pouchou J-L & Pichoir F (1984). A new model for quantitative X-ray microanalysis. Part I: homogeneous samples. *La Recherche Aérospatiale* 3, 167–192. | `1984_A NEW MODEL … PART 1.pdf` | ✅ |
| PAP84b | Pouchou J-L & Pichoir F (1984). … Part II: in-depth analysis of heterogeneous samples. *La Recherche Aérospatiale* 5, 349–367. | `1984_A NEW MODEL … PART Il.pdf` | ✅ |
| XPP88 | Pouchou J-L & Pichoir F (1988). A simplified version of the PAP model for matrix corrections in EPMA. *Microbeam Analysis*. | `J. L. Pouchou_F. M. A. Pichoir_1988_A SIMPLIFIED VERSION …pdf` | ✅ |
| Pouchou93 | Pouchou J-L (1993). X-ray microanalysis of stratified specimens. *Anal Chim Acta* 283, 81–97. (NiCr/FeGdPt/Si k-ratios) | `Jean-Louis Pouchou_1993_X-Ray microanalysis of stratified specimens.pdf` | ✅ |
| PP93 | Pouchou J-L & Pichoir F (1993). EPMA applied to thin surface films and stratified specimens. *Scanning Microsc Suppl* 7, 167–189. (modified PAP) | — | ⬜ |
| Pouchou02 | Pouchou J-L (2002). X-ray microanalysis of thin surface films and coatings. *Mikrochim Acta* 138, 133–152. (Al/SiO₂ data) | `Pouchou_2002_X-ray_microanalysis_of_thin_surface_films_and_coatings.pdf` | ✅ |
| Waldo88 | Waldo RA (1988). An iteration procedure to calculate film compositions and thicknesses in EPMA. *Microbeam Analysis*, 310–314. (GMRFilm; "RAW-MAS-88" in legacy code) | `R. A. Waldo_1988_AN ITERATION PROCEDURE …pdf` | ✅ |
| Waldo91 | Waldo RA (1991). A characteristic X-ray fluorescence correction for thin-film analysis by electron microprobe. *Microbeam Analysis*, 45–53. | — (needed for M4; hard copy only) | ⬜ |
| Waldo93 | Waldo RA et al. (1993). Quantitative thin-film analysis with an energy-dispersive X-ray detector. | `Waldo et al_1993_…pdf` | ✅ |
| Merlet93 | Merlet C (1993). Evaluation and description of a new correction procedure for quantitative EPMA (XPHI). | `C Merlet_1993_EVALUATION AND DESCRIPTION …pdf` | ✅ |
| LlovetMerlet10 | Llovet X & Merlet C (2010). EPMA of thin films and multilayers using the computer program XFILM. *Microsc Microanal* 16, 21–32. | `Llovet_Merlet_2010_…pdf` | ✅ |
| Heijligers00a | Heijligers HJM et al. (2000). A systematic database of thin-film measurements by EPMA. Part I: aluminum films. | `H. J. M Heijligers_2000_… Part 1 - Aluminum films.pdf` | ✅ |
| Heijligers00b | … Part II. | `H. J. M Heijligers_2000_… Part II.pdf` | ✅ |

## BadgerFilm and recent work

| Key | Citation | File / source | |
|---|---|---|---|
| Moy21a | Moy A & Fournelle J (2021). φ(ρz) distributions in bulk and thin film samples for EPMA. Part 1. *Microsc Microanal* 27, 266–283. doi:10.1017/S1431927620024915 | `Moy_Fournelle_2021_MM_Part1_phirz_bulk.pdf` 🔓 par.nsf.gov/servlets/purl/10222864 | ✅ |
| Moy21b | Moy A & Fournelle J (2021). … Part 2: BadgerFilm. *Microsc Microanal* 27, 284–296. doi:10.1017/S1431927620024927 | `Moy_Fournelle_2021_MM_Part2_BadgerFilm.pdf` 🔓 par.nsf.gov/servlets/purl/10222866 | ✅ |
| Moy21-SI | Supplementary material for Moy21a/b (analytic fluorescence integrals) | via the DOIs above | ⬜ |
| Llovet21 | Llovet X, Moy A, Pinard PT & Fournelle JH (2021). EPMA: a review of recent developments. *Prog Mater Sci* 116, 100673. | — | ⬜ |
| Ritchie20 | Ritchie NWM (2020). Embracing uncertainty: modeling the standard uncertainty in EPMA — Part I. *Microsc Microanal* 26, 469–483. | — | ⬜ |
| Terborg23 | Terborg R, Kim KJ & Hodoroaba V-D (2023). Composition of thin films by EPMA. *Surf Interface Anal*. | 🔓 opus4.kobv.de/opus4-bam/files/57636/SurfInterfAnal_2023_Terborg_Kim_Hodoroaba_CompoThinFilmsByEMPA.pdf | ⬜ |

## Gaussian φ(ρz) models (Packwood, Bastin)

| Key | Citation | File | |
|---|---|---|---|
| PackwoodBrown81 | Packwood RH & Brown JD (1981). A Gaussian expression to describe φ(ρz) curves. | `PACKWOOD,_RH_BROWN,_JD_1981_…pdf` | ✅ |
| BrownPackwood82 | Brown JD & Packwood RH (1982). Quantitative EPMA using Gaussian φ(ρz) curves. | `J._D._Brown_R._H._Packwood_1982_…pdf` | ✅ |
| Packwood87 | Packwood R et al. (1987). EPMA of compositionally stratified materials. | `R.Packwood et al_1987_…pdf` | ✅ |
| PackwoodGen | Packwood RH. A general equation for predicting X-ray intensities from stratified samples. | `R_H_Packwood_A_General_Equation_…pdf` | ✅ |
| Packwood91 | Packwood RH (1991). A comprehensive theory of EPMA. | `Packwood_1991_…pdf` | ✅ |
| Bastin84 | Bastin GF et al. (1984). Evaluation of the use of Gaussian φ(ρz) curves in quantitative EPMA. | `Bastin_et_al_1984_…pdf` | ✅ |
| Bastin86 | Bastin GF, Heijligers HJM & van Loo FJJ (1986). A further improvement in the Gaussian φ(ρz) approach. *Scanning* 8, 45–67. | — | ⬜ |
| Bastin90 | Bastin GF & Heijligers HJM (1990). Quantitative EPMA of ultralight elements (boron–oxygen). *Scanning* 12, 225–236. | — | ⬜ |
| BastinHeijligers86 | Bastin GF & Heijligers HJM (1986). Quantitative EPMA of carbon in binary carbides (two versions in folder). | `Bastin_Heijligers_1986_…pdf`, `…_2.pdf` | ✅ |
| PROZA96 | Bastin GF et al. (1998). An improved matrix correction (PROZA96). | `Bastin_et_al_1998_…PROZA96.pdf` | ✅ |
| Bastin01 | Bastin GF et al. (2001). Experimental determination of the surface ionization. | `Bastin_et_al_2001_…pdf` | ✅ |
| BastinModified | Bastin GF et al. The performance of the modified φ(ρz) approach compared to Love and Scott… | `G. F. Bastin et al_The Performance …pdf` | ✅ |

## Electron physics, cross sections, atomic data

| Key | Citation | File | |
|---|---|---|---|
| Hutchins74 | Hutchins GA (1974). Electron probe microanalysis. In *Characterization of Solid Surfaces*, 441–484. Plenum. | — | ⬜ |
| RusteGantois75 | Ruste J & Gantois M (1975). A quantitative analysis of very light elements by the electron probe microanalyser. *J Phys D* 8, 872–890. (J_i) | `Ruste_Gantois_1975_…pdf` | ✅ |
| Ruste79 | Ruste J (1979). Principes généraux de la microanalyse quantitative appliquée aux éléments très légers. | `Jacky RUSTE_1979_…pdf` | ✅ |
| SachsRichardson53 | Sachs & Richardson (1953). Mean excitation potentials. | `Sachs_Richardson_1953_…pdf` | ✅ |
| BakkerSegre51 | Bakker & Segrè (1951). Stopping power and energy loss for ion pair production. | `Bakker_Segre_1951_…pdf` | ✅ |
| Bambynek72 | Bambynek W et al. (1972). X-ray fluorescence yields, Auger, and Coster–Kronig transition probabilities. *Rev Mod Phys* 44, 716. (two scans in folder) | `Bambynek et al_1972_…pdf`, `…2.pdf` | ✅ |
| Heinrich87 | Heinrich KFJ (1987). Mass absorption coefficients for EPMA (MAC30). *Proc ICXOM XI*, 67–119. | — | ⬜ |
| PouchouHeijligers88 | Pouchou J-L & Heijligers HJM (1988). Determination of MACs for soft X-rays by use of the electron microprobe. | `J. L. Pouchou_H. J. M Heijligers_1988_…pdf` | ✅ |
| Kyser | Kyser DF. Experimental determination of MACs for soft X-rays. | `D. F. KYSER_…pdf` | ✅ |
| RickerbyWachter00 | Rickerby DG & Wachter N (2000). Effective L-series MACs for EDS. | `David G. Rickerby_Norbert Wachter_2000_…pdf` | ✅ |
| xraylib | Schoonjans T et al. (2011). The xraylib library for X-ray–matter interactions. *Spectrochim Acta B* 66, 776–784. | (software dependency) | — |

## Other EPMA / thin-film papers in the folder

Heinrich 1968 (NBS, *Quantitative EPMA*); Heinrich et al. 1976 (Monte Carlo in EPMA);
Love & Scott 1978; Love et al. 1978; Scott & Love 1990; Philibert 1968; Philibert & Tixier 1968;
Philibert et al. 1970 (two versions); Reuter et al. 1978 (X-ray production range);
Cosslett & Thomas 1964 (I, II) and 1965 (III); Fitting 1974; Mikheev et al. 2010
(fluorescence depth distribution); Isogai 1970; Kerrick et al. 1973 (carbon coat thickness);
Willich 1992; Willich & Schiffmann 1992; Stary & Jurek 2002; Heluani & Hoffmann 2002, 2003;
Pidun et al. 2000; Fournelle, Geology 777 course notes (2003).

## Datasets (`data/`)

| File | Source | Notes |
|---|---|---|
| `data/BadgerFilm_Al_on_SiO2_122nm_exp_data_40TOA.xlsx` | github.com/Aurelien354/BadgerFilm (UW research/edu terms) | Pouchou02 Al film on SiO₂; Al/Si/O Kα k-ratios at 5–30 kV; stds Al, Si, Y₃Fe₅O₁₂. File says TOA 40°, Moy21b text says 35° — check Pouchou02. |
