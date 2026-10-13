# References

Local library of papers and datasets cited in the code and docs.

The PDFs and third-party data live in `references/` (git-ignored here), which is its own private
repository, `openafox/epma_thins-references` (copyright: not for public distribution).
To get them: `git clone git@github.com:openafox/epma_thins-references.git references`.

Files are named `KEY_short-title.pdf`, so a citation `[KEY]` maps directly to a file. Variant scans
add a suffix (`_v2`, `_scan2`).

In-code citation format: `[KEY]p{page}e{eq}`, e.g. `[PAP91]p54e39`. Page numbers are the
journal/book page numbers printed on the page, not PDF page indices.
Legacy code uses numeric keys from the 2017 refs page: `[1]` = PAP91, `[2]` = Bastin90, `[3]` = Bastin86.

Status: ✅ in folder · ⬜ not available · 🔓 free copy online.
DOIs on ⬜ items were written from memory and are unverified.

## Core: PAP / XPP and thin-film methods

| Key | Citation | |
|---|---|---|
| PAP91 | Pouchou J-L & Pichoir F (1991). Quantitative analysis of homogeneous or stratified microvolumes applying the model "PAP". In *Electron Probe Quantitation*, Heinrich & Newbury (Eds.), pp. 31–75. Springer. | ✅ |
| PAP84a | Pouchou J-L & Pichoir F (1984). A new model for quantitative X-ray microanalysis. Part I: homogeneous samples. *La Recherche Aérospatiale* 3, 167–192. | ✅ |
| PAP84b | Pouchou J-L & Pichoir F (1984). … Part II: in-depth analysis of heterogeneous samples. *La Recherche Aérospatiale* 5, 349–367. | ✅ |
| XPP88 | Pouchou J-L & Pichoir F (1988). A simplified version of the PAP model for matrix corrections in EPMA. *Microbeam Analysis*. | ✅ |
| Pouchou93 | Pouchou J-L (1993). X-ray microanalysis of stratified specimens. *Anal Chim Acta* 283, 81–97. (NiCr/FeGdPt/Si k-ratios) | ✅ |
| PP93 | Pouchou J-L & Pichoir F (1993). EPMA applied to thin surface films and stratified specimens. *Scanning Microsc Suppl* 7, 167–189. (modified PAP) | ⬜ |
| Pouchou02 | Pouchou J-L (2002). X-ray microanalysis of thin surface films and coatings. *Mikrochim Acta* 138, 133–152. (Al/SiO₂ data) | ✅ |
| Waldo88 | Waldo RA (1988). An iteration procedure to calculate film compositions and thicknesses in EPMA. *Microbeam Analysis*, 310–314. (GMRFilm; "RAW-MAS-88" in legacy code) | ✅ |
| Waldo91 | Waldo RA (1991). A characteristic X-ray fluorescence correction for thin-film analysis by electron microprobe. *Microbeam Analysis 1991* (26th MAS conference), 45–53. 🔓 archive.org/details/1991-mas-proc-san-jose | ⬜ |
| Waldo93 | Waldo RA et al. (1993). Quantitative thin-film analysis with an energy-dispersive X-ray detector. | ✅ |
| Merlet93 | Merlet C (1993). Evaluation and description of a new correction procedure for quantitative EPMA (XPHI). | ✅ |
| LlovetMerlet10 | Llovet X & Merlet C (2010). EPMA of thin films and multilayers using the computer program XFILM. *Microsc Microanal* 16, 21–32. | ✅ |
| Heijligers00a | Heijligers HJM et al. (2000). A systematic database of thin-film measurements by EPMA. Part I: aluminum films. | ✅ |
| Heijligers00b | … Part II. | ✅ |

## BadgerFilm and recent work

| Key | Citation | |
|---|---|---|
| Moy21a | Moy A & Fournelle J (2021). φ(ρz) distributions in bulk and thin film samples for EPMA. Part 1. *Microsc Microanal* 27, 266–283. doi:10.1017/S1431927620024915 🔓 par.nsf.gov/servlets/purl/10222864 | ✅ |
| Moy21b | Moy A & Fournelle J (2021). … Part 2: BadgerFilm. *Microsc Microanal* 27, 284–296. doi:10.1017/S1431927620024927 🔓 par.nsf.gov/servlets/purl/10222866 | ✅ |
| Moy21-SI | Supplementary material for Moy21a/b (analytic fluorescence integrals), via the DOIs above | ⬜ |
| Llovet21 | Llovet X, Moy A, Pinard PT & Fournelle JH (2021). EPMA: a review of recent developments. *Prog Mater Sci* 116, 100673. | ⬜ |
| Ritchie20 | Ritchie NWM (2020). Embracing uncertainty: modeling the standard uncertainty in EPMA — Part I. *Microsc Microanal* 26, 469–483. | ⬜ |
| Terborg23 | Terborg R, Kim KJ & Hodoroaba V-D (2023). Composition of thin films by EPMA. *Surf Interface Anal*. 🔓 opus4.kobv.de/opus4-bam/files/57636/SurfInterfAnal_2023_Terborg_Kim_Hodoroaba_CompoThinFilmsByEMPA.pdf | ⬜ |

## Gaussian φ(ρz) models (Packwood, Bastin)

| Key | Citation | |
|---|---|---|
| PackwoodBrown81 | Packwood RH & Brown JD (1981). A Gaussian expression to describe φ(ρz) curves. | ✅ |
| BrownPackwood82 | Brown JD & Packwood RH (1982). Quantitative EPMA using Gaussian φ(ρz) curves. | ✅ |
| Packwood87 | Packwood R et al. (1987). EPMA of compositionally stratified materials. | ✅ |
| PackwoodGen | Packwood RH. A general equation for predicting X-ray intensities from stratified samples. | ✅ |
| Packwood91 | Packwood RH (1991). A comprehensive theory of EPMA. | ✅ |
| Bastin84 | Bastin GF et al. (1984). Evaluation of the use of Gaussian φ(ρz) curves in quantitative EPMA. | ✅ |
| Bastin86 | Bastin GF, Heijligers HJM & van Loo FJJ (1986). A further improvement in the Gaussian φ(ρz) approach. *Scanning* 8, 45–67. | ⬜ |
| Bastin90 | Bastin GF & Heijligers HJM (1990). Quantitative EPMA of ultralight elements (boron–oxygen). *Scanning* 12, 225–236. | ⬜ |
| BastinHeijligers86 | Bastin GF & Heijligers HJM (1986). Quantitative EPMA of carbon in binary carbides (two versions). | ✅ |
| PROZA96 | Bastin GF et al. (1998). An improved matrix correction (PROZA96). | ✅ |
| Bastin01 | Bastin GF et al. (2001). Experimental determination of the surface ionization. | ✅ |
| BastinModified | Bastin GF et al. The performance of the modified φ(ρz) approach compared to Love and Scott… | ✅ |

## Electron physics, cross sections, atomic data

| Key | Citation | |
|---|---|---|
| Hutchins74 | Hutchins GA (1974). Electron probe microanalysis. In *Characterization of Solid Surfaces*, 441–484. Plenum. | ⬜ |
| RusteGantois75 | Ruste J & Gantois M (1975). A quantitative analysis of very light elements by the electron probe microanalyser. *J Phys D* 8, 872–890. (J_i) | ✅ |
| Ruste79 | Ruste J (1979). Principes généraux de la microanalyse quantitative appliquée aux éléments très légers. (French) | ✅ |
| SachsRichardson53 | Sachs & Richardson (1953). Mean excitation potentials. | ✅ |
| BakkerSegre51 | Bakker & Segrè (1951). Stopping power and energy loss for ion pair production. | ✅ |
| Bambynek72 | Bambynek W et al. (1972). X-ray fluorescence yields, Auger, and Coster–Kronig transition probabilities. *Rev Mod Phys* 44, 716. (two scans) | ✅ |
| Heinrich87 | Heinrich KFJ (1987). Mass absorption coefficients for EPMA (MAC30). *Proc. 11th ICXOM*, Brown & Packwood (Eds.), Univ. Western Ontario, 67–119. OCLC 222411343. Interlibrary loan request pending. | ⬜ |
| EPQ | NIST EPQ library (DTSA-II engine), public domain. github.com/usnistgov/EPQ — `MassAbsorptionCoefficient.java` (MAC30, PAP91 MACs), `PAP1991.java`, `XPP1991.java` | 🔓 |
| PouchouHeijligers88 | Pouchou J-L & Heijligers HJM (1988). Determination of MACs for soft X-rays by use of the electron microprobe. | ✅ |
| Kyser | Kyser DF. Experimental determination of MACs for soft X-rays. | ✅ |
| RickerbyWachter00 | Rickerby DG & Wachter N (2000). Effective L-series MACs for EDS. | ✅ |
| xraylib | Schoonjans T et al. (2011). The xraylib library for X-ray–matter interactions. *Spectrochim Acta B* 66, 776–784. (software dependency) | — |

## Other EPMA background

| Key | Citation | |
|---|---|---|
| Heinrich68 | Heinrich KFJ (Ed.) (1968). *Quantitative Electron Probe Microanalysis*. NBS Special Publication 298. | ✅ |
| Heinrich76 | Heinrich KFJ et al. (1976). Use of Monte Carlo calculations in EPMA and scanning electron microscopy. | ✅ |
| LoveScott78 | Love G & Scott VD (1978). Evaluation of a new correction procedure for quantitative EPMA. | ✅ |
| Love78 | Love G et al. (1978). A versatile atomic number correction for EPMA. | ✅ |
| ScottLove90 | Scott VD & Love G (1990). The prospects of a universal correction procedure for quantitative EPMA. | ✅ |
| Philibert68 | Philibert J (1968). Quantitative analysis of thin samples. | ✅ |
| PhilibertTixier68 | Philibert J & Tixier R (1968). Electron penetration and the atomic number correction in EPMA. | ✅ |
| Philibert70 | Philibert J et al. (1970). EPMA on electron microscope thin foils (two versions). | ✅ |
| Reuter78 | Reuter W et al. (1978). X-ray production range in solids by 2–15 keV electrons. | ✅ |
| CosslettThomas64a/b, 65 | Cosslett VE & Thomas RN (1964, 1965). Multiple scattering of 5–30 keV electrons in evaporated metal films, I–III. | ✅ |
| Fitting74 | Fitting HJ (1974). Transmission, energy distribution and SE excitation of fast electrons in thin films. | ✅ |
| Mikheev10 | Mikheev NN et al. (2010). An X-ray fluorescence depth distribution function for electron beam excitation. | ✅ |
| Isogai70 | Isogai A (1970). Application of EPMA to heterogeneous structure. | ✅ |
| Kerrick73 | Kerrick DM et al. (1973). The role of carbon film thickness in electron microprobe analysis. | ✅ |
| Willich92 | Willich P (1992). EPMA — a versatile technique for the characterization of thin films and layered structures. | ✅ |
| WillichSchiffmann92 | Willich P & Schiffmann K (1992). EPMA of surface oxide films. | ✅ |
| StaryJurek02 | Starý V & Jurek K (2002). X-ray emission from thin films on a substrate — calculation and experiments. | ✅ |
| HeluaniHoffmann02 | Heluani SP & Hoffmann C (2002). Theoretical modelling of X-ray production in thin-film characterization. | ✅ |
| HeluaniHoffmann03 | Heluani SP & Hoffmann C (2003). Thin-film characterization with X-ray microanalysis. | ✅ |
| Pidun00 | Pidun M et al. (2000). Comparative analysis of a solar control coating on glass by AES, EPMA, SNMS… | ✅ |
| Fournelle03 | Fournelle J (2003). Geology 777 course notes: thin film analysis. Univ. Wisconsin. | ✅ |

## Datasets (`references/data/`)

| File | Source | Notes |
|---|---|---|
| `references/data/Pouchou02_Al-on-SiO2_kratios_from-BadgerFilm.xlsx` | github.com/Aurelien354/BadgerFilm (UW research/edu terms) | Pouchou02 Al film on SiO₂; Al/Si/O Kα k-ratios at 5–30 kV; stds Al, Si, Y₃Fe₅O₁₂. File says TOA 40°, Moy21b text says 35° — check Pouchou02. |
