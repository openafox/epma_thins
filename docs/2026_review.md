# 2026 review: state of epma_thins, BadgerFilm, and a plan

## 1. Prior art: BadgerFilm (Moy & Fournelle, UW–Madison)

- Papers (free NSF public-access PDFs):
  - Part 1 (bulk φ(ρz), fluorescence): https://par.nsf.gov/servlets/purl/10222864 — doi:10.1017/S1431927620024915
  - Part 2 (thin films, BadgerFilm): https://par.nsf.gov/servlets/purl/10222866 — doi:10.1017/S1431927620024927
  - Supplementary material (analytic fluorescence integrals, after Waldo 1991) is at the DOIs above — not yet read.
- Code: https://github.com/Aurelien354/BadgerFilm — VB.NET, Windows/.NET 4.5.2 GUI. Still maintained
  (commits through 2026-05: EPDL2023 MACs, EADL fixes, UI), but changes since 2021 are data/UI, not new models.
- License: UW "research and educational use, not in any product for sale". **Not OSI open source, not
  GPL-compatible — do not copy code.** Implement from the published equations.
- Later group output: Muley, Nachlas, Moy, Fournelle (2024, Appl. Surf. Sci., unsupported thin films with bulk
  standards); Moy, Llovet, Fournelle (2024, metallic-glass films, EPMA+XPS); a 2025 IOP conf. paper on BadgerFilm.
  Not exhaustive.

### What BadgerFilm does
- PAP φ(ρz) exactly as Pouchou & Pichoir 1991 (they could not reconstruct the 1993 modified PAP).
- Hutchins σ(E) *inside* PAP (required — the model was tuned with it); Bote–Salvat σ, EADL relaxation data and
  PENELOPE-2018 MACs only for absolute intensities (they cancel in k-ratios). MAC30 (Heinrich 1987) and EPDL97 selectable.
- Thin films via PAP weighting law → fictitious bulk composition (4 steps, see §3).
- Characteristic secondary fluorescence with the full PAP curve (Waldo 1991 approach); bremsstrahlung fluorescence
  with a per-energy "imaginary element" φ(ρz) plus two empirical α/β corrections fit to PENEPMA (coefficients only in a
  figure). Brems SF is numerically unstable near threshold (they disable it for some fits).
- Levenberg–Marquardt fit of all k-ratios at all voltages at once; optional soft Σc=1 equations per layer.
  Uncertainties propagate k-ratio errors only (they explicitly omit model/atomic-parameter uncertainty).
- Stated future work (2021): other φ(ρz) models (XPP, XPHI), better brems SF. Not visible in the repo yet.

### Where a Python library still adds value
1. Library-first, scriptable, cross-platform (notebooks, batch processing, HyperSpy/exspy ecosystem).
2. A genuinely open license.
3. Pluggable φ(ρz) models (PAP, XPP, XPHI) — their stated but unrealised future work.
4. Full uncertainty budget (model + atomic parameters, Ritchie 2020) — they explicitly don't do this.
5. Tested against the published cases below, in CI.

## 2. Current code state (as of 2026-10)

Runs (Python 3.12 + numpy/scipy): `atomic_element.py`, `film_layer.py`, `mac.py`, `analysis_sample.py` demo, `bethe.py`.
Does not run: `pap.py` (py2 print + logic bugs), `phirzeq.py`, `gmrfilm.py`, `standard.py` (indentation errors,
missing modules `initvars`, `getcfg`, `get_data`, `scale`). No main loop, no k-ratio input, no fluorescence, no tests.

Known bugs:
- `analysis_sample.py:53` — `sin(toa*180/pi)` should be `sin(toa*pi/180)`.
- `mac.py` — values off by orders of magnitude, not by a constant factor (Ti Kα in Ti ≈ 4e-7 vs ≈110 cm²/g).
- `pap.py:32,246`, `bethe.py:23` — U0 uses emitted line energy; must be the edge energy Ec.
- `pap.py:78` — `xray/1000.**p1` precedence; and R0 must integrate E0→Ec (Part 1 eq 29).
- Off-by-one `== len(layers)` / `> len(layers)` checks (`pap.py:224,262`, `analysis_sample.py:294`) — top/buried branches never run.
- `pap.py`: `clog.append[...]`, `layers.index[...]`, dict called as function, `papwt` returns None, `samp.lar_i`, `user_alert` not imported, `c1`/`ci` mixed.
- `atomic_element.py:251` string-vs-int compare (znl always 2); `:191` aliasing mutates `el.shell`.
- `input()` prompts inside model classes; mutable default args; iterator state stored on the sample.

## 3. Answers to the 2017 open questions (from the papers)

| 2017 question | Answer |
|---|---|
| J_i formula source | Zeller, cited by Ruste & Gantois 1975: J_i = 1e-3·Z(10.04 + 8.25 e^(−Z/11.22)) keV (P1 eq 10) |
| R0 with E0^(1+P) − Eq^(1+P)? | Yes, with Eq = Ec (edge): R0 = (1/M) Σ J^(1−Pk) Dk (E0^(Pk+1) − Ec^(Pk+1))/(Pk+1) (P1 eq 29) |
| mparam source | Hutchins 1974 (P1 eq 8) |
| /5, /2, /3 in `ppint1` | Antiderivative of the weighting law p = N(ρz−L)²(ρz−R)²; `ppint1` is exactly that — correct |
| L = −0.4R or −0.4Rx? | Per Moy (P2 p3): step ii Rr=0.5Rx, **Lr=−0.4Rr**; step iii Rr=0.7Rx, **Lr=−0.6Rr**. Code uses −0.4Rx/−0.6Rx — check against PAP 1991 p54 |
| Why `z_mean = (z_mean + zb)/2`? | Step iv: Rm uses Z_av = (Z_b + Z̄)/2 (P2 p3) — RAW was right |
| "F eq 13 in PAP 1991" | Misprinted in the original; use F = R·(1/S)/σ(E0) (P1 eq 17) |

Thin-film weighting procedure (P2 p3): (i) Rx from un-weighted average composition, then iterate weighting with
(L,R)=(−0.4Rx, Rx) until Rx changes <1%; (ii) (−0.4·0.5Rx, 0.5Rx) → Z_b → η, R, φ(0); (iii) (−0.6·0.7Rx, 0.7Rx) →
Z̄, M, J → 1/S → F; (iv) Z_av for Rm. Fictitious C_i = C_i·N·∫layer p dρz (limits clipped at Rx), normalised.

## 4. Validation fixtures available

| # | Source | Content | Use |
|---|---|---|---|
| F1 | P2 Table 1 | 100 nm Al₂O₃ (ρ 3.95) on Cu, 15 kV, Al Kα1: fictitious compositions per step, Rx 686.7→700.3 µg/cm², Z_b 25.94, R 0.8471, φ(0) 1.6933, Z̄ 26.8, 1/S 1.39e-4, F 8.92e-4, Z_av 26.37, Rm 106.7 µg/cm² | Exact unit tests for thin-film parameter weighting |
| F2 | P1 Tables 1–2 | Bulk Al0.398Ti0.602 (wt), Ti Kα and Al Kα k-ratios 7–40 kV: experiment + BadgerFilm, STRATAGem, CalcZAF (PAP/XPP/Armstrong), **GMRFilm**, PENEPMA. MAC30 for all but BadgerFilm/PENEPMA | Bulk forward-model regression (match GMRFilm/CalcZAF-PAP columns) |
| F3 | P1 Table 3 | Olivine k-ratios at 15 kV + PfE/BadgerFilm compositions | Bulk fit with compound standards (needs standard compositions) |
| F4 | Pouchou 2002 via P2 Fig 6; BadgerFilm repo `Al on SiO2 122nm exp data 40TOA.xlsx` | Al film on SiO₂, 5–30 kV, Al/Si/O k-ratios; all codes give ≈122 nm. Paper says TOA 35°, file name says 40° — check | First thin-film fit test |
| F5 | P2 Table 2 (k-ratios in Pouchou 1993, Anal. Chim. Acta 283, 81) | NiCr / FeGdPt / Si bilayer, 20–30 kV, RBS reference | Multilayer fit test |
| F6 | P2 Table 3 (k-ratios only in figures) | Ge on GaAs, 38 and 320 nm, Kα and Lα | Results-only comparison |

## 5. Proposed plan

Rebuild as a tested library; keep the old modules as reference until replaced.

- **M0 — structure.** `epma_thins/` package, pyproject, pytest, CI. Dataclasses for Element/Line, Layer
  (composition, mass thickness, density, fixed flags), Sample, Measurement (element, line, E0, k, σ, standard).
  No `input()` in core. Atomic data via xraylib; MAC backend pluggable with MAC30 first (for parity with F2).
- **M1 — bulk PAP forward model.** Pure functions for M, J, 1/S, R, φ(0), Rx, Rm, Rc, A1/A2/B1, low-U0 branch,
  analytic emitted integral. Pass F2 (PAP/MAC30 columns).
- **M2 — thin-film forward model.** Weighting law + layer integration with overlayer transmission (P2 eqs 1–2). Pass F1 exactly; forward-check F4.
- **M3 — fitting.** `scipy.optimize.least_squares` over compositions + thicknesses, multi-voltage, soft Σc=1,
  covariance → uncertainties. Recover ≈122 nm on F4; then F5.
- **M4 — characteristic secondary fluorescence** (Waldo 1991 + P1/P2 supplements).
- **M5 — differentiators:** XPP model, bremsstrahlung SF, full uncertainty budget, STRATAGem file import, CLI/notebook examples.

Decided (2026-10-12): keep the name `epma_thins` (PyPI `epma-thins`); license AGPL-3.0-or-later
with commercial licenses on request (contributor agreement for substantial contributions); new code
written from the literature only — GMRFilm Fortran removed from the tree, old Python kept in
`legacy/` for reference. Layer order is surface-first (legacy code was substrate-first).

Open: MAC default for non-parity use (xraylib total vs photoelectric vs PENELOPE-2018); MAC30
implementation for parity with F2 needs Heinrich 1987 (not in library — legacy `mac.py` has the
formula but returns wrong values).

### Status
- M0 done on branch `rebuild`: package skeleton, `data/` (elements, lines, MACs via xraylib with
  overrides), `sample.py` (Layer/Sample/Measurement), tests, CI, LICENSE, CITATION.cff.
