# epma_thins

Thin-film and multilayer quantification for electron probe microanalysis (EPMA): determine the
thickness and composition of films on a substrate from measured k-ratios, as a Python library
you can script.

> **Status: pre-alpha, being rebuilt (2026).** The atomic-data and sample layers work; the φ(ρz)
> models and fitting are in progress. See [docs/2026_review.md](docs/2026_review.md) for the plan.

## Goals

- φ(ρz) forward models for bulk and layered samples, starting with PAP (Pouchou & Pichoir 1991),
  with XPP and others to follow
- multi-voltage least-squares fitting of layer compositions and thicknesses, with uncertainties
- pluggable atomic data (xraylib, MAC30, user-supplied MACs)
- every model validated in CI against published data

## Quick look

```python
from epma_thins import Layer, Sample, Measurement, XrayLine

sample = Sample([
    Layer.from_nm("Al2O3", thickness_nm=100, density=3.95),  # film
    Layer({"Cu": 1.0}),                                        # substrate
])
m = Measurement(XrayLine("Al", "Ka1"), beam_energy=15.0, kratio=0.07)
```

## Development

```bash
uv venv && uv pip install -e ".[dev]"
pytest
```

The original 2015–2017 port of GMRFilm lives in `legacy/` for reference only and will be removed
as the new package replaces it.

## Related work

- [BadgerFilm](https://github.com/Aurelien354/BadgerFilm) (Moy & Fournelle 2021): a Windows
  desktop program for thin-film EPMA, free for research and education. epma_thins is an
  independent implementation written from the published literature.
- GMRFilm (Waldo 1988): the DOS program this project originally started as a port of.

## License

[AGPL-3.0-or-later](LICENSE). Commercial licenses (e.g. for embedding in closed-source products)
are available from the author: openafox@gmail.com.
