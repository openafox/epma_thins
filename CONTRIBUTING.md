# Contributing

Help making this a functional tool for EPMA thin-film analysis is very welcome.

- **Cite your equations.** Physics code must cite its source in a comment using the keys in
  [references/README.md](references/README.md), formatted `[KEY]p{page}e{eq}`, e.g. `[PAP91]p54e39`.
- **Test against published numbers.** New models or corrections need tests reproducing values
  from the literature, with the source cited in the test.
- **Write from the literature, not from other programs' source.** Do not copy code from
  BadgerFilm, GMRFilm, or other projects whose licenses are not compatible with AGPL.
- Run `ruff check src tests` and `pytest` before opening a pull request.

## Licensing of contributions

epma_thins is AGPL-3.0-or-later and is also offered under commercial licenses. To keep that
possible, substantial contributions need a contributor license agreement granting the author the
right to relicense them. Small fixes are accepted without one. Please ask before starting a
large change.
