# Contributing

Thank you for helping improve the UDB Characterization Tools. Bug reports, feature ideas
and pull requests are all welcome.

## Reporting a problem

[Open an issue](https://github.com/HirushaW/UDB-Characterization/issues) with:

- which tool (`udb_characterizer.py` or `dendrimer_param_gui_v2.py`) and version;
- the command you ran or the steps in the GUI;
- the input file if you can share it (a small MOL2 is ideal), and the full error output.

Security problems go through the private route in [SECURITY.md](SECURITY.md), not an issue.

## Making a change

1. Branch from `develop`: `git checkout -b fix/short-name develop` (or `feature/…`).
2. Keep each tool a single self-contained file, as it is now.
3. Add or update a test in `tests/` for the behaviour you change, and run:
   ```bash
   pip install -r requirements.txt pytest
   xvfb-run -a python3 -m pytest      # or plain `pytest`; the GUI test then skips
   ```
4. If you change a reported property or its unit, update `UNITS` in the code,
   [docs/PROPERTIES.md](docs/PROPERTIES.md) and the example outputs in `examples/pamam_g2/`.
5. Add a line under **Unreleased** in [CHANGELOG.md](CHANGELOG.md).
6. Open a pull request into `develop`. CI must be green before it is merged.

Releases are made by merging `develop` into `main` and tagging the version (`v1.1.0`, …).
