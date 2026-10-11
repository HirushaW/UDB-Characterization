# Changelog

All notable changes to the UDB Characterization Tools are recorded here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) ·
Versioning: [Semantic Versioning](https://semver.org).

## [Unreleased]

## [1.0.0] - 2026-10-11

First public release of the two companion tools for the Universal Dendrimer Builder.

### UDB Characteriser v1.0 (`udb_characterizer.py`)
- Standalone post-build characterisation with a GUI (single structure, batch, console,
  3-D viewer, summary cards, grouped results, light/dark themes) and a headless command
  line.
- Reads MOL2, AMBER `prmtop` + `inpcrd`/`rst7` (ASCII and NetCDF) and PDB (with `prmtop`,
  with bonds copied from a MOL2, or with bonds inferred by distance).
- Core, monomer and terminal fragment references; terminal element detected from the
  monomer's heavy leaves.
- 117 properties: size, gyration-tensor shape, exact moment of inertia and Ray's
  asymmetry, Shrake–Rupley SASA, grid and Monte Carlo van der Waals volume, void
  fraction, density, Kirkwood–Riseman and hard-sphere hydrodynamic radius with D, s and
  [η], back-folding, charge and dipole, two degree-of-branching measures.
- Quality control: connectivity, tree check, bond-length and clash audits, pass/fail gate.
- Plain-language interpretation with every threshold tagged [lit], [conv] or [!].
- CSV, per-column units file, full JSON report and PNG export.

### Universal Dendrimer Parameterization v2.0 (`dendrimer_param_gui_v2.py`)
- Sybyl → GAFF2 typing (including Si, P and S), `parmchk2` and `tleap` for AMBER, ParmEd
  for NAMD/CHARMM and `acpype` for GROMACS; save-as-script option.
- Validation of counts, atom types, charge sum, bond lengths and `frcmod` penalties.
- Extract mol2: GAFF2 → Sybyl MOL2 (and GAFF2 MOL2 and PDB) for finished dendrimers, from
  the MOL2 alone or with `prmtop`/`inpcrd`.
- Environment check for AmberTools, ParmEd and `acpype`.

### Fixed (relative to the pre-release files)
- Extract mol2 tab: the conversion worker thread was never started, so
  **Convert mol2 / Extract mol2** did nothing. It now runs.
- `eta_iv_hardsphere` was 10× too small (nm³ converted to L instead of dL). For the
  PAMAM G2 example it changes from 0.0112 to 0.1124 dL g⁻¹. The Kirkwood values were
  already correct.
- Command line: `--terminal` was accepted but ignored; the terminal fragment is now
  loaded and reported, as in the GUI.

### Also in this release
- PAMAM G2 worked example with reference outputs of both tools.
- Property reference (`docs/PROPERTIES.md`), parameterization / MOL2 guide, workflow
  diagram, example-results figure and screenshots of both tools.
- `pytest` suite and GitHub Actions CI (Python 3.10 and 3.12, GUI test under xvfb).
- MIT License, citation file and security policy.

[Unreleased]: https://github.com/HirushaW/UDB-Characterization/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/HirushaW/UDB-Characterization/releases/tag/v1.0.0
