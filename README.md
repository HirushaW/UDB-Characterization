# UDB Characterization Tools

**Version 1.0.0** · Python (Tkinter GUI + command line) · Linux · [MIT License](LICENSE)

Two companion programs for the
[Universal Dendrimer Builder (UDB)](https://github.com/HirushaW/universal-dendrimer-builder):

| Tool | File | What it does |
|---|---|---|
| **UDB Characteriser** v1.0 | `udb_characterizer.py` | Characterises dendrimers after UDB has built them (as built, after AMBER minimisation or from MD frames): size, shape, surface, volume, hydrodynamics, back-folding, charge and structure quality. |
| **Universal Dendrimer Parameterization** v2.0 | `dendrimer_param_gui_v2.py` | Converts MOL2 files between Sybyl and GAFF2 atom types and builds AMBER, NAMD/CHARMM and GROMACS topologies for core and monomer fragments. |

![UDB Characteriser](docs/images/udb_characteriser.png)

---

## Where they fit

```
Avogadro                 Parameterization GUI         Universal Dendrimer       UDB Characteriser
core / monomer   ──►   Sybyl ─► GAFF2 MOL2      ──►   Builder: build G1–G10  ──► 117 properties,
(Sybyl MOL2)           + prmtop / psf / top           + AMBER minimisation      QC and interpretation
                                                               │
                                                               ▼
                                                  Parameterization GUI, Extract mol2:
                                                  GAFF2 dendrimer ─► Sybyl MOL2 / PDB
```

## Features

### UDB Characteriser

- **One file, no other UDB files needed.** GUI by default, headless command line for
  clusters and scripts (`numpy` and `scipy` only).
- **Inputs:** MOL2; AMBER `prmtop` + `inpcrd`/`rst7` (ASCII or NetCDF); PDB with a
  `prmtop`, with bonds copied from a MOL2, or with bonds inferred by distance.
- **Fragment references:** load the core and monomer used for the build; the terminal
  element is detected from the monomer and the build is checked against core + monomers.
- **Properties (117 columns):**
  - size: Rg (all and heavy atoms), equivalent sphere and end-to-end distances;
  - shape: asphericity, acylindricity and κ² from the gyration tensor, plus the exact
    mass-weighted inertia tensor and Ray's asymmetry;
  - surface and volume: Shrake–Rupley SASA, deterministic grid van der Waals volume,
    void fraction and density;
  - hydrodynamics: Kirkwood–Riseman Rh, with diffusion coefficient, sedimentation
    coefficient and intrinsic viscosity;
  - back-folding: terminal radial distribution and radial density profiles;
  - charge: net charge, dipole moment and surface charge density.
- **Quality control:** connectivity, tree check, bond-length and steric-clash audits, a
  pass/fail gate (`--strict` drops failures).
- **Interpretation layer:** a plain-language reading of every group, with each threshold
  tagged as literature (**[lit]**), tool convention (**[conv]**) or caution (**[!]**).
- **Batch mode:** characterise a whole folder tree into one CSV.
- **Outputs:** CSV, a units file for every column, full JSON report, PNG of the 3-D view.

Every property, formula and unit: [docs/PROPERTIES.md](docs/PROPERTIES.md).

### Universal Dendrimer Parameterization (MOL2 converter)

- **Sybyl → GAFF2** typing that also handles Si, P and S (carbosilane and phosphorus
  dendrimers), then `parmchk2` and `tleap` for AMBER `.prmtop`/`.inpcrd`, ParmEd for
  NAMD/CHARMM `.psf`/`.pdb` and `acpype` for GROMACS `.top`/`.gro`.
- **Validation:** atom and bond counts, unassigned types, charge sum, bond lengths and
  `frcmod` penalty scores.
- **GAFF2 → Sybyl** conversion of finished UDB dendrimers (Extract mol2 tab), from the
  MOL2 alone or with the AMBER `prmtop`/`inpcrd`, plus PDB output for VMD.
- Environment check for AmberTools, ParmEd and `acpype`.
- Families: PAMAM, PPI, Triazine, Carbosilane, Phosphorus, Polyester, Polyglycerol, PLL.

Details and the full type maps: [docs/PARAMETERIZATION_AND_MOL2.md](docs/PARAMETERIZATION_AND_MOL2.md).

## Repository contents

| Path | What it is |
|---|---|
| `udb_characterizer.py` | UDB Characteriser (GUI and command line). |
| `dendrimer_param_gui_v2.py` | Parameterization GUI and MOL2 converter. |
| `examples/pamam_g2/` | PAMAM G2 inputs and reference outputs of both tools ([guide](examples/README.md)). |
| `docs/` | Property reference, parameterization/MOL2 guide and screenshots. |
| `tests/` | `pytest` suite (also run by GitHub Actions on every push and pull request). |
| `requirements.txt` | Python packages. |
| `CHANGELOG.md` · `CITATION.cff` · `LICENSE` · `SECURITY.md` | Changes, how to cite, licence, security policy. |

## Requirements

- Linux (Ubuntu 20.04 or newer), x86_64; a display for the GUIs
- Python 3.10 or newer with Tkinter (`sudo apt-get install python3-tk`)
- `numpy`, `scipy` (both tools); `matplotlib` (Characteriser GUI only)
- For parameterization: [AmberTools25](https://ambermd.org/AmberTools.php) (`antechamber`,
  `parmchk2`, `tleap`, `ambpdb`, ParmEd); optionally `acpype` for GROMACS

The Characteriser and Mode A of the MOL2 conversion need only Python.

## Installation

```bash
git clone https://github.com/HirushaW/UDB-Characterization.git
cd UDB-Characterization

# Use the same environment as UDB (AmberTools25), or any Python 3.10+
conda activate AmberTools25
pip install -r requirements.txt
```

## Usage

### UDB Characteriser — GUI

```bash
python3 udb_characterizer.py
```

1. **Fragment references (optional, recommended):** load the core and monomer MOL2 used
   for the build. The terminal element is picked from the monomer.
2. **Structure:** load the dendrimer (MOL2, PDB or AMBER restart; add the `prmtop` for
   PDB/restart files).
3. **Options:** family and generation (read from the file name if left blank) and the
   grid spacing for the volume.
4. **Characterise.** Results appear as summary cards and a grouped table; use
   **Explain results**, **Export CSV** or **Export JSON**. The **Batch** tab runs a folder.

### UDB Characteriser — command line

```bash
# one structure, with fragment references and a readable report
python3 udb_characterizer.py dendrimer_pamam_gen3_min_gaff.mol2 \
    --core core_typed.mol2 --monomer PAMAM_monomer_typed.mol2 --explain

# AMBER minimised structure
python3 udb_characterizer.py --prmtop dendrimer_pamam_gen3.prmtop \
    --coords dendrimer_pamam_gen3_min.rst --monomer PAMAM_monomer_typed.mol2

# every structure under a folder, one CSV plus the full JSON report
python3 udb_characterizer.py -r ./Structures/ --monomer PAMAM_monomer_typed.mol2 \
    -o all_results.csv --json all_results.json
```

| Option | Default | Meaning |
|---|---|---|
| `--core`, `--monomer`, `--terminal` | – | Fragment MOL2 references |
| `--prmtop`, `--coords` | – | AMBER topology and coordinates |
| `--bonds-from` | – | MOL2 whose bonds are copied to a PDB with the same atom order |
| `--terminal-element` | from monomer | Override the terminal element |
| `--family`, `--gen` | from file name | Override family and generation |
| `-o`, `--json` | `udb_post_results.csv` | Output CSV; optional full JSON report |
| `--grid-spacing` | 0.25 Å | Grid for the van der Waals volume |
| `--sasa-points`, `--mc-samples` | 960, 100000 | SASA sphere points; Monte Carlo volume samples |
| `--explain` / `--no-interpret` | – | Print the interpretation / leave it out of the CSV |
| `--units-row` | – | Put units on line 2 of the CSV (`*_units.csv` is always written) |
| `--no-extended`, `--strict` | – | Skip extended properties; drop structures that fail QC |

### Parameterization GUI and MOL2 conversion

```bash
python3 dendrimer_param_gui_v2.py
```

- **Sybyl → GAFF2 + topologies:** on the Parameterization tab, load the core or monomer
  MOL2 from Avogadro, choose the family, fragment, net charge and formats, and press
  **Run Parameterization** (or **Save Script Only**). Results are checked on the
  Validation tab.
- **GAFF2 → Sybyl:** on the Extract mol2 tab, load the UDB dendrimer MOL2 (Mode A) or the
  MOL2 with its `prmtop` and `inpcrd` (Mode B), tick the outputs and press **Convert mol2**.

### Worked example

[`examples/pamam_g2`](examples/README.md) has a minimised PAMAM G2 with its core and
monomer, the Characteriser results and the Sybyl conversion, with the commands to
reproduce them.

## Tests

```bash
pip install pytest
python3 -m pytest              # GUI test is skipped without a display
xvfb-run -a python3 -m pytest  # run everything headless
```

## Versions and branches

| Branch | Purpose |
|---|---|
| `main` | Stable code. Every release is tagged here (`v1.0.0`, `v1.1.0`, …). |
| `develop` | Ongoing work. Changes go into `develop` and are merged into `main` for the next release. |
| `release/x.y.z`, `fix/…`, `feature/…` | Short-lived branches merged through pull requests. |

Versions follow [Semantic Versioning](https://semver.org): new features raise the middle
number (1.1.0), fixes raise the last (1.0.1). See [CHANGELOG.md](CHANGELOG.md).

## Citation

A paper describing UDB is in preparation. Until it is published, please cite the software
(see [CITATION.cff](CITATION.cff), or use **Cite this repository** on GitHub):

```bibtex
@software{walgampaya2026udbtools,
  author  = {Walgampaya, Hirusha},
  title   = {{UDB Characterization Tools: Structural Characterisation and
             MOL2 Parameterization Companions to the Universal Dendrimer Builder}},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/HirushaW/UDB-Characterization},
  note    = {Citation pending peer-reviewed publication}
}
```

## License

Released under the [MIT License](LICENSE). AmberTools, ParmEd, acpype and the other
programs these tools call are separate software under their own licences.

## Author

**Hirusha Walgampaya**. Questions and bug reports:
[open an issue](https://github.com/HirushaW/UDB-Characterization/issues).

These tools are independent and not affiliated with AmberTools, NAMD, GROMACS or
CHARMM-GUI; they call these packages as external dependencies.
