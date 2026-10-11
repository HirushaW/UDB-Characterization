# Parameterization GUI and MOL2 conversion

`dendrimer_param_gui_v2.py` (Universal Dendrimer Parameterization v2.0) does two jobs
around the Universal Dendrimer Builder:

1. **Before a build** — turns core and monomer fragments drawn in Avogadro (Sybyl atom
   types) into GAFF2-typed MOL2 files and AMBER / NAMD-CHARMM / GROMACS topologies.
2. **After a build** — converts finished UDB dendrimers between GAFF2-typed and Sybyl MOL2
   (Extract mol2 tab), so they can be opened, edited and re-used in other programs.

![Parameterization tab](images/param_gui_parameterization.png)

## Tabs

| Tab | What it does |
|---|---|
| **Parameterization** | Load a Sybyl MOL2, choose the family and whether it is a core or monomer, set the net charge and output formats, then run the pipeline below or save it as a shell script. |
| **Environment Check** | Finds `antechamber`, `parmchk2`, `tleap`, `python3`, ParmEd and `acpype` and prints their versions. |
| **Validation** | Runs automatically after parameterization; shows the checks below and the `frcmod` with penalty lines highlighted. |
| **Extract mol2** | GAFF2 ↔ Sybyl MOL2 conversion of a whole dendrimer (Mode A or Mode B). |
| **References** | Key literature for GAFF2 and the validation sources. |

## Parameterization pipeline (Sybyl → GAFF2 → topology)

The generated script runs five stages:

| Stage | Tool | Output |
|---|---|---|
| 1. Parse and fix MOL2 | Python | Atom names made unique (`C1`, `C2`, …, `N1`, …). |
| 2. GAFF2 atom typing | Python (rules below) | `<name>_typed.mol2` |
| 3. Missing parameters | `parmchk2 -s gaff2` | `<name>.frcmod` (penalty scores printed) |
| 4. AMBER topology | `tleap` (`leaprc.gaff2`) | `<name>.prmtop`, `<name>.inpcrd`, `<name>_final.mol2` |
| 5. Format conversion | ParmEd / `acpype` | NAMD/CHARMM `.psf` + `.pdb`; GROMACS `.top` + `.gro` in `<name>.acpype/` |

Typing is done in Python rather than by `antechamber`, so silicon, phosphorus and sulfur
fragments (carbosilane, phosphorus dendrimers) are typed without antechamber failing on
them:

| Element | GAFF2 type |
|---|---|
| C | `ca` if the Sybyl type is aromatic or sp² (`C.ar`, `C.2` — including carbonyl carbons), else `c3` |
| N | `n3` (all nitrogens, including amide `N.am`) |
| O | `o` if the Sybyl type is `O.2`/`O.co2`, else `os` |
| H | `hn` on N, `ho` on O, `hs` on Si, `hp` on P, else `hc` |
| Si, P, S | `Si`, `p5`, `ss` |
| F, Cl, Br | `f`, `cl`, `br` |
| anything else | `DU` (flagged by validation) |

These rules are coarse: carbonyl carbons become `ca` rather than `c`, and
amide nitrogens `n3` rather than `n`. `parmchk2` supplies the missing terms, so check the
penalty scores on the Validation tab. Charges are carried over from the input MOL2. For the AM1-BCC fragment charges used by
UDB, run the fragments through UDB's own charge step.

### Validation checks

1. Typed MOL2 exists and parses.
2. Atom count and bond count match the input.
3. No `DU` (unassigned) atom types.
4. Charge sum within 0.5 e of the requested net charge.
5. Bond lengths (first 200 bonds) inside reference ranges, e.g. C–C 1.40–1.60 Å,
   Si–C 1.75–1.98 Å, P–O 1.45–1.70 Å.
6. `frcmod` exists and no parameter has a penalty score above 10.

## Extract mol2 (GAFF2 → Sybyl)

![Extract mol2 tab](images/param_gui_extract_mol2.png)

| Mode | Inputs | Coordinates from |
|---|---|---|
| **A — mol2 only** | GAFF2-typed dendrimer MOL2 from UDB | The MOL2 itself |
| **B — mol2 + prmtop + inpcrd** | MOL2 plus the AMBER topology and coordinates | `ambpdb` (falls back to the MOL2 if atom counts differ) |

| Output | File | Use |
|---|---|---|
| Sybyl MOL2 | `<name>_sybyl.mol2` | Avogadro, PyMOL, Open Babel, CGenFF and other Sybyl-type tools. |
| GAFF2 MOL2 | `<name>_gaff2.mol2` | Clean re-write that keeps the GAFF2 types. |
| PDB | `<name>_recovered.pdb` | VMD. Mode A uses `antechamber -fo pdb`; Mode B writes the `ambpdb` output. |

GAFF2 → Sybyl map used for the conversion:

| GAFF2 | Sybyl |
|---|---|
| `c3` | `C.3` |
| `c`, `c2` | `C.2` |
| `ca` | `C.ar` |
| `c1` | `C.1` |
| `n3`, `nt` | `N.3` |
| `n`, `nh` | `N.am` |
| `n2` | `N.2` |
| `na` | `N.ar` |
| `n4` | `N.4` |
| `o` | `O.2` |
| `os`, `oh`, `ow` | `O.3` |
| `p4`, `p5` | `P.3` |
| `ss`, `sh` | `S.3` |
| `s2` | `S.2` |
| `Si`/`si`, `cl`, `br`, `f`, `i` | `Si`, `Cl`, `Br`, `F`, `I` |
| every GAFF2 hydrogen (`hc`, `ha`, `hn`, `ho`, `h1`–`h5`, …) | `H` |

Types not in the map are written unchanged. Atom names, coordinates, substructure
fields and charges are kept; every bond is written with bond order `1`.

## Requirements

- Python 3 with Tkinter.
- AmberTools25 (`antechamber`, `parmchk2`, `tleap`, `ambpdb`, ParmEd) for parameterization,
  Mode B and PDB output. Sybyl ↔ GAFF2 MOL2 conversion in Mode A needs only Python.
- `acpype` for GROMACS output (optional): `conda install -c conda-forge acpype -y`.
