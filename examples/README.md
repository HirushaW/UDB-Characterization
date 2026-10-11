# Examples

## `pamam_g2/` — PAMAM generation 2

A PAMAM G2 dendrimer built with the Universal Dendrimer Builder and minimised with
AMBER (516 atoms, 16 terminal amines). It comes from `examples/Structures.tar.gz` in the
[UDB repository](https://github.com/HirushaW/universal-dendrimer-builder).

| File | What it is |
|---|---|
| `input/dendrimer_pamam_gen2_min_gaff.mol2` | Minimised G2 structure, GAFF2 atom types, AM1-BCC charges. |
| `input/PAMAM_core_typed.mol2` | PAMAM core fragment used for the build (GAFF2, 84 atoms). |
| `input/PAMAM_monomer_typed.mol2` | PAMAM monomer used for the build (GAFF2). |
| `input/PAMAM_monomer_sybyl.mol2` | The same monomer as drawn in Avogadro (Sybyl types), the input to the parameterization GUI. |
| `output/pamam_g2_characterisation.csv` | One row, 117 properties, from the UDB Characteriser. |
| `output/pamam_g2_characterisation_units.csv` | Unit of every column in the CSV. |
| `output/pamam_g2_characterisation.json` | Full report, including radial back-folding profiles and principal axes. |
| `output/pamam_g2_interpretation.txt` | Plain-language reading of the results (`--explain`). |
| `output/pamam_g2_sybyl.mol2` | The G2 structure converted from GAFF2 to Sybyl types (Extract mol2 tab, Mode A). |

### Reproduce the characterisation

```bash
cd examples/pamam_g2
python3 ../../udb_characterizer.py input/dendrimer_pamam_gen2_min_gaff.mol2 \
    --core input/PAMAM_core_typed.mol2 \
    --monomer input/PAMAM_monomer_typed.mol2 \
    --family PAMAM --gen 2 \
    -o output/pamam_g2_characterisation.csv \
    --json output/pamam_g2_characterisation.json --explain
```

Headline values:

| Property | Value |
|---|---|
| Atoms / heavy atoms | 516 / 228 |
| Molecular weight | 3256.24 Da |
| Terminal groups (found / expected) | 16 / 16 |
| Rg (heavy atoms) | 13.92 Å |
| Rh (Kirkwood–Riseman) | 11.20 Å |
| Rh (hard sphere, for comparison) | 17.97 Å |
| SASA | 3463.2 Å² |
| van der Waals volume (grid) | 3865.9 Å³ |
| Void fraction | 83.95 % |
| Asphericity δ | 0.512 |
| Quality gate | passed (0 % bad bonds, connected, tree) |

The structure is a single minimised conformation, so these numbers describe that
conformation only. See the *overall caveat* in the interpretation file.

### Reproduce the MOL2 conversion

Open `dendrimer_param_gui_v2.py`, go to **Extract mol2**, choose **Mode A**, select
`input/dendrimer_pamam_gen2_min_gaff.mol2`, tick **Sybyl mol2** and press
**Convert mol2**. The GAFF2 types are mapped as follows in this structure:

| GAFF2 | c | c3 | h1, hc, hn | n | n3 | o |
|---|---|---|---|---|---|---|
| Sybyl | C.2 | C.3 | H | N.am | N.3 | O.2 |
