# UDB Characteriser — property reference

Every column the Characteriser writes, how it is computed and in what unit. The same
units are written to `<output>_units.csv` on every run (`--units-row` also puts them as
line 2 of the main CSV).

Conventions used below:

- **Heavy atoms** are all atoms that are not hydrogen (name or GAFF type starting with `h`).
- **Terminal atoms** are heavy atoms of the terminal element (auto-detected from the
  monomer's heavy leaves, or `--terminal-element`) with exactly one heavy neighbour.
- **Centre** for back-folding is the graph centre of the bond tree (no assumption about
  atom order); other quantities use the centroid or the centre of mass as stated.
- Solution properties use water at **298 K**, η = 8.9 × 10⁻⁴ Pa s, and a partial
  specific volume v̄ = 0.73 cm³ g⁻¹.

Every value describes **one conformation**. For publication, characterise an ensemble
(independent builds or MD frames) and report mean ± sd.

## Identity and counts

| Column | Unit | Meaning |
|---|---|---|
| `source`, `family`, `generation`, `state` | text | File name; family and generation from `--family`/`--gen` or the file/folder name; `amber_min` if the name contains `_min`/`minim`/`.rst`, else `as_built`. |
| `input_note` | text | How atoms and bonds were obtained (MOL2 bonds, prmtop bonds, transferred, or **inferred by distance**). |
| `terminal_element`, `terminal_detection` | text | Element used for terminal detection and how it was chosen. |
| `N_atoms`, `N_bonds`, `N_heavy`, `N_H` | count | Atom, bond, heavy-atom and hydrogen counts. |
| `MW_Da` | Da | Σ standard atomic masses. |
| `n_term_detected` | count | Terminal atoms found in the structure. |
| `n_term_theoretical` | count | f · m^G, with f = 4 for PAMAM, otherwise f = m (m = 2). |
| `DB_completeness` | – | `n_term_detected / n_term_theoretical` (1.0 = complete generation). |
| `DB_structural`, `n_dendritic`, `n_linear` | –, count | Classical degree of branching 2D / (2D + L) over heavy atoms; D = ≥3 heavy neighbours, L = 2. |

## Size

| Column | Unit | Method |
|---|---|---|
| `Rg_all_A` | Å | Radius of gyration over all atoms (unweighted, about the centroid). |
| `Rg_heavy_A` | Å | Same over heavy atoms (Maiti 2004 convention). The size to compare with SAXS/SANS. |
| `R_sphere_A`, `diameter_A` | Å | Equivalent solid sphere, R = √(5/3) · Rg_heavy. |
| `ee_max_A`, `ee_mean_A` | Å | Maximum and mean terminal-to-centroid distance. |

## Shape (gyration tensor)

The gyration tensor G = (1/N) Σ rᵢ rᵢᵀ over heavy atoms, eigenvalues λ₁ ≤ λ₂ ≤ λ₃.

| Column | Unit | Formula |
|---|---|---|
| `asphericity` | – | δ = 1 − 3 I₂ / I₁², with I₁ = Σλ and I₂ = λ₁λ₂ + λ₂λ₃ + λ₁λ₃. 0 = sphere, 1 = rod. |
| `kappa2` | – | κ² = (3/2) Tr(G²)/Tr(G)² − 1/2. For a single structure this is algebraically identical to δ; the two differ only when averaged over an ensemble. |
| `acylindricity` | Å² | c = λ₂ − λ₁. |
| `shape_label` | text | δ < 0.02 near-perfect sphere; < 0.05 slightly aspherical; < 0.15 moderately; else highly aspherical. |

## Moment of inertia

Mass-weighted tensor about the centre of mass, I = Σ mᵢ (|rᵢ|² 𝟙 − rᵢ rᵢᵀ). This is
exact (no free parameters) and is a different quantity from the gyration tensor.

| Column | Unit | Meaning |
|---|---|---|
| `I1_amuA2`, `I2_amuA2`, `I3_amuA2` | amu Å² | Principal moments, I₁ ≤ I₂ ≤ I₃. |
| `I_ratio_I1_I3`, `I_ratio_I2_I3` | – | Dimensionless shape ratios. |
| `inertia_shape` | text | Spherical top / prolate / oblate / asymmetric top (5 % tolerance). |
| `Rg_axis1_A` … `Rg_axis3_A` | Å | √(Iₖ / M) for each principal axis. |
| `total_mass_amu` | amu | Σ mᵢ. |
| `rays_kappa` | – | Ray's asymmetry κ = (2B − A − C)/(A − C), A,B,C ∝ 1/I: −1 prolate, 0 asymmetric, +1 oblate. |

Principal axis vectors are in the JSON report.

## Surface

| Column | Unit | Method |
|---|---|---|
| `SASA_A2` | Å² | Shrake–Rupley (1973), probe 1.4 Å, 960 Fibonacci points per atom (`--sasa-points`), GAFF2 radii. |
| `SASA_per_term_A2` | Å² | SASA / terminals. Below ~28.4 Å² per terminal the surface is sterically congested. |

## Volume and free space

| Column | Unit | Method |
|---|---|---|
| `V_sphere_A3` | Å³ | (4/3) π R_sphere³. |
| `V_vdW_grid_A3` | Å³ | **Preferred.** Deterministic voxel occupancy of the GAFF2 vdW spheres, 0.25 Å grid (`--grid-spacing`); error well under 1 %. |
| `V_vdW_mc_A3` | Å³ | Monte Carlo estimate, 100 000 samples (`--mc-samples`), fixed seed; ±1.5 % seed-to-seed. |
| `V_enc_A3` | Å³ | V_sphere − V_vdW (MC): internal cavity estimate. |
| `void_pct_mc` | % | 100 · (1 − V_vdW / V_sphere). A static proxy for guest-loading capacity. |
| `density_g_cm3` | g cm⁻³ | 1.66054 · MW / V_sphere. |

GAFF2 radii come from the AmberTools25 `gaff2.dat` MOD4 RE block; unknown types fall
back to Bondi (1964) radii by element.

## Hydrodynamics and transport

| Column | Unit | Method |
|---|---|---|
| `Rh_kirkwood_A` | Å | **Preferred.** Kirkwood–Riseman rigid-bead estimate over heavy atoms, 1/Rh = (1/N²) Σᵢ≠ⱼ 1/rᵢⱼ. Typically within 5–10 % of HYDROPRO for globular molecules. The quantity measured by DLS. |
| `Rh_hardsphere_A` | Å | √(5/3) · Rg_heavy. Kept for comparison; unreliable for open or branched structures. |
| `Rh_method_diff_pct` | % | 100 · (Rh_hs − Rh_kw) / Rh_kw. |
| `D_kirkwood_m2s` | m² s⁻¹ | Stokes–Einstein, D = k_BT / (6πηRh). |
| `s_kirkwood_Sv` | S | Svedberg, s = M(1 − v̄ρ) / (N_A 6πηRh). |
| `Vh_kirkwood_nm3` | nm³ | (4/3) π Rh³. |
| `eta_iv_kirkwood_cm3g`, `eta_iv_kirkwood_dLg` | cm³ g⁻¹, dL g⁻¹ | Einstein / Flory–Fox, [η] = 2.5 Vh N_A / M. |
| `D_m2s_hardsphere`, `s_Sv_hardsphere`, `eta_iv_hardsphere` | as above | Same formulas from the hard-sphere Rh. |

D and s carry the same relative error as Rh; [η] scales as Rh³ and carries roughly three
times it.

## Back-folding

Distances of terminal atoms from the core centre, with R_ref = √(5/3) · Rg_heavy.

| Column | Unit | Meaning |
|---|---|---|
| `n_terminals`, `R_ref_A` | count, Å | Terminals used and the reference radius. |
| `term_r_mean_A`, `_median_`, `_min_`, `_max_`, `_std_` | Å | Terminal radial distribution. |
| `terminal_r_norm` | – | ⟨r_terminal⟩ / R_ref (1.0 = perfect outer shell). |
| `bf_frac_half_radius_pct` | % | Terminals with r < 0.5 R_ref. |
| `bf_frac_equal_volume_pct` | % | Terminals inside the volume-bisecting radius R_ref / 2^(1/3). |

Radial number-density profiles (terminals and all atoms, 10 shells) are in the JSON
report. On an **as-built** structure these numbers mostly describe the builder, which
places terminals outward by construction; back-folding is meaningful only after MD.

## Charge

From the partial charges in the input. These are **not** a zeta potential.

| Column | Unit | Meaning |
|---|---|---|
| `net_charge_e` | e | Σ qᵢ. |
| `dipole_eA`, `dipole_Debye` | e Å, D | |Σ qᵢ (rᵢ − r_COM)|; 1 e Å = 4.8032 D. |
| `dipole_origin_independent` | boolean | True when |net charge| < 0.001 e. |
| `sum_abs_charge_e` | e | Σ |qᵢ|. |
| `terminal_charge_e` | e | Σ q over terminal atoms. |
| `surface_charge_density_e_per_A2` | e Å⁻² | terminal charge / SASA. |

## Quality control

| Column | Unit | Meaning |
|---|---|---|
| `QC_PASS` | boolean | Connected **and** a tree (N_bonds = N_atoms − 1) **and** ≤ 2 % bad bonds. `--strict` drops failing structures. |
| `connected`, `n_components`, `is_tree` | boolean, count | Graph checks. |
| `centre_atom_id` | atom index | Graph centre of the bond tree. |
| `bond_min_A`, `bond_max_A`, `bond_mean_A` | Å | Bond-length statistics. |
| `n_bad_bonds`, `pct_bad_bonds` | count, % | Bonds outside 0.85–1.80 Å. |
| `min_nonbonded_A`, `worst_overlap_A`, `contacts_lt_2A` | Å, Å, count | Clash audit over pairs ≥ 4 bonds apart within 5 Å. |

## Fragment references

`core_*`, `monomer_*` and `terminalfrag_*` columns record the fragment file, atom, bond
and heavy-atom counts and heavy-leaf elements. `mass_balance_note` reports
(N − N_core) / N_monomer as monomer-equivalents.

## Interpretation columns

With interpretation on (default; `--no-interpret` to omit), `*_meaning` columns give a
plain-language reading of each group. Every threshold is tagged **[lit]** (published),
**[conv]** (this tool's convention) or **[!]** (caution). `--explain` prints the same text
as a report. Interpretations help read the table; they are not results.

## References

- Shrake, A.; Rupley, J. A. *J. Mol. Biol.* **1973**, 79, 351 — SASA.
- Bondi, A. *J. Phys. Chem.* **1964**, 68, 441 — vdW radii.
- Kirkwood, J. G.; Riseman, J. *J. Chem. Phys.* **1948**, 16, 565 — hydrodynamic radius.
- Rudnick, J.; Gaspari, G. *J. Phys. A* **1986**, 19, L191 — asphericity.
- de Gennes, P.-G.; Hervet, H. *J. Phys. Lett.* **1983**, 44, 351 — dense-shell model, g_max.
- Lescanec, R. L.; Muthukumar, M. *Macromolecules* **1990**, 23, 2280 — back-folding.
- Maiti, P. K.; Çağın, T.; Wang, G.; Goddard, W. A. *Macromolecules* **2004**, 37, 6236 — PAMAM size and shape.
