"""Tests for udb_characterizer.py (engine and command line; no GUI needed)."""
import csv
import math
import os
import subprocess
import sys

import numpy as np
import pytest

import udb_characterizer as uc
from conftest import EXAMPLE, ROOT

STRUCT = os.path.join(EXAMPLE, 'dendrimer_pamam_gen2_min_gaff.mol2')
CORE = os.path.join(EXAMPLE, 'PAMAM_core_typed.mol2')
MONO = os.path.join(EXAMPLE, 'PAMAM_monomer_typed.mol2')


def atom(i, x, y=0.0, z=0.0, name='C', typ='c3', q=0.0):
    return {'id': i, 'name': name, 'type': typ, 'x': x, 'y': y, 'z': z, 'charge': q}


# ── analytic checks ────────────────────────────────────────────────────────────
def test_radius_of_gyration_two_points():
    assert uc.radius_of_gyration([atom(1, -1.0), atom(2, 1.0)]) == pytest.approx(1.0)


def test_asphericity_limits():
    assert uc.asphericity([1.0, 1.0, 1.0]) == pytest.approx(0.0)   # sphere
    assert uc.asphericity([0.0, 0.0, 1.0]) == pytest.approx(1.0)   # rod


def test_kirkwood_two_beads():
    # 1/Rh = (1/N^2) * sum_{i!=j} 1/r_ij  ->  Rh = 2d for two beads d apart
    d = 3.0
    assert uc.hydrodynamic_radius_kirkwood([atom(1, 0.0), atom(2, d)]) == pytest.approx(2 * d)


def test_vdw_grid_volume_single_atom():
    r = uc.GAFF2_RADII['c3']
    exact = 4.0 / 3.0 * math.pi * r ** 3
    assert uc.vdw_volume_grid([atom(1, 0.0)], spacing=0.1) == pytest.approx(exact, rel=0.01)


def test_moment_of_inertia_diatomic():
    # two carbons 2 A apart: I_perp = 2 * m * 1^2, I_parallel = 0
    m = uc.MW_TABLE['C']
    mi = uc.moment_of_inertia([atom(1, -1.0), atom(2, 1.0)])
    assert mi['I1_amuA2'] == pytest.approx(0.0, abs=1e-6)
    assert mi['I3_amuA2'] == pytest.approx(2 * m, rel=1e-4)


def test_intrinsic_viscosity_units_agree():
    # legacy hard-sphere and Kirkwood paths must give the same [eta] in dL/g
    # for the same Rh (regression: legacy was 10x too small)
    rh, mw = 15.0, 5000.0
    legacy = uc.intrinsic_viscosity(uc.hydrodynamic_volume(rh), mw)
    kirkwood = uc.transport_from_Rh(rh, mw)['eta_iv_dLg']
    assert legacy == pytest.approx(kirkwood, rel=1e-3)


def test_terminal_element_tie_prefers_nitrogen():
    el, src = uc.terminal_element_from_monomer({'N': 1, 'O': 1, 'C': 1})
    assert el == 'N' and 'TIE' in src
    assert uc.terminal_element_from_monomer({}, override='o')[0] == 'O'


def test_infer_meta_from_filename():
    fam, gen, state = uc.infer_meta('/data/PAMAM/dendrimer_pamam_gen3_min.rst')
    assert (fam, gen, state) == ('PAMAM', 3, 'amber_min')


# ── worked example ─────────────────────────────────────────────────────────────
@pytest.fixture(scope='module')
def pamam_g2():
    atoms, bonds, note = uc.load_structure(STRUCT)
    _, _, core_info, _ = uc.analyse_fragment(CORE, 'core')
    _, _, mono_info, mono_el = uc.analyse_fragment(MONO, 'monomer')
    te, src = uc.terminal_element_from_monomer(mono_el)
    row, _ = uc.characterise(atoms, bonds, 'PAMAM', 2, 'amber_min', STRUCT, note,
                             te, src, core_info, mono_info, mc_samples=20000)
    return row


def test_pamam_g2_counts_and_qc(pamam_g2):
    r = pamam_g2
    assert (r['N_atoms'], r['N_bonds'], r['N_heavy']) == (516, 515, 228)
    assert r['n_term_detected'] == r['n_term_theoretical'] == 16
    assert r['QC_PASS'] and r['connected'] and r['is_tree']
    assert r['MW_Da'] == pytest.approx(3256.24, abs=0.01)


def test_pamam_g2_size_and_hydrodynamics(pamam_g2):
    r = pamam_g2
    assert r['Rg_heavy_A'] == pytest.approx(13.922, abs=1e-3)
    assert r['Rh_kirkwood_A'] == pytest.approx(11.198, abs=1e-3)
    assert r['V_vdW_grid_A3'] == pytest.approx(3865.9, abs=0.1)
    assert r['SASA_A2'] == pytest.approx(3463.2, abs=0.1)
    assert abs(r['net_charge_e']) < 1e-3


# ── command line ───────────────────────────────────────────────────────────────
def test_cli_writes_csv_and_units(tmp_path):
    out = tmp_path / 'res.csv'
    cmd = [sys.executable, os.path.join(ROOT, 'udb_characterizer.py'), STRUCT,
           '--monomer', MONO, '--terminal', MONO, '--family', 'PAMAM', '--gen', '2', '-o', str(out),
           '--mc-samples', '5000', '--sasa-points', '120', '--units-row']
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    assert res.returncode == 0, res.stderr
    rows = list(csv.DictReader(open(out, encoding='utf-8')))
    assert rows[0]['Rg_heavy_A'] == 'Å'          # units row
    assert rows[1]['N_atoms'] == '516'
    assert rows[1]['terminalfrag_n_atoms'] == '20'
    assert (tmp_path / 'res_units.csv').exists()
