"""Tests for dendrimer_param_gui_v2.py (MOL2 typing, conversion and validation)."""
import os
import time

import pytest

pytest.importorskip('tkinter')
import dendrimer_param_gui_v2 as pg  # noqa: E402
from conftest import EXAMPLE  # noqa: E402

SYBYL_MONO = os.path.join(EXAMPLE, 'PAMAM_monomer_sybyl.mol2')
GAFF2_G2 = os.path.join(EXAMPLE, 'dendrimer_pamam_gen2_min_gaff.mol2')


def mol2_types(path):
    types, sec = [], None
    for line in open(path, encoding='utf-8'):
        if line.startswith('@<TRIPOS>'):
            sec = line.strip()[9:]
        elif sec == 'ATOM' and line.split():
            types.append(line.split()[5])
    return types


def test_parse_sybyl_monomer():
    atoms, bonds, adj, _, _ = pg.parse_mol2(SYBYL_MONO)
    assert (len(atoms), len(bonds)) == (20, 19)
    assert dict(pg.detect_elements(atoms)) == {'C': 5, 'H': 12, 'N': 2, 'O': 1}
    assert not pg.needs_bypass(pg.detect_elements(atoms))


def test_sybyl_to_gaff2_typing_and_validation(tmp_path):
    atoms, bonds, adj, ms, bs = pg.parse_mol2(SYBYL_MONO)
    out = tmp_path / 'mono_typed.mol2'
    pg.write_typed_mol2(atoms, bonds, adj, ms, bs, 'PAMAM_monomer', str(out))
    types = set(mol2_types(out))
    assert types <= {'c3', 'ca', 'n3', 'o', 'os', 'hc', 'hn'}
    assert 'DU' not in types
    results = pg.validate_outputs(atoms, bonds, str(out), str(tmp_path / 'none.frcmod'), 0)
    failed = [msg for icon, msg, _ in results if icon == '❌']
    assert failed == ['frcmod not found']        # everything else passes without AmberTools


def test_gaff2_to_sybyl_map():
    m = pg.GAFF2_TO_SYBYL
    assert (m['c3'], m['c'], m['n'], m['n3'], m['o'], m['os'], m['hn']) == \
           ('C.3', 'C.2', 'N.am', 'N.3', 'O.2', 'O.3', 'H')


@pytest.mark.skipif(not os.environ.get('DISPLAY'), reason='needs a display (run under xvfb-run)')
def test_extract_tab_converts_gaff2_to_sybyl(tmp_path):
    # regression: the Extract mol2 worker thread was never started
    app = pg.ParamApp()
    app.ext_mode.set('A'); app.ext_typed.set(GAFF2_G2); app.ext_dir.set(str(tmp_path))
    app.ext_name.set('g2'); app.ext_pdb.set(False); app.ext_gaff2.set(False)
    deadline = time.time() + 60
    log = []

    def poll():
        text = app.ext_log.get('1.0', 'end')
        if 'Conversion complete' in text or 'Error' in text or time.time() > deadline:
            log.append(text); app.quit()
        else:
            app.after(100, poll)

    # the worker posts back with after(), so Tk's mainloop must be running
    app.after(0, app._run_extract)
    app.after(100, poll)
    app.mainloop()
    app.destroy()
    assert 'Conversion complete' in log[0], log[0]
    types = mol2_types(tmp_path / 'g2_sybyl.mol2')
    assert len(types) == 516
    assert set(types) == {'C.2', 'C.3', 'H', 'N.3', 'N.am', 'O.2'}
