"""End-to-end: each documented build case produces a valid, viewable system.

Every case builds the structure, writes .gro/.top/.vmd, and checks the
composition, a well-formed topology, and a VMD script with the right reps.
"""

import shutil

import pytest

from ld_maker import (build_from_reference, insert_core_into_bilayer,
                      replace_lipid, write_top, write_vmd)


def _composition(u):
    c = {}
    for r in u.residues:
        c[r.resname] = c.get(r.resname, 0) + 1
    return c


def _write_all(u, tmp_path, name):
    gro, top, vmd = (tmp_path / f"{name}.gro", tmp_path / f"{name}.top",
                     tmp_path / f"{name}.vmd")
    u.atoms.write(str(gro))
    write_top(u, str(top))
    write_vmd(str(gro), str(vmd))
    assert "[ molecules ]" in top.read_text()
    assert "mol addrep top" in vmd.read_text()
    return gro, top, vmd


def test_case_grid_build(reference, tmp_path):
    u = build_from_reference(reference, nx=3, ny=3, n_trio=6,
                             pl_composition={"POPC": 2, "DOPE": 1},
                             core_thickness=20.0, water_thickness=10.0)
    comp = _composition(u)
    assert comp["TRIO"] == 6
    assert comp["POPC"] + comp["DOPE"] == 3 * 3 * 2
    _, _, vmd = _write_all(u, tmp_path, "grid")
    assert "{resname TRIO}" in vmd.read_text()


def test_case_bilayer_pure_trio_core(reference, bilayer, tmp_path):
    u = insert_core_into_bilayer(bilayer, {"TRIO": 6}, core_source=reference,
                                 method="grid", packing_fraction=0.2, relax=False)
    assert _composition(u)["TRIO"] == 6
    _, _, vmd = _write_all(u, tmp_path, "bilcore")
    text = vmd.read_text()
    assert "{resname TRIO}" in text and "{resname POPC}" in text


def test_case_replace_surface_chyo(reference, bilayer, tmp_path):
    u = replace_lipid(bilayer, {"POPC": ("CHYO", 0.5)}, source=reference)
    comp = _composition(u)
    assert comp["CHYO"] == 4 and comp["POPC"] == 4
    _, _, vmd = _write_all(u, tmp_path, "surf")
    assert "{resname CHYO}" in vmd.read_text()


@pytest.mark.skipif(shutil.which("gmx") is None, reason="needs GROMACS gmx on PATH")
def test_case_mixed_core_insert(reference, bilayer, tmp_path):
    u = insert_core_into_bilayer(bilayer, {"TRIO": 4, "CHYO": 2},
                                 core_source=reference, method="insert")
    assert {"TRIO", "CHYO"} <= set(_composition(u))
    _write_all(u, tmp_path, "mixed")
