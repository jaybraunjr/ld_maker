import shutil

import numpy as np
import pytest

from ld_maker.bilayer import insert_core_into_bilayer, replace_lipid, _leaflet_split_z


def test_leaflet_split_z_is_membrane_centre(bilayer):
    # top heads at z=60, bottom at z=20 -> midplane ~40
    assert _leaflet_split_z(bilayer) == pytest.approx(40.0, abs=1.0)


def test_grid_core_grows_box_and_adds_trio(reference, bilayer):
    z0 = bilayer.dimensions[2]
    ld = insert_core_into_bilayer(
        bilayer, {"TRIO": 6}, core_source=reference, method="grid",
        packing_fraction=0.2, relax=False, seed=0,
    )
    comp = {}
    for r in ld.residues:
        comp[r.resname] = comp.get(r.resname, 0) + 1
    assert comp["TRIO"] == 6
    assert comp["POPC"] == 8 and comp["DOPE"] == 4      # leaflets unchanged
    assert ld.dimensions[2] > z0                         # box grew in z for the core


def test_replace_lipid_composition(reference, bilayer):
    ld = replace_lipid(bilayer, {"POPC": ("CHYO", 0.5)}, source=reference, seed=0)
    comp = {}
    for r in ld.residues:
        comp[r.resname] = comp.get(r.resname, 0) + 1
    # 8 POPC total, half -> 4 CHYO placed, 4 POPC remain; DOPE untouched
    assert comp["POPC"] == 4
    assert comp["CHYO"] == 4
    assert comp["DOPE"] == 4


def test_int_core_is_trio_shorthand(reference, bilayer):
    ld = insert_core_into_bilayer(bilayer, 5, core_source=reference,
                                  method="grid", relax=False)
    assert sum(r.resname == "TRIO" for r in ld.residues) == 5


def test_replace_lipid_rejects_missing_lipid(reference, bilayer):
    with pytest.raises(ValueError, match="no residue named"):
        replace_lipid(bilayer, {"NOPE": ("CHYO", 0.5)}, source=reference)


def test_replace_lipid_rejects_bad_fraction(reference, bilayer):
    with pytest.raises(ValueError, match="fraction"):
        replace_lipid(bilayer, {"POPC": ("CHYO", 1.5)}, source=reference)


@pytest.mark.skipif(shutil.which("gmx") is None, reason="needs GROMACS gmx on PATH")
def test_insert_method_runs_with_gmx(reference, bilayer):
    ld = insert_core_into_bilayer(bilayer, {"TRIO": 4, "CHYO": 2},
                                  core_source=reference, method="insert", seed=0)
    comp = {r.resname for r in ld.residues}
    assert {"TRIO", "CHYO"} <= comp
