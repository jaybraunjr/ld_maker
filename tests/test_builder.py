import numpy as np

from ld_maker.builder import build_from_reference, _leaflet_composition_counts


def test_composition_counts_sum_exactly():
    counts = _leaflet_composition_counts(100, {"POPC": 176, "DOPE": 74, "SAPI": 20})
    assert sum(counts.values()) == 100
    assert counts["POPC"] > counts["DOPE"] > counts["SAPI"]


def test_grid_build_counts_and_box(reference):
    u = build_from_reference(
        reference, nx=3, ny=3, n_trio=5,
        pl_composition={"POPC": 2, "DOPE": 1},
        core_thickness=20.0, water_thickness=10.0, seed=0,
    )
    comp = {}
    for r in u.residues:
        comp[r.resname] = comp.get(r.resname, 0) + 1
    assert comp["TRIO"] == 5
    assert comp["POPC"] + comp["DOPE"] == 3 * 3 * 2      # two leaflets, nx*ny each
    assert all(d > 0 for d in u.dimensions[:3])


def test_grid_build_reproducible(reference):
    kw = dict(nx=3, ny=3, n_trio=4, pl_composition={"POPC": 1},
              core_thickness=20.0, water_thickness=10.0, seed=7)
    a = build_from_reference(reference, **kw)
    b = build_from_reference(reference, **kw)
    assert np.allclose(a.atoms.positions, b.atoms.positions)
