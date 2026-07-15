"""Shared pytest fixtures: tiny in-memory structures so tests need no big files.

The maker's logic is geometric/bookkeeping, so synthetic residues (a phosphate
head plus a short carbon tail for lipids, an atom blob for neutral lipids) are
enough to exercise template extraction, grid building, topology, index groups and
the bilayer geometry without shipping a real .gro.
"""

import os
import sys

import numpy as np
import MDAnalysis as mda
import pytest

# make the package importable when running pytest from the repo dir without an
# install (the package lives at the repo root, so its parent must be on the path)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def make_universe(residues, box=(120.0, 120.0, 120.0)):
    """Build an in-memory Universe from ``[(resname, names, positions), ...]``."""
    counts = [len(names) for _, names, _ in residues]
    n_res = len(residues)
    u = mda.Universe.empty(
        n_atoms=sum(counts), n_residues=n_res,
        atom_resindex=np.repeat(np.arange(n_res), counts), trajectory=True,
    )
    u.add_TopologyAttr("names", [nm for _, names, _ in residues for nm in names])
    u.add_TopologyAttr("resnames", [rn for rn, _, _ in residues])
    u.add_TopologyAttr("resids", np.arange(1, n_res + 1))
    u.add_TopologyAttr("masses", np.full(sum(counts), 12.0))   # COM/wrap need masses
    u.atoms.positions = np.vstack([p for _, _, p in residues])
    u.dimensions = [*box, 90.0, 90.0, 90.0]
    return u


def lipid(resname, x, y, z_head, points_down, n_tail=6):
    """A linear lipid: P head at ``z_head``, tail extending away from it."""
    names = ["P"] + [f"C{i}" for i in range(n_tail)]
    step = -1.0 if points_down else 1.0
    z = z_head + step * np.arange(len(names), dtype=float)
    pos = np.column_stack([np.full(len(names), x), np.full(len(names), y), z])
    return (resname, names, pos)


def blob(resname, center, n=10, seed=0):
    rng = np.random.default_rng(seed)
    names = [f"A{i}" for i in range(n)]
    return (resname, names, np.asarray(center) + rng.uniform(-4, 4, size=(n, 3)))


def water(x, y, z):
    return ("TIP3", ["OH2", "H1", "H2"],
            np.array([[x, y, z], [x + 1, y, z], [x, y + 1, z]], dtype=float))


@pytest.fixture
def reference():
    """A reference structure with >=1 of each residue the maker uses."""
    res, s = [], 0
    for i in range(3):
        res.append(lipid("POPC", i * 8.0, 0.0, 50.0, points_down=False))
    for i in range(2):
        res.append(lipid("DOPE", i * 8.0, 12.0, 50.0, points_down=False))
    for i in range(2):
        res.append(lipid("SAPI", i * 8.0, 24.0, 50.0, points_down=False))
    for i in range(6):
        res.append(blob("TRIO", [i * 8.0, 36.0, 50.0], n=12, seed=(s := s + 1)))
    for i in range(4):
        res.append(blob("CHYO", [i * 8.0, 48.0, 50.0], n=10, seed=(s := s + 1)))
    return make_universe(res)


@pytest.fixture
def bilayer():
    """A tiny solvated bilayer: two POPC/DOPE leaflets, water above and below."""
    res = []
    # top leaflet: heads high (~60), tails pointing down toward the midplane
    for i in range(4):
        res.append(lipid("POPC", i * 12.0, 0.0, 60.0, points_down=True))
    for i in range(2):
        res.append(lipid("DOPE", i * 12.0, 12.0, 60.0, points_down=True))
    # bottom leaflet: heads low (~20), tails pointing up
    for i in range(4):
        res.append(lipid("POPC", i * 12.0, 0.0, 20.0, points_down=False))
    for i in range(2):
        res.append(lipid("DOPE", i * 12.0, 12.0, 20.0, points_down=False))
    # water slabs outside the membrane
    for i in range(6):
        res.append(water(i * 8.0, 30.0, 75.0))
        res.append(water(i * 8.0, 30.0, 5.0))
    return make_universe(res, box=(48.0, 48.0, 90.0))


@pytest.fixture
def noncontig_universe():
    """A universe where TRIO reappears after POPC (two separate [molecules] runs)."""
    return make_universe([blob("TRIO", [0, 0, 0]), blob("POPC", [9, 0, 0]),
                          blob("TRIO", [18, 0, 0])])
