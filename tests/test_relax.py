import numpy as np

from ld_maker.relax import relax_clashes


def _min_inter_distance(mols, box):
    from MDAnalysis.lib.distances import self_capped_distance
    coords = np.vstack([m[2] for m in mols]).astype(np.float32)
    sizes = [len(m[1]) for m in mols]
    mol_of = np.repeat(np.arange(len(mols)), sizes)
    pairs, d = self_capped_distance(coords, max_cutoff=5.0, box=np.asarray(box, np.float32),
                                    return_distances=True)
    inter = [dd for (a, b), dd in zip(pairs, d) if mol_of[a] != mol_of[b]]
    return min(inter) if inter else 9.0


def test_separates_overlapping_molecules():
    # two identical atom blobs sitting almost on top of each other
    names = [f"A{i}" for i in range(6)]
    base = np.random.default_rng(0).uniform(-2, 2, size=(6, 3))
    mols = [("LIG", names, base.copy()),
            ("LIG", names, base.copy() + 0.3)]
    box = [40.0, 40.0, 40.0, 90.0, 90.0, 90.0]
    assert _min_inter_distance(mols, box) < 1.0
    relax_clashes(mols, box=box, d_min=2.0, iterations=200)
    assert _min_inter_distance(mols, box) > 1.0


def test_lateral_only_keeps_z_fixed():
    names = [f"A{i}" for i in range(6)]
    base = np.zeros((6, 3))
    base[:, 2] = np.arange(6)
    mols = [("PL", names, base.copy()), ("PL", names, base.copy() + [0.2, 0, 0])]
    z_before = [m[2][:, 2].copy() for m in mols]
    relax_clashes(mols, box=[30.0, 30.0, 30.0, 90.0, 90.0, 90.0], d_min=2.0,
                  lateral_only_resnames=("PL",))
    for m, z0 in zip(mols, z_before):
        assert np.allclose(m[2][:, 2], z0)   # z untouched for lateral-only residues
