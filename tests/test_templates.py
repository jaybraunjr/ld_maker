import numpy as np
import pytest

from ld_maker.templates import extract_templates


def test_extracts_one_per_resname(reference):
    tpls = extract_templates(reference, ["POPC", "TRIO", "CHYO"])
    assert set(tpls) == {"POPC", "TRIO", "CHYO"}
    assert tpls["POPC"].resname == "POPC"


def test_atom_count_and_centering(reference):
    tpl = extract_templates(reference, ["TRIO"])["TRIO"]
    assert tpl.n_atoms == 12
    # centred on its geometric centre
    assert np.allclose(tpl.positions.mean(axis=0), 0.0, atol=1e-4)


def test_lipid_oriented_head_up(reference):
    # reference POPC heads point up already; extraction should keep P on the +z side
    tpl = extract_templates(reference, ["POPC"], align_axis=True)["POPC"]
    assert tpl.head_up
    p_z = tpl.positions[tpl.names.index("P"), 2]
    assert p_z > tpl.positions[:, 2].mean()


def test_missing_residue_raises(reference):
    with pytest.raises(ValueError, match="no residue named"):
        extract_templates(reference, ["NOPE"])


def test_random_rotation_preserves_shape(reference):
    tpl = extract_templates(reference, ["TRIO"])["TRIO"]
    rng = np.random.default_rng(0)
    rotated = tpl.randomly_rotated(rng)
    # a rigid rotation preserves all pairwise distances
    def pdist(p):
        d = p[:, None, :] - p[None, :, :]
        return np.sqrt((d ** 2).sum(-1))
    assert np.allclose(pdist(rotated), pdist(tpl.positions), atol=1e-6)
