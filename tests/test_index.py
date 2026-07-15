from ld_maker.index import group_definitions, write_index_ndx


def test_core_and_membrane_split(reference):
    groups = group_definitions(reference)
    assert "TRIO" in groups["SOLU"] and "CHYO" in groups["SOLU"]   # neutral core
    assert "POPC" in groups["MEMB"] and "DOPE" in groups["MEMB"]   # phospholipids
    assert groups["SYSTEM"] == "all"


def test_membrane_extra_moves_chyo(reference):
    groups = group_definitions(reference, membrane_extra=("CHYO",))
    assert "CHYO" not in groups["SOLU"]     # no longer core
    assert "CHYO" in groups["MEMB"]         # now a membrane lipid


def test_write_index_ndx_counts_and_format(reference, tmp_path):
    out = tmp_path / "index.ndx"
    counts = write_index_ndx(reference, str(out))
    text = out.read_text()
    assert "[ SOLU ]" in text and "[ MEMB ]" in text and "[ SYSTEM ]" in text
    assert counts["SYSTEM"] == reference.atoms.n_atoms
    # indices are 1-based for GROMACS
    first_idx = int(text.split("[ SYSTEM ]")[1].split()[0])
    assert first_idx == 1
