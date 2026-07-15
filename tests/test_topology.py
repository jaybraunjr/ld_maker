from ld_maker.topology import molecule_counts, molecules_section, write_top


def test_counts_group_consecutive_runs(reference):
    counts = molecule_counts(reference)
    # file order in the reference fixture: POPC, DOPE, SAPI, TRIO, CHYO
    assert counts == [("POPC", 3), ("DOPE", 2), ("SAPI", 2), ("TRIO", 6), ("CHYO", 4)]


def test_counts_split_when_not_contiguous(noncontig_universe):
    # a resname that reappears after another must produce two runs (GROMACS order)
    assert molecule_counts(noncontig_universe) == [("TRIO", 1), ("POPC", 1), ("TRIO", 1)]


def test_molecules_section_text(reference):
    text = molecules_section(reference)
    assert text.startswith("[ molecules ]")
    assert "POPC" in text and "CHYO" in text


def test_write_top(reference, tmp_path):
    out = tmp_path / "sys.top"
    write_top(reference, str(out), includes=["toppar/forcefield.itp"])
    content = out.read_text()
    assert '#include "toppar/forcefield.itp"' in content
    assert "[ system ]" in content
    assert "[ molecules ]" in content
