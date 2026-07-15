from ld_maker.vmd import write_vmd


def _write_gro(universe, tmp_path, name="sys.gro"):
    path = tmp_path / name
    universe.atoms.write(str(path))
    return str(path)


def test_vmd_has_reps_for_present_components(reference, tmp_path):
    gro = _write_gro(reference, tmp_path)
    out = tmp_path / "view.vmd"
    write_vmd(gro, str(out))
    text = out.read_text()
    assert "mol new" in text and "mol addrep top" in text
    for sel in ("resname TRIO", "resname CHYO", "resname POPC", "name P"):
        assert "{" + sel + "}" in text            # a representation for each present kind


def test_vmd_omits_absent_components(bilayer, tmp_path):
    gro = _write_gro(bilayer, tmp_path, "bil.gro")
    out = tmp_path / "bil.vmd"
    write_vmd(gro, str(out))
    text = out.read_text()
    assert "resname TRIO" not in text             # bilayer has no core
    assert "POPE" not in text                      # the "other lipids" rep is skipped
    assert "{resname POPC}" in text


def test_vmd_trajectory_line(reference, tmp_path):
    gro = _write_gro(reference, tmp_path)
    out = tmp_path / "traj.vmd"
    write_vmd(gro, str(out), trajectory="run.xtc")
    assert 'mol addfile "run.xtc"' in out.read_text()


def test_vmd_in_memory_universe_rejected(reference, tmp_path):
    import pytest
    with pytest.raises(ValueError, match="file path"):
        write_vmd(reference, str(tmp_path / "x.vmd"))
