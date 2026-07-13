"""Generate the GROMACS ``[ molecules ]`` section for a built system.

Mirrors the counting logic in ``analysis/make_big.ipynb``: walk residues in file
order, collapse consecutive runs of the same residue name, and emit one line per
run.  A full .top still needs the matching ``#include`` lines for your force
field / lipid itp files (see :func:`write_top`).
"""

import MDAnalysis as mda


def molecule_counts(universe):
    """Return an ordered list of ``(resname, count)`` runs, in file order."""
    runs = []
    for res in universe.residues:
        rn = res.resname
        if runs and runs[-1][0] == rn:
            runs[-1][1] += 1
        else:
            runs.append([rn, 1])
    return [(rn, n) for rn, n in runs]


def molecules_section(universe):
    """Return the text of a GROMACS ``[ molecules ]`` block."""
    lines = ["[ molecules ]", "; Compound        #mols"]
    for rn, n in molecule_counts(universe):
        lines.append(f"{rn:<16} {n}")
    return "\n".join(lines) + "\n"


def write_top(universe, path, includes=None, system_name="Planar LD"):
    """Write a minimal but complete .top file.

    Parameters
    ----------
    includes : list[str] or None
        Force-field / itp paths to ``#include``.  Defaults to the CHARMM36 set
        used in this project's ``toppar`` directory.
    """
    if includes is None:
        includes = [
            "toppar/forcefield.itp",
            "toppar/POPC.itp",
            "toppar/DOPE.itp",
            "toppar/SAPI.itp",
            "toppar/TIP3.itp",
            "toppar/SOD.itp",
            "toppar/CLA.itp",
        ]
    with open(path, "w") as fh:
        for inc in includes:
            fh.write(f'#include "{inc}"\n')
        fh.write(f"\n[ system ]\n{system_name}\n\n")
        fh.write(molecules_section(universe))
    return path
