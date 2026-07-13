"""Generate a GROMACS index.ndx with the thermostat/restraint groups the
CHARMM-GUI equilibration protocol expects (``SOLU``, ``MEMB``, ``SOLV`` ...).

For a lipid-droplet system the groups map to the physical layers:
    SOLU  = TRIO neutral-lipid core
    MEMB  = phospholipids (POPC/DOPE/SAPI/...)
    SOLV  = water + ions
    SOLU_MEMB = SOLU + MEMB   (comm_grps)
    SYSTEM    = everything
This keeps all three thermostat groups non-empty even without a protein.
"""

import sys
import numpy as np
import MDAnalysis as mda

# residue-name buckets (extend as needed for new lipids/sterols)
CORE = ["TRIO", "CHYO"]
PHOSPHOLIPIDS = ["POPC", "DOPE", "SAPI", "POPE", "POPG", "POPS", "PSM", "CHL1"]
SOLVENT = ["TIP3", "SOL", "SOD", "CLA", "POT", "K", "NA"]


def group_definitions(universe):
    """Return an ordered dict of group-name -> MDAnalysis selection string,
    including only residue kinds actually present."""
    present = set(universe.residues.resnames)
    core = [r for r in CORE if r in present]
    memb = [r for r in PHOSPHOLIPIDS if r in present]
    solv = [r for r in SOLVENT if r in present]
    sel = lambda names: "resname " + " ".join(names) if names else "resid -1"
    return {
        "SOLU": sel(core),
        "MEMB": sel(memb),
        "SOLV": sel(solv),
        "SOLU_MEMB": sel(core + memb),
        "SYSTEM": "all",
    }


def write_index_ndx(structure, path):
    """Write ``path`` (index.ndx) for ``structure`` (a .gro path or Universe)."""
    u = structure if isinstance(structure, mda.Universe) else mda.Universe(structure)
    groups = group_definitions(u)
    with open(path, "w") as fh:
        for name, selection in groups.items():
            idx = u.select_atoms(selection).indices + 1        # 1-based for GROMACS
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(map(str, idx[i:i + 15])) + "\n")
    return {name: int(len(u.select_atoms(s))) for name, s in groups.items()}


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python -m ld_maker.index <structure.gro> <index.ndx>")
        raise SystemExit(1)
    counts = write_index_ndx(sys.argv[1], sys.argv[2])
    for name, n in counts.items():
        print(f"  {name:10s} {n}")
