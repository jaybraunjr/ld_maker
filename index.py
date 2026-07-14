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


def group_definitions(universe, membrane_extra=()):
    """Return an ordered dict of group-name -> MDAnalysis selection string,
    including only residue kinds actually present.

    ``membrane_extra`` forces resnames into MEMB (and out of the core SOLU
    group) -- use it for surface cholesteryl ester, where CHYO is a membrane
    lipid rather than a core lipid.
    """
    present = set(universe.residues.resnames)
    extra = set(membrane_extra)
    core = [r for r in CORE if r in present and r not in extra]
    memb = [r for r in PHOSPHOLIPIDS if r in present] + \
           [r for r in extra if r in present]
    solv = [r for r in SOLVENT if r in present]
    sel = lambda names: "resname " + " ".join(names) if names else "resid -1"
    return {
        "SOLU": sel(core),
        "MEMB": sel(memb),
        "SOLV": sel(solv),
        "SOLU_MEMB": sel(core + memb),
        "SYSTEM": "all",
    }


def write_index_ndx(structure, path, membrane_extra=()):
    """Write ``path`` (index.ndx) for ``structure`` (a .gro path or Universe).

    ``membrane_extra`` forces resnames (e.g. surface ``CHYO``) into MEMB.
    """
    u = structure if isinstance(structure, mda.Universe) else mda.Universe(structure)
    groups = group_definitions(u, membrane_extra=membrane_extra)
    with open(path, "w") as fh:
        for name, selection in groups.items():
            idx = u.select_atoms(selection).indices + 1        # 1-based for GROMACS
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(map(str, idx[i:i + 15])) + "\n")
    return {name: int(len(u.select_atoms(s))) for name, s in groups.items()}


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    membrane = ()
    for a in sys.argv[1:]:
        if a.startswith("--membrane="):
            membrane = tuple(a.split("=", 1)[1].split(","))
    if len(args) != 2:
        print("usage: python -m ld_maker.index <structure.gro> <index.ndx> "
              "[--membrane=CHYO,...]")
        raise SystemExit(1)
    counts = write_index_ndx(args[0], args[1], membrane_extra=membrane)
    for name, n in counts.items():
        print(f"  {name:10s} {n}")
