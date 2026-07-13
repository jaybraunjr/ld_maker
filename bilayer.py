"""Convert a phospholipid bilayer into a lipid droplet by inserting a TRIO core.

This is the physically-motivated route (as opposed to the from-scratch grid
builder in :mod:`ld_maker.builder`): start from an equilibrated, solvated,
ionized CHARMM-GUI bilayer, pry the two leaflets apart along z, and fill the
vacated midplane with neutral lipid (TRIO).  The bilayer's own water and ions are
kept exactly where they are (they sit outside the membrane); only the hydrophobic
interior is touched.  The result is the LD "trilayer": PL monolayer / TRIO core /
PL monolayer, ready for the standard minimization + equilibration pipeline.

    bilayer                          lipid droplet
  ~~~~~~~~~~~~  heads              ~~~~~~~~~~~~   monolayer
  ||||||||||||  tails      -->     ||||||||||||
  ||||||||||||  tails              ############   <- inserted TRIO core
  ~~~~~~~~~~~~  heads              ||||||||||||
                                   ~~~~~~~~~~~~   monolayer
"""

import numpy as np
import MDAnalysis as mda

from .templates import extract_templates
from .relax import relax_clashes

# approximate molecular volume of triolein (g/mol 885, rho ~0.91 g/cm^3)
TRIO_VOLUME_A3 = 1630.0


def _leaflet_split_z(universe, head_name="P"):
    """Return the membrane mid-plane z (mean of all phosphate atoms)."""
    heads = universe.select_atoms(f"name {head_name}")
    if not len(heads):
        raise ValueError(f"no '{head_name}' atoms found - is this a phospholipid bilayer?")
    return float(heads.positions[:, 2].mean())


def insert_core_into_bilayer(
    bilayer,
    n_trio,
    trio_source="modules/run.gro",
    core_thickness=None,
    packing_fraction=0.35,
    head_name="P",
    seed=0,
    relax=True,
    d_min=1.6,
):
    """Insert ``n_trio`` TRIO molecules into a bilayer's midplane.

    Parameters
    ----------
    bilayer : str or mda.Universe
        A solvated phospholipid bilayer (e.g. a CHARMM-GUI ``.gro``).
    n_trio : int
        Number of TRIO molecules to insert (the oil core).
    trio_source : str or mda.Universe
        Structure to take a single TRIO template from (must contain TRIO).
    core_thickness : float or None
        Thickness of the inserted core slab in angstrom.  If ``None`` it is
        derived from ``n_trio`` and the box area via ``packing_fraction``.
    packing_fraction : float
        Target initial fill of the core (0-1).  Lower = looser start = safer
        minimization; NPT then compresses the box to the real density.
    head_name : str
        Atom name marking phospholipid head groups (``P`` for PC/PE/PI/PS/PG).
    relax : bool
        Optional rigid-body clash relaxation after insertion.  Usually
        unnecessary — the loose core + GROMACS minimization handle contacts.

    Returns
    -------
    mda.Universe
        The LD system (bilayer + TRIO core), box extended along z.
    """
    u = bilayer if isinstance(bilayer, mda.Universe) else mda.Universe(bilayer)
    # make molecules whole and inside the box so leaflet/water assignment is clean
    u.atoms.wrap(compound="residues")
    box = u.dimensions.copy()
    area = box[0] * box[1]

    center = _leaflet_split_z(u, head_name)

    if core_thickness is None:
        core_thickness = n_trio * TRIO_VOLUME_A3 / (area * packing_fraction)
    t = float(core_thickness)

    # --- assign every residue to the top or bottom half by its COM z, then
    #     translate the two halves apart to open a gap at the midplane ---
    res_com_z = np.array([r.atoms.center_of_mass()[2] for r in u.residues])
    top_mask = res_com_z >= center
    top_atoms = u.residues[top_mask].atoms
    bot_atoms = u.residues[~top_mask].atoms
    top_atoms.positions += np.array([0.0, 0.0, +t / 2.0])
    bot_atoms.positions += np.array([0.0, 0.0, -t / 2.0])

    # --- place TRIO on a jittered 3-D grid in the vacated gap ---
    trio = extract_templates(trio_source, ["TRIO"], align_axis=True)["TRIO"]
    rng = np.random.default_rng(seed)
    span = trio.positions.max(axis=0) - trio.positions.min(axis=0)
    cell = 0.85 * float(max(span[0], span[1]))
    nx = max(1, int(box[0] // cell))
    ny = max(1, int(box[1] // cell))
    nz = max(1, int(np.ceil(n_trio / (nx * ny))))
    sx, sy = box[0] / nx, box[1] / ny
    inset = min(0.30 * t, 0.25 * float(span[2]))
    z_lo, z_hi = center - t / 2 + inset, center + t / 2 - inset
    sz = (z_hi - z_lo) / nz if nz else 0.0
    jit = 0.2 * min(sx, sy, max(sz, 1.0))

    slots = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    rng.shuffle(slots)
    trio_mols = []
    for (i, j, k) in slots[:n_trio]:
        coords = trio.randomly_rotated(rng)
        centre = np.array([
            (i + 0.5) * sx + rng.uniform(-jit, jit),
            (j + 0.5) * sy + rng.uniform(-jit, jit),
            z_lo + (k + 0.5) * sz + rng.uniform(-jit, jit),
        ])
        trio_mols.append(coords + centre)

    ld = _merge_bilayer_and_core(u, trio, trio_mols, box, t)

    # recenter so the core (originally at `center`) sits at the new box midplane,
    # keeping the whole system inside [0, box_z] with symmetric water slabs
    ld.atoms.translate([0.0, 0.0, ld.dimensions[2] / 2.0 - center])

    if relax:
        # relax only the core + lipids (skip bulk water), then write coords back
        mols, idx = _core_and_lipid_molecules(ld, head_name)
        relax_clashes(mols, box=ld.dimensions, d_min=d_min,
                      lateral_only_resnames=_phospholipid_resnames(ld, head_name))
        relaxed = np.vstack([m[2] for m in mols])
        ld.atoms[np.concatenate(idx)].positions = relaxed
    return ld


def _merge_bilayer_and_core(u, trio_tpl, trio_mols, box, t):
    """Combine the shifted bilayer with the TRIO molecules into one Universe."""
    n_trio = len(trio_mols)
    trio_natoms = trio_tpl.n_atoms
    total = u.atoms.n_atoms + n_trio * trio_natoms

    n_res = u.residues.n_residues + n_trio
    # residue atom-count array: existing residues then TRIO residues
    exist_counts = [r.atoms.n_atoms for r in u.residues]
    counts = exist_counts + [trio_natoms] * n_trio
    atom_resindex = np.repeat(np.arange(n_res), counts)

    new = mda.Universe.empty(n_atoms=total, n_residues=n_res,
                             atom_resindex=atom_resindex, trajectory=True)
    new.add_TopologyAttr("names",
                         list(u.atoms.names) + list(trio_tpl.names) * n_trio)
    new.add_TopologyAttr("resnames",
                         list(u.residues.resnames) + ["TRIO"] * n_trio)
    new.add_TopologyAttr("resids", np.arange(1, n_res + 1))

    coords = np.vstack([u.atoms.positions] + trio_mols)
    new.atoms.positions = coords
    box = box.copy()
    box[2] += t                     # grow z to accommodate the core
    new.dimensions = box
    return new


def _phospholipid_resnames(universe, head_name="P"):
    """Resnames that carry a head-group atom (the phospholipids)."""
    names = set()
    for r in universe.residues:
        if head_name in r.atoms.names:
            names.add(r.resname)
    return tuple(names)


def _core_and_lipid_molecules(universe, head_name="P"):
    """Return (molecules, indices) for lipids + TRIO, for relaxation + write-back.

    ``molecules`` is a list of ``(resname, names, coords)`` tuples (coords are
    copies that relaxation mutates); ``indices`` are the matching atom indices in
    ``universe`` so the relaxed coordinates can be written back.
    """
    keep = _phospholipid_resnames(universe, head_name) + ("TRIO",)
    mols, idx = [], []
    for r in universe.residues:
        if r.resname in keep:
            mols.append((r.resname, list(r.atoms.names), r.atoms.positions.copy()))
            idx.append(r.atoms.indices)
    return mols, idx


def bilayer_to_ld(bilayer, n_trio, out_gro=None, out_top=None,
                  trio_source="modules/run.gro", toppar_includes=None, **kwargs):
    """Build an LD from a bilayer and (optionally) write .gro + .top.

    Extra keyword arguments are passed to :func:`insert_core_into_bilayer`
    (``core_thickness``, ``packing_fraction``, ``seed``, ``relax`` ...).
    Returns the resulting Universe.
    """
    ld = insert_core_into_bilayer(bilayer, n_trio, trio_source=trio_source, **kwargs)
    if out_gro:
        ld.atoms.write(out_gro)
    if out_top:
        _write_top(ld, out_top, toppar_includes)
    return ld


def _write_top(universe, path, includes=None, system_name="Lipid droplet (bilayer + TRIO)"):
    """Write topol.top with [ molecules ] counted in file order (make_big method)."""
    from .topology import molecules_section
    if includes is None:
        # include forcefield + one itp per residue kind present
        order = ["POPC", "DOPE", "SAPI", "POPE", "POPG", "POPS", "PSM", "CHL1",
                 "TRIO", "CHYO", "SOD", "CLA", "POT", "TIP3"]
        present = [r for r in order if r in set(universe.residues.resnames)]
        includes = ["toppar/charmm36-mar2019.ff/forcefield.itp"] + \
                   [f"toppar/{r}.itp" for r in present]
    with open(path, "w") as fh:
        for inc in includes:
            fh.write(f'#include "{inc}"\n')
        fh.write(f"\n[ system ]\n{system_name}\n\n")
        fh.write(molecules_section(universe))
    return path
