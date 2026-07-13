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

import os

import numpy as np
import MDAnalysis as mda

from .templates import extract_templates
from .relax import relax_clashes

# packaged single-molecule core templates (one TRIO + one CHYO), so mixed
# cores work without pointing at an external structure
DEFAULT_CORE_SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "templates", "core_templates.gro")

# approximate molecular volume of triolein (g/mol 885, rho ~0.91 g/cm^3)
TRIO_VOLUME_A3 = 1630.0

# approximate molecular volumes (A^3) of neutral-lipid core species,
# used only to size the initial core gap (NPT then sets the real density)
CORE_VOLUME_A3 = {"TRIO": 1630.0, "CHYO": 1180.0}

# resname -> topology include filename, where it differs from "<RESNAME>.itp"
ITP_FILENAME = {"CHYO": "chyo.itp"}


def _leaflet_split_z(universe, head_name="P"):
    """Return the membrane mid-plane z (mean of all phosphate atoms)."""
    heads = universe.select_atoms(f"name {head_name}")
    if not len(heads):
        raise ValueError(f"no '{head_name}' atoms found - is this a phospholipid bilayer?")
    return float(heads.positions[:, 2].mean())


def insert_core_into_bilayer(
    bilayer,
    core,
    core_source=None,
    core_thickness=None,
    packing_fraction=0.35,
    head_name="P",
    seed=0,
    relax=True,
    d_min=1.6,
):
    """Insert a neutral-lipid core between a bilayer's leaflets.

    Parameters
    ----------
    bilayer : str or mda.Universe
        A solvated phospholipid bilayer (e.g. a CHARMM-GUI ``.gro``).
    core : int or dict
        The core composition.  An int is shorthand for ``{"TRIO": int}``.
        A dict mixes species, e.g. ``{"TRIO": 120, "CHYO": 30}`` for a
        triolein / cholesteryl-ester droplet.
    core_source : str or mda.Universe
        Structure to take single-molecule templates from; must contain every
        species named in ``core`` (e.g. a 50:50 TRIO/CHYO system for CHYO).
    core_thickness : float or None
        Thickness of the inserted core slab in angstrom.  If ``None`` it is
        derived from the composition and box area via ``packing_fraction``.
    packing_fraction : float
        Target initial fill of the core (0-1).  Lower = looser start = safer
        minimization; NPT then compresses the box to the real density.
    head_name : str
        Atom name marking phospholipid head groups (``P`` for PC/PE/PI/PS/PG).
    relax : bool
        Optional rigid-body clash relaxation after insertion.

    Returns
    -------
    mda.Universe
        The LD system (bilayer + neutral-lipid core), box extended along z.
    """
    if isinstance(core, (int, np.integer)):
        core = {"TRIO": int(core)}
    core = {k: int(v) for k, v in core.items() if int(v) > 0}
    if not core:
        raise ValueError("core must name at least one species, e.g. {'TRIO': 150}")
    n_total = sum(core.values())
    if core_source is None:
        core_source = DEFAULT_CORE_SOURCE

    u = bilayer if isinstance(bilayer, mda.Universe) else mda.Universe(bilayer)
    # make molecules whole and inside the box so leaflet/water assignment is clean
    u.atoms.wrap(compound="residues")
    box = u.dimensions.copy()
    area = box[0] * box[1]

    center = _leaflet_split_z(u, head_name)

    if core_thickness is None:
        vol = sum(n * CORE_VOLUME_A3.get(sp, TRIO_VOLUME_A3) for sp, n in core.items())
        core_thickness = vol / (area * packing_fraction)
    t = float(core_thickness)

    # --- assign every residue to the top or bottom half by its COM z, then
    #     translate the two halves apart to open a gap at the midplane ---
    res_com_z = np.array([r.atoms.center_of_mass()[2] for r in u.residues])
    top_mask = res_com_z >= center
    top_atoms = u.residues[top_mask].atoms
    bot_atoms = u.residues[~top_mask].atoms
    top_atoms.positions += np.array([0.0, 0.0, +t / 2.0])
    bot_atoms.positions += np.array([0.0, 0.0, -t / 2.0])

    # --- place the core species on a shared jittered 3-D grid in the gap ---
    templates = extract_templates(core_source, list(core.keys()), align_axis=True)
    rng = np.random.default_rng(seed)
    lateral_span = max(float(max((tpl.positions.max(0) - tpl.positions.min(0))[:2]))
                       for tpl in templates.values())
    z_span = max(float((tpl.positions.max(0) - tpl.positions.min(0))[2])
                 for tpl in templates.values())
    cell = 0.85 * lateral_span
    nx = max(1, int(box[0] // cell))
    ny = max(1, int(box[1] // cell))
    nz = max(1, int(np.ceil(n_total / (nx * ny))))
    sx, sy = box[0] / nx, box[1] / ny
    inset = min(0.30 * t, 0.25 * z_span)
    z_lo, z_hi = center - t / 2 + inset, center + t / 2 - inset
    sz = (z_hi - z_lo) / nz if nz else 0.0
    jit = 0.2 * min(sx, sy, max(sz, 1.0))

    # interleave the species across grid slots so the mix is spatially uniform
    species_bag = [sp for sp, n in core.items() for _ in range(n)]
    rng.shuffle(species_bag)
    slots = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    rng.shuffle(slots)

    core_mols = []  # (resname, names, coords)
    for (i, j, k), sp in zip(slots[:n_total], species_bag):
        tpl = templates[sp]
        coords = tpl.randomly_rotated(rng)
        centre = np.array([
            (i + 0.5) * sx + rng.uniform(-jit, jit),
            (j + 0.5) * sy + rng.uniform(-jit, jit),
            z_lo + (k + 0.5) * sz + rng.uniform(-jit, jit),
        ])
        core_mols.append((sp, list(tpl.names), coords + centre))

    # group by species in file order so [ molecules ] stays contiguous per type
    species_order = list(core.keys())
    core_mols.sort(key=lambda m: species_order.index(m[0]))

    ld = _merge_bilayer_and_core(u, core_mols, box, t)

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


def _merge_bilayer_and_core(u, core_mols, box, t):
    """Combine the shifted bilayer with the inserted core molecules.

    ``core_mols`` is a list of ``(resname, names, coords)`` tuples, already
    grouped by species so [ molecules ] stays contiguous.
    """
    n_core = len(core_mols)
    core_counts = [len(names) for _, names, _ in core_mols]
    total = u.atoms.n_atoms + sum(core_counts)

    n_res = u.residues.n_residues + n_core
    counts = [r.atoms.n_atoms for r in u.residues] + core_counts
    atom_resindex = np.repeat(np.arange(n_res), counts)

    new = mda.Universe.empty(n_atoms=total, n_residues=n_res,
                             atom_resindex=atom_resindex, trajectory=True)
    core_names = [nm for _, names, _ in core_mols for nm in names]
    new.add_TopologyAttr("names", list(u.atoms.names) + core_names)
    new.add_TopologyAttr("resnames",
                         list(u.residues.resnames) + [rn for rn, _, _ in core_mols])
    new.add_TopologyAttr("resids", np.arange(1, n_res + 1))

    coords = np.vstack([u.atoms.positions] + [c for _, _, c in core_mols])
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


def bilayer_to_ld(bilayer, core, out_gro=None, out_top=None,
                  core_source=None, toppar_includes=None, **kwargs):
    """Build an LD from a bilayer and (optionally) write .gro + .top.

    ``core`` is an int (TRIO only) or a dict, e.g. ``{"TRIO": 120, "CHYO": 30}``.
    Extra keyword arguments are passed to :func:`insert_core_into_bilayer`
    (``core_thickness``, ``packing_fraction``, ``seed``, ``relax`` ...).
    Returns the resulting Universe.
    """
    ld = insert_core_into_bilayer(bilayer, core, core_source=core_source, **kwargs)
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
                   [f"toppar/{ITP_FILENAME.get(r, r + '.itp')}" for r in present]
    with open(path, "w") as fh:
        for inc in includes:
            fh.write(f'#include "{inc}"\n')
        fh.write(f"\n[ system ]\n{system_name}\n\n")
        fh.write(molecules_section(universe))
    return path
