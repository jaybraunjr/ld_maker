"""Construct a planar lipid-droplet model (PL leaflet / TRIO core / PL leaflet).

The geometry is the flat "trilayer" slab used throughout this project: a neutral
lipid core (TRIO) sandwiched between two phospholipid leaflets whose head groups
face outward toward the (optional) water.  Molecules are stamped onto grids from
the single-molecule templates in :mod:`ld_maker.templates`; the result is meant
to be energy-minimised in GROMACS before use.
"""

import numpy as np
import MDAnalysis as mda

from .templates import extract_templates
from .relax import relax_clashes


def _leaflet_composition_counts(total, composition):
    """Turn a {resname: weight} dict into integer counts summing to ``total``."""
    names = list(composition)
    weights = np.array([composition[n] for n in names], dtype=float)
    weights /= weights.sum()
    counts = np.floor(weights * total).astype(int)
    # hand out the remainder to the largest-weight residues
    for i in np.argsort(-weights)[: total - counts.sum()]:
        counts[i] += 1
    return dict(zip(names, counts.tolist()))


class PlanarLDBuilder:
    """Build a planar LD slab from molecule templates.

    Parameters
    ----------
    templates : dict[str, Template]
    apl : float
        Target area per phospholipid in angstrom^2 (sets the leaflet grid pitch).
    seed : int
        RNG seed, so a given set of parameters reproduces the same structure.
    """

    def __init__(self, templates, apl=65.0, seed=0):
        self.templates = templates
        self.apl = float(apl)
        self.rng = np.random.default_rng(seed)

    # -- helpers ----------------------------------------------------------
    def _pl_spacing(self):
        # grid pitch from APL, but never tighter than the widest lipid footprint
        pitch = np.sqrt(self.apl)
        widest = max(
            self.templates[r].lateral_size
            for r in self.templates
            if r != "TRIO" and r in self.templates
        )
        return max(pitch, 0.55 * widest)

    def _stamp_leaflet(self, nx, ny, spacing, composition, z_interface, top):
        """Place one leaflet; returns (list_of_(resname, names, coords))."""
        counts = _leaflet_composition_counts(nx * ny, composition)
        bag = []
        for rn, c in counts.items():
            bag.extend([rn] * c)
        self.rng.shuffle(bag)

        placed = []
        jitter = 0.15 * spacing
        for idx, rn in enumerate(bag):
            gx, gy = idx % nx, idx // nx
            tpl = self.templates[rn]
            # top leaflet keeps head-up pose; bottom leaflet is flipped head-down
            coords = tpl.oriented(flip=not top, rng=self.rng)
            # align the tail end (facing the core) to the core interface
            if top:
                shift_z = z_interface - coords[:, 2].min()
            else:
                shift_z = z_interface - coords[:, 2].max()
            cx = (gx + 0.5) * spacing + self.rng.uniform(-jitter, jitter)
            cy = (gy + 0.5) * spacing + self.rng.uniform(-jitter, jitter)
            coords = coords + np.array([cx, cy, 0.0])
            coords[:, 2] += shift_z
            placed.append((rn, tpl.names, coords))
        return placed

    def _stamp_core(self, box_x, box_y, z_lo, z_hi, n_trio):
        """Place ``n_trio`` TRIO molecules on a jittered 3D grid inside the core."""
        tpl = self.templates["TRIO"]
        # cell size from the molecular volume proxy (cube root of bounding box)
        span = tpl.positions.max(axis=0) - tpl.positions.min(axis=0)
        # base the lateral cell on the molecule's in-plane size (with headroom)
        # so randomly-rotated rods start only mildly overlapped
        cell = 0.85 * float(max(span[0], span[1]))
        nx = max(1, int(box_x // cell))
        ny = max(1, int(box_y // cell))
        nz = max(1, int(np.ceil(n_trio / (nx * ny))))
        sx, sy = box_x / nx, box_y / ny
        # inset the centres so molecules don't poke far past the core interfaces
        inset = min(0.35 * (z_hi - z_lo), 0.25 * float(span[2]))
        z_lo_c, z_hi_c = z_lo + inset, z_hi - inset
        sz = (z_hi_c - z_lo_c) / nz if nz else 0.0
        jitter = 0.2 * min(sx, sy, max(sz, 1.0))

        placed = []
        slots = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
        self.rng.shuffle(slots)
        for (i, j, k) in slots[:n_trio]:
            coords = tpl.randomly_rotated(self.rng)
            centre = np.array([
                (i + 0.5) * sx + self.rng.uniform(-jitter, jitter),
                (j + 0.5) * sy + self.rng.uniform(-jitter, jitter),
                z_lo_c + (k + 0.5) * sz + self.rng.uniform(-jitter, jitter),
            ])
            placed.append(("TRIO", tpl.names, coords + centre))
        return placed

    # -- public API -------------------------------------------------------
    def build(
        self,
        nx,
        ny,
        n_trio,
        pl_composition=None,
        core_thickness=40.0,
        water_thickness=0.0,
        relax=True,
        d_min=1.6,
        relax_iterations=150,
    ):
        """Assemble the slab and return an ``mda.Universe``.

        Parameters
        ----------
        nx, ny : int
            Lipids per leaflet along x and y (leaflet holds ``nx * ny`` each).
        n_trio : int
            TRIO molecules in the core.
        pl_composition : dict[str, float] or None
            Relative amounts of each phospholipid, e.g.
            ``{"POPC": 176, "DOPE": 74, "SAPI": 20}``.  Defaults to POPC only.
        core_thickness : float
            Thickness of the TRIO core slab in angstrom.
        water_thickness : float
            Padding of empty space (for solvation) added above and below in
            angstrom.  No water molecules are added; solvate downstream.
        relax : bool
            Run rigid-body clash relaxation so the structure survives MD
            minimisation (recommended; expands the box slightly).
        d_min : float
            Target minimum inter-molecular contact for relaxation, angstrom.
        """
        if pl_composition is None:
            pl_composition = {"POPC": 1.0}
        missing = [r for r in list(pl_composition) + ["TRIO"] if r not in self.templates]
        if missing:
            raise ValueError(f"missing templates for: {missing}")

        spacing = self._pl_spacing()
        box_x, box_y = nx * spacing, ny * spacing

        # leaflet height = tallest PL template along z; leaves room so the
        # outward-facing heads sit above the bottom water pad (no negative z).
        leaflet_h = max(
            float(self.templates[r].positions[:, 2].ptp()) for r in pl_composition
        )
        core_lo = water_thickness + leaflet_h
        core_hi = core_lo + core_thickness

        molecules = []
        molecules += self._stamp_core(box_x, box_y, core_lo, core_hi, n_trio)
        molecules += self._stamp_leaflet(nx, ny, spacing, pl_composition, core_hi, top=True)
        molecules += self._stamp_leaflet(nx, ny, spacing, pl_composition, core_lo, top=False)

        if relax:
            # oversize the periodic box during relaxation so molecules have room
            # to spread; the final box is recomputed from the relaxed extent.
            pad = np.array([box_x, box_y, core_hi + leaflet_h + 2 * water_thickness])
            relax_clashes(
                molecules,
                box=[pad[0] * 1.3, pad[1] * 1.3, pad[2] * 1.3, 90.0, 90.0, 90.0],
                d_min=d_min,
                iterations=relax_iterations,
                lateral_only_resnames=tuple(pl_composition),
            )

        return self._to_universe(molecules, box_x, box_y, water_thickness)

    def _to_universe(self, molecules, box_x, box_y, water_thickness):
        # order residues so like molecules are contiguous (TRIO, then each PL)
        order = ["TRIO"] + [r for r in self.templates if r != "TRIO"]
        molecules.sort(key=lambda m: order.index(m[0]) if m[0] in order else 99)

        total_atoms = sum(len(m[1]) for m in molecules)
        n_res = len(molecules)
        u = mda.Universe.empty(
            n_atoms=total_atoms,
            n_residues=n_res,
            atom_resindex=np.repeat(np.arange(n_res), [len(m[1]) for m in molecules]),
            trajectory=True,
        )
        u.add_TopologyAttr("names", [n for m in molecules for n in m[1]])
        u.add_TopologyAttr("resnames", [m[0] for m in molecules])
        u.add_TopologyAttr("resids", np.arange(1, n_res + 1))

        coords = np.vstack([m[2] for m in molecules])
        # relaxation may have pushed molecules laterally; derive the periodic box
        # from the actual extent and shift so the lower corner sits at the origin
        # (with a water pad kept in z).  A half-grid margin keeps xy tiling clean.
        lo = coords.min(axis=0)
        coords[:, 0] -= lo[0]
        coords[:, 1] -= lo[1]
        coords[:, 2] -= lo[2] - water_thickness
        u.atoms.positions = coords

        margin = 0.5 * self._pl_spacing()
        box_x = max(box_x, coords[:, 0].max() + margin)
        box_y = max(box_y, coords[:, 1].max() + margin)
        zmax = coords[:, 2].max() + water_thickness
        u.dimensions = [box_x, box_y, zmax, 90.0, 90.0, 90.0]
        return u


def build_from_reference(reference, nx, ny, n_trio, pl_composition=None,
                         core_thickness=40.0, water_thickness=15.0, apl=65.0,
                         seed=0, relax=True, d_min=1.6):
    """Convenience: extract templates from ``reference`` and build in one call."""
    if pl_composition is None:
        pl_composition = {"POPC": 176, "DOPE": 74, "SAPI": 20}
    resnames = list(pl_composition) + ["TRIO"]
    templates = extract_templates(reference, resnames, align_axis=True)
    builder = PlanarLDBuilder(templates, apl=apl, seed=seed)
    return builder.build(nx, ny, n_trio, pl_composition, core_thickness,
                         water_thickness, relax=relax, d_min=d_min)
