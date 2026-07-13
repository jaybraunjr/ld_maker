"""Single-molecule coordinate templates extracted from an equilibrated structure.

The maker doesn't ship coordinate files for TRIO/POPC/etc.  Instead it pulls one
copy of each residue straight out of one of your own equilibrated systems (e.g.
``modules/run.gro``), so the conformations are already physical.  A template is
just the atoms of a single residue, centred on its geometric centre and, for the
lipids, rotated so the long molecular axis points along +z.
"""

import numpy as np
import MDAnalysis as mda


# Atom used to find the "head" end of each lipid so we can orient it head-up.
# Anything not listed falls back to the principal-axis sign heuristic.
HEAD_ATOMS = {
    "POPC": "P",
    "DOPE": "P",
    "SAPI": "P",
}


def _principal_axis(positions):
    """Return the unit eigenvector of the largest moment of the coordinate set."""
    centred = positions - positions.mean(axis=0)
    cov = np.cov(centred.T)
    evals, evecs = np.linalg.eigh(cov)
    return evecs[:, np.argmax(evals)]


def _rotation_between(a, b):
    """Rotation matrix that maps unit vector ``a`` onto unit vector ``b``."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.linalg.norm(v) < 1e-8:
        # already (anti)parallel
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


class Template:
    """A single molecule's atoms, ready to be stamped into a grid.

    Attributes
    ----------
    resname : str
    names : list[str]      atom names, in file order
    positions : np.ndarray (N, 3)  centred coordinates in angstrom
    head_up : bool         True if oriented head toward +z (lipids only)
    """

    def __init__(self, resname, names, positions, head_up=False):
        self.resname = resname
        self.names = list(names)
        self.positions = np.asarray(positions, dtype=float)
        self.head_up = head_up

    @property
    def n_atoms(self):
        return len(self.names)

    @property
    def lateral_size(self):
        """Max in-plane (xy) extent, used to pick a safe grid spacing."""
        xy = self.positions[:, :2]
        return float((xy.max(axis=0) - xy.min(axis=0)).max())

    def oriented(self, flip=False, rng=None, random_z_spin=True):
        """Return a fresh copy of the coordinates for stamping.

        flip
            Mirror through the xy plane (z -> -z) to make a bottom-leaflet copy.
        random_z_spin
            Apply a random rotation about z so neighbouring copies don't line up.
        """
        pos = self.positions.copy()
        if random_z_spin and rng is not None:
            theta = rng.uniform(0, 2 * np.pi)
            ct, st = np.cos(theta), np.sin(theta)
            rz = np.array([[ct, -st, 0], [st, ct, 0], [0, 0, 1]])
            pos = pos @ rz.T
        if flip:
            pos = pos * np.array([1.0, 1.0, -1.0])
        return pos

    def randomly_rotated(self, rng):
        """Return coordinates under a uniformly random 3D rotation (for the core)."""
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        w, x, y, z = q
        rot = np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ])
        return self.positions @ rot.T


def extract_templates(reference, resnames, align_axis=True):
    """Pull one template molecule per residue name from a reference structure.

    Parameters
    ----------
    reference : str or mda.Universe
        Path to a .gro/.pdb (or a Universe) that contains at least one copy of
        each requested residue.
    resnames : iterable[str]
    align_axis : bool
        Rotate each lipid so its principal axis is along +z and its head atom
        (see ``HEAD_ATOMS``) points toward +z.  Leave False to keep the raw pose.

    Returns
    -------
    dict[str, Template]
    """
    u = reference if isinstance(reference, mda.Universe) else mda.Universe(reference)
    templates = {}
    for rn in resnames:
        sel = u.select_atoms(f"resname {rn}")
        if not len(sel):
            raise ValueError(f"reference has no residue named {rn!r}")

        # scan every copy and keep the most compact (smallest lateral footprint
        # once stood upright) so grid stamping produces the fewest clashes
        best = None
        candidates = sel.residues[: min(len(sel.residues), 200)]
        for res in candidates:
            ag = res.atoms
            pos = ag.positions - ag.positions.mean(axis=0)
            head_up = False
            if align_axis:
                axis = _principal_axis(pos)
                pos = pos @ _rotation_between(axis, np.array([0.0, 0.0, 1.0])).T
                head = HEAD_ATOMS.get(rn)
                if head is not None and head in ag.names:
                    if pos[list(ag.names).index(head), 2] < pos[:, 2].mean():
                        pos = pos * np.array([1.0, 1.0, -1.0])
                    head_up = True
            footprint = float((pos[:, :2].max(axis=0) - pos[:, :2].min(axis=0)).max())
            if best is None or footprint < best[0]:
                best = (footprint, ag.names, pos, head_up)

        _, names, pos, head_up = best
        templates[rn] = Template(rn, names, pos, head_up=head_up)
    return templates
