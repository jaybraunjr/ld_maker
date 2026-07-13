"""Rigid-body clash relaxation for freshly stamped grids.

Real lipids pulled from a bilayer are bent and laterally splayed, so stamping
them onto an area-per-lipid grid leaves hard atomic overlaps that would blow up
an MD minimiser.  This module nudges whole molecules apart (translating each as a
rigid body, never distorting the conformation) until the closest inter-molecular
contact is comfortably above ``d_min``.  The box expands slightly; a normal NPT
run compresses it back to the target density.
"""

import numpy as np
from MDAnalysis.lib.distances import capped_distance


def relax_clashes(molecules, box, d_min=1.6, iterations=60, step=0.4,
                  lateral_only_resnames=(), verbose=False):
    """In-place rigid-body separation of clashing molecules.

    Parameters
    ----------
    molecules : list[(resname, names, coords ndarray)]
        As produced by the builder; ``coords`` arrays are modified in place.
    box : array-like (6,)
        MDAnalysis box (a, b, c, alpha, beta, gamma) for minimum-image contacts.
    d_min : float
        Target minimum inter-molecular atom-atom distance (angstrom).
    iterations : int
        Maximum relaxation sweeps.
    step : float
        Fraction of each overlap distance to move per sweep (0-1); damped.
    lateral_only_resnames : tuple[str]
        Residues that may only be pushed in x/y (keeps leaflets layered).
    """
    n = len(molecules)
    sizes = np.array([len(m[1]) for m in molecules])
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    mol_of_atom = np.repeat(np.arange(n), sizes)
    lateral = np.array(
        [m[0] in lateral_only_resnames for m in molecules], dtype=bool
    )
    box = np.asarray(box, dtype=np.float32)

    for it in range(iterations):
        coords = np.vstack([m[2] for m in molecules]).astype(np.float32)
        pairs, dists = capped_distance(
            coords, coords, max_cutoff=d_min, box=box, return_distances=True
        )
        # keep genuine inter-molecular contacts once (a < b)
        mask = mol_of_atom[pairs[:, 0]] != mol_of_atom[pairs[:, 1]]
        mask &= pairs[:, 0] < pairs[:, 1]
        pairs, dists = pairs[mask], dists[mask]
        if len(pairs) == 0:
            if verbose:
                print(f"  relaxed: clash-free after {it} sweeps")
            break

        shift = np.zeros((n, 3), dtype=np.float64)
        ai, bi = pairs[:, 0], pairs[:, 1]
        ma, mb = mol_of_atom[ai], mol_of_atom[bi]
        d = np.maximum(dists, 1e-3)
        overlap = (d_min - d)
        # push each molecule along the vector from the contact partner
        direction = coords[ai] - coords[bi]
        direction -= box[:3] * np.round(direction / box[:3])   # min image
        direction /= np.linalg.norm(direction, axis=1, keepdims=True) + 1e-9
        push = (overlap * step)[:, None] * direction
        # accumulate (sum, not mean) so buried molecules still escape; the total
        # per-molecule displacement is capped for stability
        np.add.at(shift, ma, push)
        np.add.at(shift, mb, -push)
        shift[lateral, 2] = 0.0
        mag = np.linalg.norm(shift, axis=1)
        cap = 1.5
        over = mag > cap
        shift[over] *= (cap / mag[over])[:, None]

        moved = mag > 1e-6
        for mi in np.nonzero(moved)[0]:
            molecules[mi][2][:] = molecules[mi][2] + shift[mi]

        if verbose and (it % 10 == 0):
            print(f"  sweep {it}: {len(pairs)} clashing pairs, min d={d.min():.2f}")
    return molecules
