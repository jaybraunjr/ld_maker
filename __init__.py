"""ld_maker -- build planar lipid-droplet (trilayer) starting structures.

A neutral-lipid TRIO core sandwiched between two phospholipid leaflets, stamped
onto grids from single-molecule templates taken from your own equilibrated
systems.  Output a .gro plus a GROMACS ``[ molecules ]`` topology; minimise in
GROMACS before running.

Quick start
-----------
>>> from ld_maker import build_from_reference, write_top
>>> u = build_from_reference("modules/run.gro", nx=12, ny=12, n_trio=200)
>>> u.atoms.write("planar_ld.gro")
>>> write_top(u, "planar_ld.top")
"""

from .templates import Template, extract_templates
from .builder import PlanarLDBuilder, build_from_reference
from .bilayer import insert_core_into_bilayer, bilayer_to_ld, replace_lipid
from .topology import molecule_counts, molecules_section, write_top
from .vmd import write_vmd

__all__ = [
    "Template",
    "extract_templates",
    "PlanarLDBuilder",
    "build_from_reference",
    "insert_core_into_bilayer",
    "bilayer_to_ld",
    "replace_lipid",
    "molecule_counts",
    "molecules_section",
    "write_top",
    "write_vmd",
]
