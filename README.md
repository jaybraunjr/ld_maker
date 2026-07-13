# ld_maker — planar lipid-droplet builder

Builds a **planar LD starting structure**: a neutral-lipid **TRIO** core sandwiched
between two **phospholipid leaflets** (POPC / DOPE / SAPI) whose head groups face
outward — the flat "trilayer" model used throughout this project.

Molecules are stamped onto grids using single-molecule templates taken straight
from one of your own equilibrated structures (e.g. `modules/run.gro`), so the
conformations are physical and no external tools (Packmol, CHARMM-GUI) are needed.
A rigid-body clash relaxation then nudges molecules apart until the structure is
safe to minimise in GROMACS.

## Install / requirements
Pure Python; needs only `MDAnalysis` and `numpy` (already in this env).

## Two ways to build

### 1. Bilayer → LD (recommended): insert TRIO into a real bilayer
Start from a solvated CHARMM-GUI phospholipid bilayer, pry the leaflets apart, and
fill the midplane with TRIO. The bilayer's water/ions are kept as-is. **Composition
-agnostic** — whatever lipids/sterols (cholesterol, mixed PLs) are in the bilayer
just carry through, so build exotic compositions in CHARMM-GUI and convert here.

```bash
python -m ld_maker \
    --bilayer modules/move_prot/drude/bil/bil.gro \
    --n-trio 150 \
    --trio-source modules/run.gro \
    --out bil2ld
```
Options: `--packing-fraction` (initial core fill, lower = looser/safer start that
NPT compresses; default 0.35), `--core-thickness` (override the auto-sized gap),
`--no-relax`, `--seed`. Python: `ld_maker.bilayer_to_ld(...)`.

### 2. From-scratch grid build
Synthesize the whole PL/TRIO/PL sandwich on a grid from single-molecule templates
(no bilayer needed; produces lipids+core only, solvate downstream).

```bash
python -m ld_maker \
    --reference modules/run.gro \
    --nx 12 --ny 12 \
    --n-trio 200 \
    --composition POPC:176,DOPE:74,SAPI:20 \
    --core-thickness 40 \
    --water 15 \
    --out planar_ld
```
Writes `<out>.gro` (coordinates) and `<out>.top` (GROMACS `[ molecules ]` section,
with CHARMM36 `#include`s pointing at `toppar/`).

Key options: `--apl` target area per phospholipid (Å²), `--d-min` minimum
inter-molecular contact target for relaxation (Å), `--no-relax` to skip
relaxation (fast but clashy), `--seed` for reproducibility.

## Python API
```python
from ld_maker import build_from_reference, write_top

u = build_from_reference(
    "modules/run.gro",
    nx=12, ny=12, n_trio=200,
    pl_composition={"POPC": 176, "DOPE": 74, "SAPI": 20},
    core_thickness=40.0, water_thickness=15.0,
)
u.atoms.write("planar_ld.gro")
write_top(u, "planar_ld.top")
```

For finer control, use the classes directly:
```python
from ld_maker import extract_templates, PlanarLDBuilder
templates = extract_templates("modules/run.gro", ["POPC", "DOPE", "SAPI", "TRIO"])
builder = PlanarLDBuilder(templates, apl=65.0, seed=0)
u = builder.build(nx=12, ny=12, n_trio=200,
                  pl_composition={"POPC": 176, "DOPE": 74, "SAPI": 20})
```

## What the geometry looks like
- Leaflet lipids: `nx * ny` per leaflet, heads pointing outward (+z top, −z bottom).
- Core: `n_trio` TRIO on a jittered 3-D grid between the leaflet tail regions.
- `water_thickness`: empty pad added above and below for you to **solvate**
  afterward (e.g. `gmx solvate` / add ions). No water or ions are placed here.

Validate with `validation_zprofile.png` (regenerate any time): two phospholipid
head-group peaks at the outer surfaces with the TRIO core between them.

## Recommended downstream workflow (GROMACS)
1. `gmx editconf` if you need to adjust the box / add vacuum for solvation.
2. `gmx solvate` to add TIP3 water in the pads; then add SOD/CLA ions.
3. Update the `[ molecules ]` counts in the `.top` for added water/ions.
4. **Energy minimise** (steepest descent) before any MD — the build is a grid
   start; minimisation removes residual close contacts.
5. Equilibrate (NVT then NPT with semi-isotropic pressure coupling); the box
   compresses from the relaxed grid to the target density.

## Files
| file | purpose |
|------|---------|
| `bilayer.py`   | **bilayer → LD**: insert a TRIO core between a bilayer's leaflets |
| `templates.py` | extract compact single-molecule templates from a reference |
| `builder.py`   | `PlanarLDBuilder` — stamp leaflets + core, set box (grid build) |
| `relax.py`     | rigid-body clash relaxation |
| `topology.py`  | GROMACS `[ molecules ]` / `.top` generation (from `make_big`) |
| `__main__.py`  | command-line interface (`--bilayer` or `--reference` mode) |

## Notes & limitations
- The maker places **lipids + core only**; water and ions are added downstream.
- Grid + relaxation gives a clash-light start (min contact ≈ target `--d-min`),
  **not** an equilibrated structure — always minimise first.
- The reference must contain at least one copy of every requested residue.
