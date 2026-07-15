# ld_maker

[![tests](https://github.com/jaybraunjr/ld_maker/actions/workflows/ci.yml/badge.svg)](https://github.com/jaybraunjr/ld_maker/actions/workflows/ci.yml)

Build **lipid-droplet starting structures** for GROMACS / CHARMM36 molecular
dynamics. An LD is a neutral-lipid core (TRIO ± cholesteryl ester) wrapped by two
phospholipid monolayers; `ld_maker` builds that geometry and hands off to the
standard six-stage CHARMM-GUI equilibration.

Templates come straight from your own equilibrated structures, so conformations
are physical and no external builder is required for the common cases. Docs:
**https://jaybraunjr.github.io/ld_maker/**

## Install

```bash
pip install -e .            # from the repo root
# or with test deps:
pip install -e ".[test]"
```
Requires Python ≥3.9, `numpy`, `scipy`, `MDAnalysis`. The mixed-core and
lipid-replacement packing steps additionally call `gmx` (GROMACS) if it is on
`PATH`; everything else is pure Python.

## What it builds

### 1. Bilayer → LD  (recommended)
Take a solvated CHARMM-GUI bilayer, pry the leaflets apart, and fill the midplane
with neutral lipid. The bilayer's water and ions are kept as-is.

```bash
# pure triolein core
python -m ld_maker --bilayer bil.gro --n-trio 150 --out runs/ld

# mixed triolein / cholesteryl-ester core (packed with gmx insert-molecules)
python -m ld_maker --bilayer bil.gro --core "TRIO:120,CHYO:30" --out runs/ld_ce
```
Composition-agnostic: whatever lipids/sterols are in the bilayer carry through, so
build exotic monolayers in CHARMM-GUI and convert here.

### 2. Swap lipids in the leaflets
Replace a fraction of a leaflet lipid with another in place — e.g. a surface
cholesteryl-ester model.

```bash
python -m ld_maker --bilayer bil.gro --replace "POPC:CHYO:0.5" --out runs/surface_ce
```
Bulky substitutions over-pack the leaflet; `ld_maker` warns when the result would
blow up minimization (see Limitations).

### 3. From-scratch grid build
Synthesize the whole PL / core / PL slab on a grid (lipids + core only; solvate
downstream).

```bash
python -m ld_maker --reference run.gro --nx 12 --ny 12 --n-trio 200 \
    --composition POPC:176,DOPE:74,SAPI:20 --out runs/planar
```

Each mode writes `<out>.gro` and `<out>.top` (CHARMM36 `#include`s + a
`[ molecules ]` section counted in file order).

## Python API

```python
from ld_maker import bilayer_to_ld, replace_lipid, build_from_reference

bilayer_to_ld("bil.gro", core={"TRIO": 120, "CHYO": 30}, out_gro="ld.gro", out_top="ld.top")
replace_lipid("bil.gro", {"POPC": ("CHYO", 0.5)}, out_gro="ce.gro", out_top="ce.top")
```

## End-to-end pipeline

`pipeline.sh` builds, stages a self-contained run directory (coordinates,
topology, `index.ndx`, CHARMM36 `toppar`, `step6/step7` `.mdp`s), and optionally
minimizes:

```bash
bash pipeline.sh bil.gro 150 runs/ld1 --minimize
cd runs/ld1 && sbatch run_equil.slurm      # or: bash run_equil.sh
```

For the **whole run in one shot** — build the LD from a bilayer, then EM →
equilibration → NPT production:

```bash
bash equil/run_ld.sh bil.gro runs/ld1 150      # BILAYER OUTDIR N_TRIO [CHYO_FRAC]
```
Pass a CHYO/DOPE/POPC bilayer (CHYO in the leaflets) with `CHYO_FRAC=0`; the
script inserts the TRIO core and runs `step6.0` → `step6.6` → `step7`. (In-place
`CHYO_FRAC>0` over-packs the leaflet — build that composition in CHARMM-GUI.)

Generate the thermostat/restraint index directly with
`python -m ld_maker.index system.gro index.ndx [--membrane=CHYO]`
(`--membrane` keeps *surface* CHYO in `MEMB` rather than the core `SOLU` group).

## Visualize (VMD)

Every build can emit a `.vmd` view script (`--vmd`), or generate one for any
structure after the fact:

```bash
python -m ld_maker --bilayer bil.gro --core "TRIO:120,CHYO:30" --out ld --vmd
python -m ld_maker.vmd ld.gro -o ld.vmd -t traj.xtc     # + optional trajectory
vmd -e ld.vmd
```
It sets one representation per component present — TRIO core (orange), CHYO
(purple), POPC/DOPE/SAPI, phosphates (silver spheres), ions, water — and wraps
the trajectory with PBCTools. `write_vmd(...)` is also available from the API.

## Tests

```bash
pytest -q
```
Fast, synthetic fixtures (no large files). Tests that need GROMACS skip
automatically when `gmx` isn't installed. CI runs the suite on Python 3.9/3.11/3.12.

## Package layout

| file | purpose |
|------|---------|
| `bilayer.py`   | bilayer→LD core insertion, mixed cores, `replace_lipid` |
| `builder.py`   | `PlanarLDBuilder` — from-scratch grid build |
| `templates.py` | extract compact single-molecule templates from a reference |
| `relax.py`     | rigid-body clash relaxation |
| `topology.py`  | GROMACS `[ molecules ]` / `.top` generation |
| `index.py`     | SOLU / MEMB / SOLV `index.ndx` generator |
| `pipeline.sh`  | one-command build → stage → minimize |
| `equil/`       | CHARMM-GUI `step6/step7` protocol + run scripts |
| `templates/`   | packaged single-molecule core templates (TRIO, CHYO) |

## Limitations

- Output is a **starting structure** — always energy-minimize before MD.
- Mixed / CHYO cores are packed with `gmx insert-molecules`, which needs a roomy
  gap (the tool loosens the packing automatically); NPT then compresses to density.
- `replace_lipid` places the substitute in place. A **high fraction of a bulky
  lipid** (e.g. 50% CHYO) over-packs a phospholipid leaflet and cannot minimize —
  the tool warns, and the right route for that density is CHARMM-GUI Membrane
  Builder with correct per-lipid areas.
