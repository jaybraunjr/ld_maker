# Pure lipid droplet — minimization + equilibration package

A protein-free lipid droplet built from your equilibrated `modules/run.gro`, ready
for the CHARMM-GUI 6-stage equilibration. Follows your `make_big` method (reuse an
already-solvated/ionized system, regenerate topology by residue order) rather than
solvating from scratch.

## How it was built
1. **Strip protein** from `modules/run.gro` (kept `not protein`) → 233,441 atoms,
   real TIP3 water + SOD/CLA ions retained (no re-solvation).
2. **Regenerate `topol.top`** by counting residues in file order (your make_big cell).
3. **Neutralize**: removing the −1 protein left net +1, so one Na⁺ was removed
   → `pure_ld_neutral.gro` (grompp confirms zero total charge).
4. **`index.ndx`** with thermostat groups mapped to the droplet's physical layers:
   - `SOLU` = TRIO core, `MEMB` = POPC/DOPE/SAPI, `SOLV` = TIP3/SOD/CLA
   - (protein-free, so the usual SOLU=protein group is repurposed — no empty groups)

## Verified locally (GROMACS 2022, WSL)
- **step6.0 minimization: converged**, PE = −2.26×10⁶ kJ/mol, Fmax < 1000 in 1547 steps.
- **step6.1 equilibration: MD-stable** (restrained NVT ran cleanly, no LINCS blow-ups).
- Local throughput ~0.94 ns/day (CPU-only) → full 1.5 ns eq ≈ 1.6 days locally,
  hence run it on the cluster.

## Run it

### On CHPC (recommended — GPU, finishes in <1 h)
```bash
sbatch run_equil.slurm      # fix the `module load gromacs` line first
```

### Locally (CPU, slow)
```bash
bash run_equil.sh           # sources GMXRC, runs min -> step6.1..6.6
```

Both chain: `-c` = previous stage, `-r` = `pure_ld_neutral.gro` (restraint reference),
`-n index.ndx`, exactly like the CHARMM-GUI README. Final structure: `step6.6.gro`,
then start production (`step7_production.mdp`) from it.

## File inventory
| file | role |
|------|------|
| `pure_ld_neutral.gro` | starting coords (protein-stripped, neutral) — also the `-r` reference |
| `topol.top` | regenerated topology (`[ molecules ]` in file order) |
| `index.ndx` | SOLU/MEMB/SOLV/SOLU_MEMB/SYSTEM groups |
| `toppar/` | CHARMM36 force field + lipid/ion/water itps (from big_v2) |
| `step6.0_minimization.mdp` | minimization |
| `step6.1..6.6_equilibration.mdp` | 6-stage restraint ramp-down |
| `step7_production.mdp` | production template |
| `run_equil.slurm` | CHPC GPU submission (min → 6.6) |
| `run_equil.sh` | local CPU driver (min → 6.6) |
| `min.gro`, `min.tpr`, `min.log` | completed minimization output |
| `step6.1_demo.*` | short 5 ps local equilibration snapshot (validation) |

## Notes
- Stripping the protein leaves a small cavity; the restraint ramp-down + semi-isotropic
  NPT closes it during equilibration.
- To make a *bigger* droplet, tile `pure_ld_neutral.gro` 2×2 in XY first
  (your `make_big` `tile_system_xy`), then regenerate `topol.top`/`index.ndx`.
