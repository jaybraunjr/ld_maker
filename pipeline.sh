#!/bin/bash
# =============================================================================
# ld_maker pipeline: bilayer  ->  lipid droplet  ->  minimized + equilibration-ready
# -----------------------------------------------------------------------------
# One command builds the LD, stages a self-contained run directory (coords,
# topology, index, mdps, toppar), and (optionally) runs energy minimization.
# Equilibration is then a single `bash run_equil.sh` (local) or
# `sbatch run_equil.slurm` (CHPC GPU).
#
# Usage:
#   bash pipeline.sh BILAYER.gro N_TRIO OUTDIR [--minimize]
#
# Example:
#   bash pipeline.sh modules/move_prot/drude/bil/bil.gro 150 runs/ld1 --minimize
# =============================================================================
set -euo pipefail

BILAYER=${1:?need a bilayer .gro (e.g. modules/move_prot/drude/bil/bil.gro)}
NTRIO=${2:?need number of TRIO to insert (e.g. 150)}
OUTDIR=${3:?need an output directory (e.g. runs/ld1)}
DO_MIN=${4:-}

HERE=$(cd "$(dirname "$0")" && pwd)          # ld_maker/
ROOT=$(cd "$HERE/.." && pwd)                  # project root
PYTHON=${PYTHON:-python}                       # python with MDAnalysis (override if needed)
TRIO_SRC=${TRIO_SRC:-$ROOT/modules/run.gro}   # TRIO template donor
TOPPAR_SRC=${TOPPAR_SRC:-$ROOT/analysis/dcd_big/big_v2/toppar}
MDP_SRC=${MDP_SRC:-$ROOT/ld_maker/equil}      # step6.*/step7 mdps live here

echo "[1/4] building lipid droplet: insert $NTRIO TRIO into $BILAYER"
mkdir -p "$OUTDIR"
# the equilibration scripts expect the reference named 'pure_ld_neutral' + topol.top
$PYTHON -m ld_maker --bilayer "$BILAYER" --n-trio "$NTRIO" \
    --trio-source "$TRIO_SRC" --out "$OUTDIR/pure_ld_neutral"
mv "$OUTDIR/pure_ld_neutral.top" "$OUTDIR/topol.top"

echo "[2/4] staging run directory in $OUTDIR"
cp -r "$TOPPAR_SRC" "$OUTDIR/toppar"
cp "$MDP_SRC"/step6.*.mdp "$MDP_SRC"/step7_production.mdp "$OUTDIR/" 2>/dev/null || true
cp "$HERE/equil/run_equil.sh" "$HERE/equil/run_equil.slurm" "$OUTDIR/" 2>/dev/null || true

echo "[3/4] writing index.ndx (SOLU/MEMB/SOLV groups)"
$PYTHON -m ld_maker.index "$OUTDIR/pure_ld_neutral.gro" "$OUTDIR/index.ndx"

if [ "$DO_MIN" == "--minimize" ]; then
    echo "[4/4] energy minimization (step6.0)"
    ( cd "$OUTDIR"
      set +u; source /usr/local/gromacs/bin/GMXRC; set -u
      gmx grompp -f step6.0_minimization.mdp -c pure_ld_neutral.gro -r pure_ld_neutral.gro \
          -p topol.top -n index.ndx -o min.tpr -maxwarn 5
      gmx mdrun -deffnm min -ntmpi 1 -ntomp "$(nproc)" )
    echo "    -> $OUTDIR/min.gro"
else
    echo "[4/4] skipped minimization (pass --minimize to run it)"
fi

cat <<EOF

Done. Run directory: $OUTDIR
  pure_ld_neutral.gro       built lipid droplet (also the restraint reference)
  topol.top                 topology
  index.ndx                 thermostat/restraint groups
  step6.*.mdp, step7_*.mdp  CHARMM-GUI protocol
  run_equil.sh              local equilibration  (bash run_equil.sh)
  run_equil.slurm           CHPC GPU submission  (sbatch run_equil.slurm)

Next: cd $OUTDIR && bash run_equil.sh      # or sbatch run_equil.slurm
EOF
