#!/bin/bash
# ============================================================================
# Full lipid-droplet workflow: bilayer -> TRIO-core LD -> EM -> equilibration -> NPT
#
# Builds a lipid droplet from a (CHYO/DOPE/POPC) phospholipid bilayer by
# inserting a TRIO core, then runs the CHARMM-GUI protocol end to end:
# energy minimization (step6.0), six-stage equilibration (step6.1-6.6) and
# NPT production (step7).
#
# Usage:
#   bash run_ld.sh BILAYER OUTDIR [N_TRIO] [CHYO_FRAC]
#
#     BILAYER    solvated bilayer .gro (CHARMM36). If it already contains CHYO
#                in the leaflets, leave CHYO_FRAC=0. Otherwise CHYO_FRAC>0 first
#                swaps that fraction of leaflet POPC for CHYO.
#     OUTDIR     run directory to create
#     N_TRIO     TRIO molecules in the core         (default 150)
#     CHYO_FRAC  fraction of leaflet POPC -> CHYO   (default 0)
#
# Env overrides: GMXRC, TOPPAR, NT (mdrun thread flags)
#
# Note: in-place POPC->CHYO replacement over-packs the leaflet even at low
# fractions and will NOT minimize. For CHYO in the leaflets, build the
# CHYO/DOPE/POPC bilayer in CHARMM-GUI (which sizes the box for the real
# per-lipid areas) and pass it here with CHYO_FRAC=0. On CPU the
# equilibration+NPT take hours -- use run_equil.slurm on a GPU node for production.
# ============================================================================
set -euo pipefail

BILAYER=${1:?need a bilayer .gro}
OUTDIR=${2:?need an output directory}
NTRIO=${3:-150}
CHYO_FRAC=${4:-0}

HERE=$(cd "$(dirname "$0")" && pwd)                 # ld_maker/equil
ROOT=$(cd "$HERE/../.." && pwd)                      # project root (parent of ld_maker)
: "${GMXRC:=/usr/local/gromacs/bin/GMXRC}"
: "${TOPPAR:=$ROOT/analysis/dcd_big/big_v2/toppar}"
: "${NT:=-ntmpi 1 -ntomp $(nproc)}"
export PYTHONPATH="$ROOT:${PYTHONPATH:-}"           # make 'ld_maker' importable
source "$GMXRC"

mkdir -p "$OUTDIR"; OUTDIR=$(cd "$OUTDIR" && pwd)    # absolutize

# --- 1. build the LD --------------------------------------------------------
BIL="$BILAYER"
if [ "$CHYO_FRAC" != "0" ]; then
  echo "[warn] in-place POPC->CHYO ($CHYO_FRAC) over-packs the leaflet and usually"
  echo "[warn] fails to minimize; prefer a CHARMM-GUI CHYO bilayer with CHYO_FRAC=0."
  echo "[build] swapping $CHYO_FRAC of leaflet POPC -> CHYO"
  python -m ld_maker --bilayer "$BILAYER" --replace "POPC:CHYO:$CHYO_FRAC" \
      --out "$OUTDIR/chyo_bilayer"
  BIL="$OUTDIR/chyo_bilayer.gro"
fi
echo "[build] inserting $NTRIO TRIO core"
python -m ld_maker --bilayer "$BIL" --core "TRIO:$NTRIO" --out "$OUTDIR/system" --vmd
mv "$OUTDIR/system.top" "$OUTDIR/topol.top"

# --- 2. stage run directory -------------------------------------------------
cp -r "$TOPPAR" "$OUTDIR/toppar"
cp "$HERE"/step6.*.mdp "$HERE"/step7_production.mdp "$OUTDIR/"
python -m ld_maker.index "$OUTDIR/system.gro" "$OUTDIR/index.ndx" --membrane=CHYO
cp "$OUTDIR/system.gro" "$OUTDIR/ref.gro"           # position-restraint reference

cd "$OUTDIR"
run() { echo "[$1] $2"; gmx grompp -f "$3" -c "$4" -r ref.gro ${5:+-t "$5"} \
          -p topol.top -n index.ndx -o "$6.tpr" -maxwarn 5 >/dev/null 2>&1;
        gmx mdrun -deffnm "$6" $NT; }

# --- 3. energy minimization -------------------------------------------------
run EM  "minimization"       step6.0_minimization.mdp  system.gro  ""        em

# --- 4. equilibration step6.1 -> step6.6 ------------------------------------
run EQ1 "NVT (gen. vel.)"    step6.1_equilibration.mdp em.gro      ""        eq1
prev=eq1
for s in 2 3 4 5 6; do
  run "EQ$s" "step6.$s"      step6.${s}_equilibration.mdp ${prev}.gro ${prev}.cpt eq${s}
  prev=eq${s}
done

# --- 5. NPT production (step7) ----------------------------------------------
run NPT "production"         step7_production.mdp      eq6.gro     eq6.cpt   npt

echo "DONE -> $OUTDIR/npt.gro   (view: cd $OUTDIR && vmd -e system.vmd)"
