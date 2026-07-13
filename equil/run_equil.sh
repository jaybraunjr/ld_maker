#!/bin/bash
# CHARMM-GUI 6-stage equilibration for the pure lipid droplet.
# Chain: min.gro -> step6.1 -> ... -> step6.6, restraints referenced to the
# initial structure (pure_ld_neutral.gro), exactly like the CHARMM-GUI README.
set -e
source /usr/local/gromacs/bin/GMXRC
cd "$(dirname "$0")"          # run from wherever this script lives

ref=pure_ld_neutral          # restraint reference (initial structure)
NT=$(nproc)
[ "$NT" -gt 8 ] && NT=8       # cap threads; diminishing returns + stability

# step6.0 minimization first (idempotent: skipped if min.gro already exists)
if [ ! -f min.gro ]; then
    echo "=================  step6.0 minimization  ($(date))  ================="
    gmx grompp -f step6.0_minimization.mdp -o min.tpr \
        -c ${ref}.gro -r ${ref}.gro -p topol.top -n index.ndx -maxwarn 5
    gmx mdrun -deffnm min -ntmpi 1 -ntomp ${NT}
fi

prev=min
for i in 1 2 3 4 5 6; do
    step=step6.$i
    echo "=================  $step  ($(date))  ================="
    gmx grompp -f ${step}_equilibration.mdp -o ${step}.tpr \
        -c ${prev}.gro -r ${ref}.gro -p topol.top -n index.ndx -maxwarn 5
    gmx mdrun -deffnm ${step} -ntmpi 1 -ntomp ${NT}
    prev=${step}
    echo "-----------------  $step DONE  ($(date))  -----------------"
done
echo "EQUILIBRATION_COMPLETE $(date)"
