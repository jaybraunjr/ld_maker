#!/bin/bash
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ld_equil
#SBATCH --out=ld_equil.out
#SBATCH --gres=gpu:2080ti:1
#SBATCH --account=swanson-gpu-np
#SBATCH --partition=swanson-gpu-np

# --- GROMACS module: adjust to whatever is available on CHPC ---
# (your CHARMM jobs use `module load charmm/46b2`; do `module spider gromacs`
#  to find the exact GROMACS module name/version on the cluster)
module load gcc/8.5.0 openmpi/4.1.1
module load gromacs            # e.g. gromacs/2022.x  -- FIX to real module name
# If GROMACS isn't a module, source your build's GMXRC instead:
# source /path/to/gromacs/bin/GMXRC

set -e
cd "$SLURM_SUBMIT_DIR"

ref=pure_ld_neutral
NT=${SLURM_CPUS_PER_TASK:-8}

# ---- step6.0 minimization ----
gmx grompp -f step6.0_minimization.mdp -o min.tpr \
    -c ${ref}.gro -r ${ref}.gro -p topol.top -n index.ndx -maxwarn 5
gmx mdrun -deffnm min -ntmpi 1 -ntomp ${NT} -nb gpu

# ---- step6.1 -> step6.6 equilibration ----
prev=min
for i in 1 2 3 4 5 6; do
    step=step6.$i
    gmx grompp -f ${step}_equilibration.mdp -o ${step}.tpr \
        -c ${prev}.gro -r ${ref}.gro -p topol.top -n index.ndx -maxwarn 5
    # GPU offload: nonbonded + PME on the GPU (single rank)
    gmx mdrun -deffnm ${step} -ntmpi 1 -ntomp ${NT} -nb gpu -pme gpu
    prev=${step}
done

echo "EQUILIBRATION COMPLETE -> final structure: step6.6.gro"
# Production (step7) can then be started from step6.6.gro / step6.6.cpt
