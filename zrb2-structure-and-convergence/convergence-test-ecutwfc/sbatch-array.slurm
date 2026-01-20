#!/bin/bash
#SBATCH --job-name=zrb2-ecutwfc-array
#SBATCH --output=%j_%a.slurm-out
#SBATCH --error=%j_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=14
#SBATCH --array=0-7
module load QuantumESPRESSO/7.4-foss-2024a
# Get the ecutwfc value for this array task
ecut=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ecut_list.tmp)
cd ecutwfc-$ecut
mpirun -np 14 pw.x -npool 14 -in qe.in > qe.out
echo ecutwfc-$ecut
