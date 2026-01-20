#!/bin/bash
#SBATCH --job-name=zrb2-kpts-array
#SBATCH --output=%j_%a.slurm-out
#SBATCH --error=%j_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --array=0-7
module load QuantumESPRESSO/7.4-foss-2024a
# Get the kpts value for this array task
kpts=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" kpts_list.tmp)
cd kpts-$kpts
pw.x -in qe.in > qe.out
echo kpts-$kpts
