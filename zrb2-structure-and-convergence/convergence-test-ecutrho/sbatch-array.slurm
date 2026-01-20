#!/bin/bash
#SBATCH --job-name=zrb2-ecutrho-mult-array
#SBATCH --output=%j_%a.slurm-out
#SBATCH --error=%j_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --array=0-8
module load QuantumESPRESSO/7.4-foss-2024a
# Get the ecutrho value for this array task
ecutrho=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ecutrhomult_list.tmp)
cd ecutrho-mult-$ecutrho
pw.x -in qe.in > qe.out
echo ecutrho-mult-$ecutrho
