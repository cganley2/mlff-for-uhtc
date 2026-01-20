#!/bin/bash
#SBATCH --job-name=zrb2-o-H-relax
#SBATCH --partition=normal_q
#SBATCH --account=aiaa2026
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --mail-user=cganley2@vt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load QuantumESPRESSO/7.4-foss-2024a

mpirun -np 64 pw.x -npool 2 -in zrb2-o-H-relax.in > zrb2-o-H-relax.out