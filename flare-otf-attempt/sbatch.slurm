#!/bin/bash
#SBATCH --job-name=flare_training
#SBATCH --partition=normal_q
#SBATCH --account=aiaa2026
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --mail-user=cganley2@vt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

source ~/.bashrc
module load CMake/3.18.4 imkl GCC/12.3.0 QuantumESPRESSO/7.4-foss-2024a

conda init
conda activate flare

#export OMP_NUM_THREADS=32
python sandbox-script-otf.py
