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

module reset
module load QuantumESPRESSO/7.4-foss-2024a
module load Miniconda3/24.7.1-0
module load imkl

source activate flare

which python
python sandbox-script-otf.py 