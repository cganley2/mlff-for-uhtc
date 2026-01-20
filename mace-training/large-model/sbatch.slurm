#!/bin/bash
#SBATCH --job-name=mace_training
#SBATCH --partition=normal_q
#SBATCH --account=aiaa2026
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --constraint=amd
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --mail-user=cganley2@vt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

ml Miniconda3
source activate mace

mkdir models
mace_run_train --config mace-config.yml
mkdir evals
mace_eval_configs --config /home/cganley2/scitech2026/pre-training-data/aimd123_train_1200frames.extxyz --model models/zr-b-o_run-123_stagetwo.model --output evals/mace-train-predictions.xyz
mace_eval_configs --config /home/cganley2/scitech2026/pre-training-data/aimd123_test_150frames.extxyz --model models/zr-b-o_run-123_stagetwo.model --output evals/mace-test-predictions.xyz
