#!/bin/bash
#SBATCH --account=aiaa2026                # your account allocation
#SBATCH --job-name=2000K  # name of slurm job
#SBATCH --partition=a30_normal_q        # partition
#SBATCH --nodes=2                       # node count
#SBATCH --ntasks-per-node=8             # number of tasks per node
#SBATCH --cpus-per-task=1               # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=256G                       # total memory per node
#SBATCH --gres=gpu:4                    # number of gpus per node
#SBATCH --time=48:00:00                 # total run time limit (HH:MM:SS)

module reset
module load QuantumESPRESSO/7.4.1-NVHPC-24.9-CUDA-12.6.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpiexec -np $SLURM_NTASKS -x OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK pw.x -input zrb2-slab-001-o-aimd-2000K.in -npool 4 > zrb2-slab-001-o-aimd-2000K.out
