#!/bin/bash
#SBATCH -n 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gpus=1
#SBATCH --job-name=rotacf_$1
#SBATCH --output=rotacf.out
#SBATCH --error=rotacf.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1

module load gcc-11.2.1/gcc gcc-11.2.1/ompi-4.1.2 slurm miniconda3/python3.10

export OMP_NUM_THREADS=16

conda init bash
conda activate sim

python run_one_analysis.py "$1" dummy.mdp

