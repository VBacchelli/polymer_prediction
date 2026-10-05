#!/bin/bash
#SBATCH -n 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --gpus=1
#SBATCH --job-name=analysis_$1
#SBATCH --error=analysisbulk.err
#SBATCH --output=analysisbulk.out

module load gcc-11.2.1/gcc gcc-11.2.1/ompi-4.1.2 slurm miniconda3/python3.10

export OMP_NUM_THREADS=16

source ~/.bashrc
conda activate sim

PID=$1

echo "======================================"
echo " Running analysis for $PID"
echo "======================================"

python run_bulk_analysis.py "$PID"
