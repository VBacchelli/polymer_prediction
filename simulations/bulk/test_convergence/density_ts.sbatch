#!/bin/bash
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

module load gcc-11.2.1/gcc gcc-11.2.1/ompi-4.1.2 slurm miniconda3/python3.10

source ~/.bashrc
conda activate sim

PID=$1

echo "======================================"
echo " Running density extraction for $PID"
echo "======================================"

python extract_density_timeseries.py $PID
