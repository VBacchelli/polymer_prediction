#!/bin/bash
#SBATCH -n 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --gpus=1
#SBATCH --job-name=poly_job_$1
#SBATCH --error=job.error
#SBATCH --output=job.output
#---------------------------------------------------------------------------------------

module load gcc-11.2.1/gcc gcc-11.2.1/ompi-4.1.2 slurm miniconda3/python3.10

export OMP_NUM_THREADS=16

conda init bash
conda activate sim

NAME=$1
SMILES=$2

echo "Running $NAME"

python simulate_polymer.py --name "$NAME" --smiles "$SMILES"
