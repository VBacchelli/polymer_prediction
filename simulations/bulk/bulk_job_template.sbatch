#!/bin/bash
#SBATCH -n 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --gpus=1

module load gcc-11.2.1/gcc gcc-11.2.1/ompi-4.1.2 slurm miniconda3/python3.10

export OMP_NUM_THREADS=16

source ~/.bashrc
conda activate sim

NAME=$1
SMILES=$2

echo "Running $NAME"
echo "SMILES: $SMILES"

python run_bulk_from_smiles.py --name "$NAME" --smiles "$SMILES"


echo "======================================"
echo " Cleaning unnecessary files for $PID"
echo "======================================"

BASE=/home/vbacchelli/BIO-SUSHY-workflows/FOOD_v2/simulations/$PID

# previous checkpoints (keeps last)
find $BASE -name "*_prev.cpt" -delete

