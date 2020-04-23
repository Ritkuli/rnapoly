#!/bin/bash
#SBATCH -J $NAME
#SBATCH -t 0-4:00:00
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint='kepler|pascal|volta'

srun --gres=gpu $OXDNA $INPUT
