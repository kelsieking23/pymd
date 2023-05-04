#!/bin/bash
#SBATCH -N1 --ntasks-per-node=32 --gres=gpu:1
#SBATCH -p t4_normal_q
#SBATCH -t 72:00:00
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --job-name=

module load fosscuda
module load GROMACS
