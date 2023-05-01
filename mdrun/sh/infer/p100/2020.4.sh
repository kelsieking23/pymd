#!/bin/bash
#SBATCH --nodes=
#SBATCH --ntasks-per-node=
#SBATCH -p p100_normal_q
#SBATCH -t 
#SBATCH --gres=gpu:1
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --job-name=
#SBATCH -o %x-slurm-%j.out

export MODULEPATH=$MODULEPATH:/projects/bevanlab/software/infer/modules/modules/infer-skylake/all
module load gromacs-p100/2020.4

