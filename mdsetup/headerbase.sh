#!/bin/bash
#SBATCH --nodes=
#SBATCH --ntasks-per-node=24
#SBATCH --gres=gpu:
#SBATCH --account=cascades
#SBATCH -p 
#SBATCH -t 
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --job-name=
#SBATCH --export=NONE

# load modules, export library path, load gromacs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/groups/bevanlab/software/cascades/gromacs/2019.3/lib64
