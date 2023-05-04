#!/bin/bash
#SBATCH --nodes=
#SBATCH --ntasks-per-node=
#SBATCH --account=cascades
#SBATCH -p v100_normal_q
#SBATCH --gres=gpu:2
#SBATCH -t 
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --job-name=
#SBATCH --export=NONE

module load gcc/7.3.0 
module load cuda/10.0.130 
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/groups/bevanlab/software/cascades/gromacs/2020.3/lib64
source groups/bevanlab/software/cascades/gromacs/2020.3/bin/GMXRC