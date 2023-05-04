#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --gres=gpu:
#SBATCH --account=cascades
#SBATCH -p v100_normal_q
#SBATCH -t 
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=kelsieking23@vt.edu
#SBATCH --job-name=
#SBATCH --export=NONE

# load modules, export library path, load gromacs
module load gcc/7.3.0
module load cuda/10.0.130
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/groups/bevanlab/software/cascades/gromacs/2019.3/lib64
source /groups/bevanlab/software/cascades/gromacs/2019.3/bin/GMXRC
