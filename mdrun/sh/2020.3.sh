#!/bin/bash
#SBATCH --nodes=
#SBATCH --ntasks-per-node=
#SBATCH -p normal_q
#SBATCH -t 
#SBATCH --account=cascades
#SBATCH -A bevanlab
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --job-name=
#SBATCH --export=NONE
# load modules, export library path, load gromacs

module load gcc/7.3.0
module load cuda/9.2.148
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kelsieking23/cascades/software/fftw/3.3.5/lib:/groups/bevanlab/software/cascades/software/gromacs/2020.3/lib64
source /groups/bevanlab/software/cascades/gromacs/2020.3/bin/GMXRC