#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 
#SBATCH --account==cascades
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --job-name=
#SBATCH --mail-user=kelsieking23@vt.edu
#SBATCH --mail-type=all
#SBATCH --export=NONE

module load gcc/7.3.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/groups/bevanlab/software/cascades/gromacs/2018.1/lib64
export GMX=/groups/bevanlab/software/cascades/gromacs/2018.1/bin
source /groups/bevanlab/software/cascades/gromacs/2018.1/bin/GMXRC