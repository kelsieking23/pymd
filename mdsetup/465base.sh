#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH -t 144:00:000
#SBATCH --account==cascades
#SBATCH -p normal_q
#SBATCH -A bevanlab
#SBATCH --job-name=tax5_400_450
#SBATCH --mail-user=kelsieking23@vt.edu
#SBATCH --mail-type=all
#SBATCH --export=NONE
module load intel
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/bevanlab/software/cascades/fftw/3.3.8/lib:/home/kelsieking23/software/gromacs/4.6.5/bin
source /home/kelsieking23/software/gromacs/4.6.5/bin/GMXRC
