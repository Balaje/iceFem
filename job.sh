#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=10:00:00
#PBS -k oe

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
#mpirun -np 4 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 1 -hsize 0.02 -nborders 6 -Tr 10000 -Ti 0 -iter 0 -notchWidth 0.002 -notchHeight 0.0
mpirun -np 8 FreeFem++-mpi -nw -v 0 simple5.edp -Tr 50 -Ti 0 -iter 0 -hsize 0.01 -N 10 -isSplit 1 -notchWidth 0.005 -notchHeight 0.0
exit 0
