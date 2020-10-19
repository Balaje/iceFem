#!/bin/bash

#PBS -l select=1:ncpus=16:mpiprocs=16:mem=16GB
#PBS -l walltime=10:00:00
#PBS -k oe
#PBS -q testq

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 1 -hsize 0.02 -nborders 6 -Tr 10000 -Ti 0 -iter 0 -notchWidth 0.005 -notchHeight 0.0
#mpirun -np 16 FreeFem++-mpi -nw -v 0 simple5.edp -Tr 2000 -Ti 0 -iter 0 -hsize 0.005 -N 10 -isSplit 1 -notchWidth 0.005 -notchHeight 0.0
exit 0
