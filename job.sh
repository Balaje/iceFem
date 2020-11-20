#!/bin/bash

#PBS -l select=1:ncpus=16:mpiprocs=16:mem=16GB
#PBS -l walltime=10:00:00
#PBS -k oe

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 16 FreeFem++-mpi -v 0 simple5.edp -isSplit 1 -iter 0 -isUniRef 1 -hsize 0.01 -nev 64 -N 12
#mpirun -np 16 FreeFem++-mpi -v 0 solveBEDMAP2.edp -nborders 6 -isMesh 1 -isSplit 1 -isUniRef 1 -iter 0 -hsize 0.05 -nev 64
exit 0
