#!/bin/bash

#PBS -l select=1:ncpus=2:mpiprocs=2:mem=16GB
#PBS -l walltime=10:00:00
#PBS -k oe

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 2 FreeFem++-mpi -v 0 eigenSolve.edp -hsize 0.01 -N 12 -nev 64 -isUniRef 1
#mpirun -np 2 FreeFem++-mpi -v 0 eigenSolve.edp -hsize 0.01 -N 12 -nev 64 -isBEDMAP 1 -isMesh 1

