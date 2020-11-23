#!/bin/bash

#PBS -l select=1:ncpus=14:mpiprocs=14:mem=16GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 2-14

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np $PBS_ARRAY_INDEX FreeFem++-mpi -v 0 simple5.edp -isSplit 1 -iter $(( $PBS_ARRAY_INDEX+1 )) -isUniRef 1 -hsize 0.005 -nev 64 -N 12
#mpirun -np $PBS_ARRAY_INDEX FreeFem++-mpi -v 0 solveBEDMAP2.edp -isSplit 1 -iter $(( $PBS_ARRAY_INDEX+1 )) -hsize 0.01 -nborders 6 -isMesh 0 -nev 32
exit 0
