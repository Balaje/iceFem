#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 0-49

JOBS=50
TMIN=$(echo 1/0.0001 | bc -l)
TMAX=$(echo 1/0.01 | bc -l)

STEP=$(echo \($TMAX-$TMIN\)/\($JOBS-1\) | bc -l)
TVAL=$(echo $TMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 6 -Tr $TVAL -Ti 0 -iter $(( $PBS_ARRAY_INDEX+1 ))
exit 0
