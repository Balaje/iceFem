#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=32GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 0-50

JOBS=50
TMIN=1000
TMAX=10000

STEP=$(echo \($TMAX-$TMIN\)/$JOBS | bc -l)
TVAL=$(echo $TMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 6 -Tr $TVAL -Ti 0 -iter $PBS_ARRAY_INDEX
exit 0
