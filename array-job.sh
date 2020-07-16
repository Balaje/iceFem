#!/bin/bash

#PBS -l select=2:ncpus=10:mpiprocs=10:mem=32GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 0-20

JOBS=20
TMIN=20
TMAX=200

STEP=$(echo \($TMAX-$TMIN\)/$JOBS | bc -l)
TVAL=$(echo $TMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 20 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 4 -Tr $TVAL -Ti 0 -iter $PBS_ARRAY_INDEX
exit 0
