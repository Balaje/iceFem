#!/bin/bash

#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=1:00:00
#PBS -k oe
#PBS -J 0-10

JOBS=10
TMIN=20
TMAX=200

STEP=$(echo \($TMAX-$TMIN\)/$JOBS | bc -l)
TVAL=$(echo $TMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
echo mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 4 -Tr $TVAL -i 0 -iter $PBS_ARRAY_INDEX
exit 0
