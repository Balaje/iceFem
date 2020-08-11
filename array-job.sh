#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 0-199

JOBS=200

OMEGAMIN=0.0005
OMEGAMAX=0.01
PI=$(echo 4*a\(1\) | bc -l)

STEP=$(echo \($OMEGAMAX-$OMEGAMIN\)/\($JOBS-1\) | bc -l)
OMEGAVAL=$(echo $OMEGAMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)
TVAL=$(echo 2*$PI/$OMEGAVAL | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 6 -Tr $TVAL -Ti 0 -iter $(( $PBS_ARRAY_INDEX+1 ))
exit 0
