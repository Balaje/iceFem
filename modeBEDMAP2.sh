#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 1000-1007

IND=$(( $PBS_ARRAY_INDEX-999 ))
TVAL=$(sed "${IND}q;d" waveperiods.txt)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 6 -Tr $TVAL -Ti 0 -iter $(( $PBS_ARRAY_INDEX+1 ))
exit 0


#iter=1000
#for line in `cat waveperiods.txt`
#do    
#    mpirun -np 4 FreeFem++-mpi -v 0 solveBEDMAP2.edp -Tr $line -Ti 0 -isMesh 0 -nborders 6 -iter $iter > $iter.txt
#    iter=$(( $iter + 1))
#done

   
