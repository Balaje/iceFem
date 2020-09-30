#!/bin/bash

#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=10:00:00
#PBS -k oe
#PBS -J 0-14

JOBS=15

NWMIN=0.05
NWDELTA=0.05
NW=$(echo \($NWMIN+$PBS_ARRAY_INDEX*$NWDELTA\) | bc -l)


source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"

echo mpirun -np 8 FreeFem++-mpi -nw -v 0 simple5.edp -isSplit 1 -Tr 2000 -iter $(( $PBS_ARRAY_INDEX+1 )) -notchWidth 0.005 -notchHeight $NW

exit 0
