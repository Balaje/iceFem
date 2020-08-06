#!/bin/bash

#PBS -l select=1:ncpus=10:mpiprocs=10:mem=32GB
#PBS -l walltime=10:00:00
#PBS -k oe

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
mpirun -np 10 FreeFem++-mpi -nw -v 0 timeDomain.edp
exit 0
