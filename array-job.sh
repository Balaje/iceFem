#!/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=8:mem=16GB
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -J 0-49
#PBS -q testq

JOBS=50
TMIN=$(echo 1/0.0001 | bc -l)
TMAX=$(echo 1/0.01 | bc -l)

STEP=$(echo \($TMAX-$TMIN\)/\($JOBS-1\) | bc -l)
TVAL=$(echo $TMIN+$PBS_ARRAY_INDEX*$STEP | bc -l)

source /etc/profile.d/modules.sh
module load freefem
cd $PBS_O_WORKDIR
export FF_INCLUDEPATH="$PBS_O_WORKDIR/include"
PbsJobID==$(echo $PBS_JOBID | cut -d "[" -f1) # Extract the job-id

# Function to check if all jobs are over.
checkDoneScript(){    
    declare -a DONE_ARRAY=( $(for i in `seq 0 $(( $njobs - 1 ))`; do echo ${DONE_ARRAY[i]:-bar}; done ) );
    N_JOBS_DONE=0;
    while [ true ]; do
	for i in `seq 0 $(( $JOBS -1 ))`; do
	    FILENAME="array-job.sh.o$PbsJobID.$i"
	    if [ -f "$FILENAME" ]; then # If file exists ...
		tag=$(tail -n 1 $FILENAME) # Get the last line
		if [ "$tag" == "**done**" ] && [[ "${DONE_ARRAY[i]}" -eq "bar" ]]; then
		    echo "Program $i done"
		    DONE_ARRAY[i]=1;
		    N_JOBS_DONE=$(( $N_JOBS_DONE + 1 ))
		fi	    
	    fi	
	done
	
	sleep 5
	if [ $N_JOBS_DONE -eq $JOBS ]; then
	    break;
	fi    
    done
}

echo "Run FreeFem scripts ..."
echo mpirun -np 8 FreeFem++-mpi -nw -v 0 solveBEDMAP2.edp -isMesh 0 -nborders 6 -Tr $TVAL -Ti 0 -iter $(( $PBS_ARRAY_INDEX+1 ))
echo $PBS_JOBID
echo "**done**"

# Check if all the jobs are done.
checkDoneScript;

echo "All FreeFem Scripts Finished ... Running MATLAB Code"
echo "**done**"

# Check if all the jobs are done.
checkDoneScript;

echo "All MATLAB Scripts Finished ... Running Time domain Code"
echo "**done**"

exit 0
