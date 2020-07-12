echo "#!/bin/bash" > job.sh
echo "" >> job.sh
echo "#PBS -l select=$1:ncpus=$2:mem=16GB" >> job.sh
echo "#PBS -l walltime=1:00:00" >> job.sh
echo "#PBS -k oe" >> job.sh
echo "" >> job.sh
echo "source /etc/profile.d/modules.sh" >> job.sh
echo "module load freefem" >> job.sh
echo "cd \$PBS_O_WORKDIR" >> job.sh
echo "export FF_INCLUDEPATH=\"\$PBS_O_WORKDIR/include\"" >> job.sh
var=$(echo "$1*$2" | bc -l)
echo "mpirun -np $var FreeFem++-mpi -nw -v 0 $3 -L $4 -H $5 -h $6 -Tr $7 -Ti $8 -iter $9" >> job.sh
echo "exit 0" >> job.sh

cat job.sh
