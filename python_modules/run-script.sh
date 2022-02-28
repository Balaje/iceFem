#!/bin/zsh

##################################################################################
# Script to solve the iceberg vibration problem using iceFEM.
#
# INSTRUCTIONS:
# - The numerical model is a modification of the model used to simulate the
# wave-induced motion of iceshelf. Make sure to uncomment the line:
#       iceshelf2iceberg;
# in the script.
#
# - This script solves the frequency domain problem on the real axis for
# two values of Young's moduli.
#
# - Then the program getSolutionInterp.edp interpolates the frequency domain
# on the finer space and computes the strain and the displacement as a function
# of frequency. Written in
#           $solDir/strainVsFreq.dat
#           $solDir/dispVsFreq.dat
#
# - Then use the Python Scripts to perform the post-processing.
#
# NOTE:
# A similar script can be written for the problem on the complex frequency domain
# See ComplexRefIce.py for the direct script. However, the getSolutionInterp.edp
# program does not work since FreeFem does not support 2D interpolation, and
# cannot be used to obtain the strainVsFreq.dat, dispVsFreq.dat files.
##################################################################################

for y in $(seq 1 0.5 1.5);
do
    youngs=$(echo $y\*10^9 | bc -l)
    solDir=ICEBERG$y/
    echo -n "Solving "
    for i in $(seq 1 51);
    do
        omegar=$(echo 2\*a\(1.\)\*4\*\(0.01 + \($i-1\)\*0.0023\) | bc -l);
        echo -n "$i..."
        mpirun -np 2 FreeFem++-mpi -v 0 iceshelf2d.edp -L 3630 -H 500 -h 280 \
               -Youngs $youngs -nev 10 -hsize 0.01 -N 12 -isUniRef 1 -isSplit 1 \
               -NModes 3 -iter $i -omegar $omegar -solDir $solDir > /dev/null;
    done
    omega0=$(echo 2\*a\(1.\)\*4\*0.01 | bc -l)
    omega1=$(echo 2\*a\(1.\)\*4\*0.125 | bc -l)
    mpirun -np 2 FreeFem++-mpi -v 0 getSolutionInterp.edp -L 3630 -H 500 -h 280 \
           -Youngs $youngs -nev 10 -hsize 0.01 -N 12 -iter test -isUniRef 0 \
           -omega0 $omega0 -omega1 $omega1 -nfreqcoarse 51 -nfreqfine 400 \
           -startindx 1 -endindx 51 -SolDir $solDir

    echo "Done Youngs = $youngs. Results in $solDir"
done
