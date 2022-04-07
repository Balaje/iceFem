#############################################################################################
# Python script to compute the analytic extension of the solution
# Plots the reflection and transmission coefficients in the complex plane
# as a function of frequency.
#
# Run:
#           python3 ComplexRefIce.py ICEBERG_COMPLEX/
#
# Build the initial complex frequency grid.
# Externally call the FreeFem script iceshelf2d.edp to solve the problem on the 11x11 grid.
# Interpolate the solution to a 300x300 grid.
# Compute the reflection and transmission coefficients on the complex plane.
# Plot the solution.
#
#
# Running notes:
#         Before running, obtain the paths for mpirun and FreeFem++-mpi
#         Open terminal and run
#                   which mpirun
#                   which FreeFem++-mpi
#         Append the results into the cmd variable inside the for loop. The same can be done
#         using a bash script and then ignoring the for-loop part in this code completely.
##############################################################################################

import numpy as np
from modules.interpolateFreq import *
import cplot
import matplotlib.pyplot as plt
from os import system
from matplotlib.ticker import FormatStrFormatter
import sys

SolutionDir=sys.argv[1];

npts=11
pi=np.pi

a=2*pi/100
b=2*pi/8
c=-0.08
d=0.08

x=np.linspace(a,b,npts)
y=np.linspace(c,d,npts)
X,Y=np.meshgrid(x,y)
omega=X+1j*Y

count=1
for m in range(0,npts):
    for n in range(0,npts):
        cmd=('/opt/homebrew/bin/mpirun -np 2 /usr/local/ff++/mpich3/bin/FreeFem++-mpi '
             '-v 0 iceshelf2d.edp -N 12 -omegar ' +str(omega[m,n].real)+
             ' -omegai '+str(omega[m,n].imag)+
             ' -Youngs 2e9 -L 3630 -H 500 -h 280 -nev 10 -iter ' +str(count)+
             ' -solDir '+SolutionDir+' -hsize 0.01 -isSplit 1 > /dev/null')
        #print(cmd)
        a=system(cmd)
        print("Done running",count)
        count=count+1

NModes=3
nptsNew=200
xq=np.linspace(a,b,nptsNew)
yq=np.linspace(c,d,nptsNew)
Xq,Yq=np.meshgrid(xq,yq)
omeganew=Xq+1j*Yq;

interpolateCoeffsFreqComplex(a,b,c,d,npts,10,SolutionDir+"2_ModesMatrix/",nptsNew)
LAM=buildLam(SolutionDir)

interpolateRefCoeffComplex(a,b,c,d,npts,10,SolutionDir+"2_RefCoeff/","C",nptsNew,NModes)
RC=buildRMat(LAM,SolutionDir,"C",0)

plt.figure(figsize=[10,4])
ax=plt.subplot(2,1,1)
genComplexPlot(RC,omeganew/(2*pi))

interpolateRefCoeffComplex(a,b,c,d,npts,10,SolutionDir+"2_RefCoeff/","T",nptsNew,NModes)
RT=buildRMat(LAM,SolutionDir,"T",0)
ax=plt.subplot(2,1,2)
genComplexPlot(RT,omeganew/(2*pi))

plt.show()
