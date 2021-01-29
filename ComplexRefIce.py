import numpy as np
from modules.interpolateFreq import *
import cplot
import matplotlib.pyplot as plt
from os import system
from matplotlib.ticker import FormatStrFormatter
import sys

SolutionDir=sys.argv[1];

npts=11;
pi=np.pi

a=2*pi/80
b=2*pi/10
c=-0.08
d=0.08

x=np.linspace(a,b,npts)
y=np.linspace(c,d,npts)
X,Y=np.meshgrid(x,y)
omega=X+1j*Y
T=2*pi/omega
RC=np.zeros((npts**2,4))
count=1
for m in range(0,npts):
    for n in range(0,npts):
        cmd="/usr/local/bin/mpirun -np 2 /usr/local/ff++/mpich3/bin/FreeFem++-mpi -v 0 iceberg.edp -N 10 -Tr "+str(T[m,n].real)+" -Ti "+str(T[m,n].imag)+" -Youngs 2e8 -L 3000 -H 1000 -h 200 -nev 8 -iter "+str(count)+" -SolDir 3_ICEBERG/ -hsize 0.04 -isSplit 1 > /dev/null"
        #print(cmd)
        #a=system(cmd)
        print("Done running",count)
        count=count+1
for m in range(0,npts**2):
    RC[m,:]=np.loadtxt(SolutionDir+"2_RefCoeff/rc"+str(m+1)+".dat")
RCOld=RC[:,0]+1j*RC[:,1]
RCOld=np.reshape(RCOld,(npts,npts))
TCOld=RC[:,2]+1j*RC[:,3]
TCOld=np.reshape(TCOld,(npts,npts))

nptsNew=300;
xq=np.linspace(2*pi/80,2*pi/10,nptsNew)
yq=np.linspace(-0.08,0.08,nptsNew)
Xq,Yq=np.meshgrid(xq,yq)
omeganew=Xq+1j*Yq;
Tnew=2*pi/omeganew

interpolateCoeffsFreqComplex(a,b,c,d,npts,8,SolutionDir+"2_ModesMatrix/",nptsNew)
LAM=buildLam(SolutionDir)

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","C",nptsNew)
RC=buildRMat(LAM,SolutionDir,"C")

plt.figure(figsize=[10,4])
ax=plt.subplot(2,1,1)
genComplexPlot(RC,omeganew/(2*pi))

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","T",nptsNew)
RT=buildRMat(LAM,SolutionDir,"T")
ax=plt.subplot(2,1,2)
genComplexPlot(RT,omeganew/(2*pi))

plt.savefig("ComplexPlot.pdf",bbox_inches='tight')
