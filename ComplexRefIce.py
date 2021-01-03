import numpy as np
from modules.interpolateFreq import *
import cplot
import matplotlib.pyplot as plt
from os import system

SolutionDir="1_ICEBERG/"

npts=10;
pi=np.pi

a=2*pi/80
b=2*pi/15
c=-0.06
d=0.06

x=np.linspace(a,b,npts)
y=np.linspace(c,d,npts)
X,Y=np.meshgrid(x,y)
omega=X+1j*Y
T=2*pi/omega
RC=np.zeros((npts**2,4))
count=1
for m in range(0,npts):
    for n in range(0,npts):
        cmd="/usr/local/bin/mpirun -np 2 /usr/local/ff++/mpich3/bin/FreeFem++-mpi -v 0 iceberg.edp -N1 20 -N2 80 -Tr "+str(T[m,n].real)+" -Ti "+str(T[m,n].imag)+" -L 1000 -H 5000 -h 200 -nev 8 -iter "+str(count)+" > /dev/null";
        #print(cmd)
        #a=system(cmd)
        #print("Done running",count)
        count=count+1
for m in range(0,npts**2):
    RC[m,:]=np.loadtxt(SolutionDir+"2_RefCoeff/rc"+str(m+1)+".dat")
RCOld=RC[:,0]+1j*RC[:,1]
RCOld=np.reshape(RCOld,(npts,npts))
TCOld=RC[:,2]+1j*RC[:,3]
TCOld=np.reshape(TCOld,(npts,npts))

nptsNew=100;
xq=np.linspace(2*pi/80,2*pi/15,nptsNew)
yq=np.linspace(-0.06,0.06,nptsNew)
Xq,Yq=np.meshgrid(xq,yq)
omeganew=Xq+1j*Yq;
Tnew=2*pi/omeganew

interpolateCoeffsFreqComplex(a,b,c,d,npts,8,SolutionDir+"2_ModesMatrix/",nptsNew)
LAM=buildLam(SolutionDir)

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","C",nptsNew)
RC=buildRMat(LAM,SolutionDir,"C")
RCNew=np.reshape(RC,np.shape(omeganew))
RCNew=RCNew.T
VAL=cplot.get_srgb1(RCNew,colorspace="hsl",alpha=0)
plt.subplot(2,1,1)
plt.title("Reflection Coefficient")
plt.imshow(VAL,aspect=0.4)

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","T",nptsNew)
RT=buildRMat(LAM,SolutionDir,"T")
RTNew=np.reshape(RT,np.shape(omeganew))
RTNew=RTNew.T
VAL=cplot.get_srgb1(RTNew,colorspace="hsl",alpha=0)
plt.subplot(2,1,2)
plt.title("Transmission Coefficient")
plt.imshow(VAL,aspect=0.4)
plt.show()
