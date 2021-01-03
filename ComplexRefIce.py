import numpy as np
from modules.interpolateFreq import interpolateCoeffsFreqComplex, interpolateRefCoeffComplex
import cplot
import matplotlib.pyplot as plt
from os import system

SolutionDir="1_ICEBERG/"

npts=10;
pi=np.pi;

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
        a=system(cmd)
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
LAMRe=np.loadtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
LAMIm=np.loadtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
LAM=LAMRe+1j*LAMIm

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","C",nptsNew)
RefCDifRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refCDifRe.dat");
RefCDifIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refCDifIm.dat");
RefCRadRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refCRadRe.dat");
RefCRadIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refCRadIm.dat");
RCDif=RefCDifRe+1j*RefCDifIm
RCRad=RefCRadRe+1j*RefCRadIm
RC=np.zeros(len(RCDif),dtype=complex)
for indx in range(0,len(RC)):
    L=LAM[indx,:]
    RC[indx]=RCDif[indx]+np.dot(RCRad[indx,:],L)

RCNew=np.reshape(RC,np.shape(omeganew))
VAL=cplot.get_srgb1(RCNew)
plt.subplot(1,2,1)
plt.imshow(VAL)

interpolateRefCoeffComplex(a,b,c,d,npts,8,SolutionDir+"2_RefCoeff/","T",nptsNew)
RefTDifRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refTDifRe.dat");
RefTDifIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refTDifIm.dat");
RefTRadRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refTRadRe.dat");
RefTRadIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/refTRadIm.dat");
RTDif=RefTDifRe+1j*RefTDifIm
RTRad=RefTRadRe+1j*RefTRadIm
RT=np.zeros(len(RTDif),dtype=complex)
for indx in range(0,len(RC)):
    L=LAM[indx,:]
    RT[indx]=RTDif[indx]+np.dot(RTRad[indx,:],L)

RTNew=np.reshape(RT,np.shape(omeganew))
VAL=cplot.get_srgb1(RTNew)
plt.subplot(1,2,2)
plt.imshow(VAL)
plt.show()
