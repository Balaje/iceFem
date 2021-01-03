import numpy as np
from modules.interpolateFreq import interpolateCoeffsFreq, interpolateRefCoeff
import matplotlib.pyplot as plt

filePath="2_ICEBERG/"

T=np.linspace(15,80,27)
pi=np.pi
omega=2*pi/T
npts=299
omeganew=interpolateCoeffsFreq(2*pi/80,2*pi/15,omega,8,filePath+"2_ModesMatrix/",npts,1)

LAMRe=np.loadtxt(filePath+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
LAMIm=np.loadtxt(filePath+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
LAM=LAMRe+1j*LAMIm

#for m in range(0,4):
#    L=LAM[:,m]
#    plt.plot(omeganew,abs(L),linewidth=2)
    #plt.show()

## Interpolating reflection coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","C")
RefCDifRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCDifRe.dat");
RefCDifIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCDifIm.dat");
RefCRadRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCRadRe.dat");
RefCRadIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCRadIm.dat");
RCDif=RefCDifRe+1j*RefCDifIm
RCRad=RefCRadRe+1j*RefCRadIm
RC=np.zeros((1,len(RCDif)),dtype=complex)

for indx in range(0,npts+1):
    L=LAM[indx,:]
    RC[0,indx]=RCDif[indx]+np.dot(L,RCRad[indx,:])
plt.plot(omeganew,np.transpose(abs(RC)),linewidth=2)

V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","T")
RefTDifRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTDifRe.dat");
RefTDifIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTDifIm.dat");
RefTRadRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTRadRe.dat");
RefTRadIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTRadIm.dat");
RTDif=RefTDifRe+1j*RefTDifIm
RTRad=RefTRadRe+1j*RefTRadIm
RT=np.zeros((1,len(RTDif)),dtype=complex)

for indx in range(0,npts+1):
    L=LAM[indx,:]
    RT[0,indx]=RTDif[indx]+np.dot(L,RTRad[indx,:])
plt.plot(omeganew,np.transpose(abs(RT)),linewidth=2)
plt.plot(omeganew,np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2)
plt.show()
