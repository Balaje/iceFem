import numpy as np
from modules.interpolateFreq import interpolateCoeffsFreq, interpolateRefCoeff
import matplotlib.pyplot as plt

filePath="2_ICEBERG/"

T=np.linspace(15,80,27)
pi=np.pi
omega=2*pi/T
RC=np.zeros((27,4),dtype=complex)
for m in range(1,27):
    RC[m,:]=np.loadtxt(filePath+"2_RefCoeff/rc"+str(m)+".dat")
plt.plot(omega,abs(RC[:,0]+1j*RC[:,1]),'r+')

npts=299
omeganew=interpolateCoeffsFreq(2*pi/80,2*pi/15,omega,8,filePath+"2_ModesMatrix/",npts,1)

LAMRe=np.loadtxt(filePath+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
LAMIm=np.loadtxt(filePath+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
LAM=LAMRe+1j*LAMIm

plt.subplot(2,1,1)
for m in range(0,4):
    L=LAM[:,m]
    plt.plot(omeganew,abs(L),linewidth=2)

## Interpolating reflection coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","C")
RefCDifRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCDifRe.dat");
RefCDifIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCDifIm.dat");
RefCRadRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCRadRe.dat");
RefCRadIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refCRadIm.dat");
RCDif=RefCDifRe+1j*RefCDifIm
RCRad=RefCRadRe+1j*RefCRadIm
RC=np.zeros(len(RCDif),dtype=complex)

for indx in range(0,len(RC)):
    L=LAM[indx,:]
    RC[indx]=RCDif[indx]+np.dot(RCRad[indx,:],L)
plt.subplot(2,1,2)
plt.plot(omeganew,np.transpose(abs(RC)),linewidth=2)

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","T")
RefTDifRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTDifRe.dat");
RefTDifIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTDifIm.dat");
RefTRadRe=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTRadRe.dat");
RefTRadIm=np.loadtxt(filePath+"2_RefCoeff/Interpolated_R/refTRadIm.dat");
RTDif=RefTDifRe+1j*RefTDifIm
RTRad=RefTRadRe+1j*RefTRadIm
RT=np.zeros(len(RTDif),dtype=complex)

for indx in range(0,len(RT)):
    L=LAM[indx,:]
    RT[indx]=RTDif[indx]+np.dot(RTRad[indx,:],L)
plt.plot(omeganew,np.transpose(abs(RT)),linewidth=2)
plt.plot(omeganew,np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2)
plt.show()
