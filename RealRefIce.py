import numpy as np
from modules.interpolateFreq import *
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
LAM=buildLam(filePath)

plt.subplot(2,1,1)
for m in range(0,4):
    L=LAM[:,m]
    plt.plot(omeganew,abs(L),linewidth=2)

## Interpolating reflection coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","C")
RC=buildRMat(LAM,filePath,"C")
plt.subplot(2,1,2)
plt.plot(omeganew,np.transpose(abs(RC)),linewidth=2)

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","T")
RT=buildRMat(LAM,filePath,"T")
plt.plot(omeganew,np.transpose(abs(RT)),linewidth=2)

# Check Energy Conservation
plt.plot(omeganew,np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2)
plt.show()
