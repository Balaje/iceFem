import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys

filePath=sys.argv[1];

T=np.linspace(10,80,27)
pi=np.pi
omega=2*pi/T
RC=np.zeros((27,4),dtype=complex)
for m in range(1,10):
    RC[m,:]=np.loadtxt(filePath+"2_RefCoeff/rc"+str(m)+".dat")
plt.figure(figsize=(15,8))
plt.plot(omega/(2*pi),abs(RC[:,0]+1j*RC[:,1]),'r+')

npts=599
omeganew=interpolateCoeffsFreq(2*pi/80,2*pi/10,omega,8,filePath+"2_ModesMatrix/",npts,1)
LAM=buildLam(filePath)

plt.subplot(2,1,1)
plt.title("Coefficients vs Frequency")
for m in np.arange(0,6):
    L=LAM[:,m]
    plt.plot(omeganew/(2*pi),abs(L),linewidth=2,label="$|\lambda_"+str(m+1)+"|$")
plt.legend()
plt.xlim([1/80,0.2])
plt.ylim([0,4])

## Interpolating reflection coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","C")
RC=buildRMat(LAM,filePath,"C")
plt.subplot(2,1,2)
plt.plot(omeganew/(2*pi),np.transpose(abs(RC)),linewidth=2)

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","T")
RT=buildRMat(LAM,filePath,"T")
plt.title("Transmission and Reflection Coefficients vs Frequency")
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)),linewidth=2)

# Check Energy Conservation
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2)
plt.xlim([1/80,0.2])
plt.ylim([0,1.1])
print(filePath)
plt.savefig(sys.argv[2])
