import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys

filePath=sys.argv[1];

T=np.linspace(10,80,21)
pi=np.pi
omega=2*pi/T
RC=np.zeros((21,4),dtype=complex)

npts=299
nev=8
omeganew=interpolateCoeffsFreq(2*pi/80,2*pi/10,omega,nev,filePath+"2_ModesMatrix/",npts,1)
LAM=buildLam(filePath)

plt.figure(figsize=[10,4])
plt.title("Modal Amplitudes")
for m in np.arange(0,6):
    L=LAM[:,m]
    plt.plot(omeganew/(2*pi),abs(L),linewidth=2,label="$|\lambda_"+str(m+1)+"|$")
plt.legend()
plt.xlabel('$\omega$')
plt.savefig(sys.argv[2]+"_1.pdf",bbox_inches='tight')

## Interpolating reflection coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","C")
RC=buildRMat(LAM,filePath,"C")
plt.figure(figsize=[10,4])
plt.plot(omeganew/(2*pi),np.transpose(abs(RC)),linewidth=2,label="R($\omega$)")

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,8,filePath+"2_RefCoeff/","T")
RT=buildRMat(LAM,filePath,"T")
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)),linewidth=2,label="T($\omega$)")
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2,label="$1.0$")
plt.title("$T(\omega)$ and $R(\omega)$")
plt.legend()
plt.xlabel('$\omega$')
plt.savefig(sys.argv[2]+"_2.pdf",bbox_inches='tight')

# Check Energy Conservation
#plt.xlim([1/80,0.2])
plt.ylim([0,1.1])
print(filePath)
