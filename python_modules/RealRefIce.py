#####################################################################
# Python script to compute and plot the reflection coefficients.
#
# Run:
#        python3 RealRefIce.py ICEBERG1/
#
# The initial frequency space: omega
# Interpolated onto a finer space omeganew with endpoints of omega
# Compute the modal coefficients and plot
# Compute the reflection and transmission coefficient and plot
# Compute the energy conservation result.
######################################################################

import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys
import time

filePath=sys.argv[1];

pi=np.pi
omega=2*pi*np.linspace(0.01,0.125,51)
T=2*pi/omega
RC=np.zeros((51,4),dtype=complex)

npts=399
nev=10
omeganew=interpolateCoeffsFreq(2*pi*0.01,2*pi*0.125,omega,nev,filePath+"2_ModesMatrix/",npts,1)
LAM=buildLam(filePath)

plt.figure(figsize=[10,4])
plt.title("Modal Amplitudes")
for m in np.arange(3,6):
    L=LAM[:,m]
    plt.plot(omeganew/(2*pi),abs(L),linewidth=2,label="$|\lambda_"+str(m+1)+"|$")
plt.legend()
plt.xlabel('$\omega/(2\pi)$')
plt.ylim([0,5])

## Interpolating reflection coefficients
NModes = 3
V=interpolateRefCoeff(omega,omeganew,nev,filePath+"2_RefCoeff/",NModes,"C")
RC=buildRMat(LAM,filePath,"C",0)
plt.figure(figsize=[10,4])
plt.plot(omeganew/(2*pi),np.transpose(abs(RC)),linewidth=2,label="R($\omega$)")

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,nev,filePath+"2_RefCoeff/",NModes,"T")
RT=buildRMat(LAM,filePath,"T",0)
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)),linewidth=2,label="T($\omega$)")

# Check Energy Conservation
plt.plot(omeganew/(2*pi),np.transpose(abs(RT)**2+abs(RC)**2),linewidth=2,label="$1.0$")

plt.title("$T(\omega)$ and $R(\omega)$")
plt.legend()
plt.xlabel('$\omega/(2\pi)$')
plt.ylim([0,1.1])

plt.show()
