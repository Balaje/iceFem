import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys

filePath=sys.argv[1];

pi=np.pi;
omega=np.linspace(2*pi*0.0001,2*pi*0.01,10)
T=2*pi/omega;

npts=199
nev=16
omeganew=interpolateCoeffsFreq(2*pi*0.0001,2*pi*0.01,omega,nev,filePath+"2_ModesMatrix/",npts,1)
LAM=buildLam(filePath)

plt.subplot(2,1,1)
plt.title("Coefficients vs Frequency")
for m in range(0,5):
    L=LAM[:,m]
    plt.plot(omeganew,abs(L),linewidth=2)

#plt.show()
