import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys

plt.figure(figsize=[10,6])
for i in range(1,len(sys.argv)):
    filePath=sys.argv[i];
    pi=np.pi;
    T=np.linspace(20,1160,25)
    omega=2*pi/T
    npts=399
    nev=16
    omeganew=interpolateCoeffsFreq(2*pi/1160,2*pi/20,omega,nev,filePath+"2_ModesMatrix/",npts,1)
    LAM=buildLam(filePath)
    ax=plt.subplot(len(sys.argv)-1,1,i)
    heading="Coefficients vs Frequency ("+filePath+")";
    ax.text(.5,.9,heading,
        horizontalalignment='center',
        transform=ax.transAxes)
    for m in range(0,4):
        L=LAM[:,m]
        plt.semilogy(omeganew,abs(L),linewidth=2,label='Mode '+str(m+1))

    plt.xlabel('$\omega$')
    plt.legend(loc='upper right')

plt.show()
