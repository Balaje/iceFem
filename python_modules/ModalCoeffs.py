##########################################################################
# Python Script to compute the modal coefficients \lambda_j:
#             [H]{\lambdaj} = {f}
# Run:
#        python3 ModalCoeffs.py ICEBERG1/
#        python3 ModalCoeffs.py ICEBERG1/ ICEBERG1.5/
#
# The initial frequency space: omega
# Interpolated onto a finer space omeganew with the endpoints of omega.
# Compute the modal coefficients and plot
##########################################################################

import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys

pi=np.pi
omega = 2*pi*np.linspace(0.01,0.125,51)

plt.figure(figsize=[10,6])
for i in range(1,len(sys.argv)):
    filePath=sys.argv[i];
    npts=399
    nev=10
    omeganew=interpolateCoeffsFreq(2*pi*0.01, 2*pi*0.125, omega, nev, filePath+"2_ModesMatrix/",npts,1)
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
