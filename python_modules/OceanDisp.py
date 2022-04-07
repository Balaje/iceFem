#################################################################################################
# Python script to compute the displacement of the ocean at x=-x0, x=L+x0
# This yields the wave heights of the incident and transmitted wave
# Run:
#         python3 OceanDisp.py ICEBERG1/ 1e9
#         python3 OceanDisp.py ICEBERG1.5/ 1.5e9
# where 1.5e9 is the Young's modulus
#
# The initial frequency space omega
# The dimensional parameters L, H, h to compute the roots of the dispersion equation
# Interpolate the solution to a finer space.
# Obtain the reflection coefficients - This is just the solution vector in the open-ocean part.
# Compute the displacement and plot.
#################################################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
from modules.interpolateFreq import *
from modules.nonDimensional import *
from modules.dispersionFreeSurface import *

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

SolutionDir = sys.argv[1]

# Dimensional parameters of the ice and cavity
L = 3630
H = 500
h = 280
Youngs = float(sys.argv[2])
pi = np.pi
omega = 2*pi*np.linspace(0.01,0.125,51)

NModes = 3
omegaNew = interpolateCoeffsFreq(2*pi*0.01, 2*pi*0.125, omega, 10, SolutionDir+"2_ModesMatrix/", 299, 1)
alphaNew = np.zeros(np.size(omegaNew), dtype=complex)

LAMRe = np.genfromtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
LAMIm = np.genfromtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
LAM = LAMRe+1j*LAMIm
interpolateRefCoeff(omega, omegaNew, np.size(LAM,1), SolutionDir+"2_RefCoeff/", NModes, "T")
interpolateRefCoeff(omega, omegaNew, np.size(LAM,1), SolutionDir+"2_RefCoeff/", NModes, "C")

Cp = np.zeros((np.size(LAM,0),NModes+1), dtype=complex)
Cm = np.zeros((np.size(LAM,0),NModes+1), dtype=complex)
k = np.zeros((np.size(LAM,0),NModes+1), dtype=complex)
for m in np.arange(0,len(alphaNew)):
    [LL,HH,hh,alphaNew[m]] = nonDimensionalize(L, H, h, Youngs, omegaNew[m])
    k[m,:] = dispersion_free_surface(alphaNew[m], NModes, HH).ravel()

RC = buildRMat(LAM, SolutionDir, "C", 0)
RT = buildRMat(LAM, SolutionDir, "T", 0)

#plt.figure(figsize=[10,4])
fig, ax = plt.subplots(2,1,figsize=[10,6])
fig1, ax1 = plt.subplots(2,1,figsize=[10,6])
for m in np.arange(0,NModes+1):
    Cp[:,m] = buildRMat(LAM, SolutionDir, "C", m)
    Cm[:,m] = buildRMat(LAM, SolutionDir, "T", m)

    ax[0].semilogy(omegaNew/(2*pi), abs(Cm[:,m]), label="Mode = "+str(m))
    ax1[0].semilogy(omegaNew/(2*pi), abs(Cp[:,m]), label="Mode = "+str(m))

for x in [0, 0.01*LL, 0.05*LL, 0.1*LL]:
    sumR = np.exp(k[:,0]*x)
    sumT = 0*RT
    for m in np.arange(0,NModes+1):
        sumT = sumT + Cm[:,m]*np.exp(k[:,m]*(LL-(x+LL)))
        sumR = sumR + Cp[:,m]*np.exp(k[:,m]*(-x))

    ax[1].semilogy(omegaNew/(2*pi), abs(sumT), label="$x = L + $"+str(x/LL)+"$L$")
    ax1[1].semilogy(omegaNew/(2*pi), abs(sumR), label="$x = L + $"+str(x/LL)+"$L$")

ax[0].set_ylabel("$T(\omega)$")
ax[1].set_ylabel("$||\eta(\omega)||_\infty$")
ax1[0].set_ylabel("$R(\omega)$")
ax1[1].set_ylabel("$||\eta(\omega)||_\infty$")

ax[1].set_xlabel("$\omega/(2\pi)$")
ax1[1].set_xlabel("$\omega/(2\pi)$")

ax[0].legend()
ax1[1].legend()
ax[1].legend()
ax1[0].legend()

plt.show()
