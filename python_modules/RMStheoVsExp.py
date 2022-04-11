#######################################################
# Python script to plot the RMS strain vs omega
# Run
#           python3 RMStheoVsExp.py [DIRECTORY]
#
# Before running the script, solve a set of frequency
# domain problems for a single Young's modulus.
#
# Compute the RMS strain and compare with the experimental
# data:
#      Kristensen, M., Squire, V.A., Moore, S.C., 1982. Tabular icebergs in ocean waves 297, Nature 3.
#######################################################

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import sys

def rolling_rms(x, y, N):
    yc = np.cumsum(abs(y));
    xc = np.cumsum(abs(x));
    return np.sqrt((yc[N:] - yc[:-N])*(xc[N:] - xc[:-N]) / N)

fig = plt.figure(figsize=[10,8])
g=9.8
omegap=0.092*2*np.pi
omega0=1.14*omegap
alpha=0.0081
beta=0.74

num_samples = 5
for i in range(1, len(sys.argv)):
# Load the simulation data and compute strain^2*PSD
    filePath=sys.argv[i]
    strainVsFreq = np.genfromtxt(filePath+"/strainVsFreq.dat")
    omeganew = 2*np.pi*np.linspace(0.01,0.125,len(strainVsFreq))
    S=alpha*g**2/omeganew**5*np.exp(-beta*(omega0/omeganew)**4) + 1e-10;
    strain2 = (strainVsFreq[:,0]/1e-6)**2
    strain2S = (strain2)*S
    strain2Srms = rolling_rms(omeganew/(2*np.pi), strain2S, num_samples)
    if(i==1):
        ax1 = plt.subplot(2,1,1)
        plt.plot(omeganew/(2*np.pi), strain2S, 'r-', linewidth=2, label="E = "+sys.argv[i][7:]+" GPa $(E_{opt})$")
        ax2 = plt.subplot(2,1,2)
        plt.plot(omeganew[num_samples::]/(2*np.pi), strain2Srms, 'r-', linewidth=2, label="E = "+sys.argv[i][7:]+" GPa $(E_{opt})$")
    else:
        ax1 = plt.subplot(2,1,1)
        plt.plot(omeganew/(2*np.pi), strain2S, label="E = "+sys.argv[i][7:]+" GPa", linewidth=1)
        ax2 = plt.subplot(2,1,2)
        plt.plot(omeganew[num_samples::]/(2*np.pi), strain2Srms, label="E = "+sys.argv[i][7:]+" GPa")


# Load the experimental strain data
D2=np.loadtxt('strain.txt',delimiter=',');
f=abs(D2[:,0])
psd=abs(D2[:,1])
ax1=plt.subplot(2,1,1)
plt.plot(f, psd, '--k', linewidth=1, label="Exp.")
plt.axvline(x=1/10, color='black', linewidth=0.5, linestyle='--')
plt.axvline(x=1/15, color='black', linewidth=0.5, linestyle='--')
plt.axvline(x=1/35, color='black', linewidth=0.5, linestyle='--')
ax1.text(1/10-0.006, 0.5e10, "$T=10$ s")
ax1.text(1/15-0.006, 0.5e-7, "$T=15$ s")
ax1.text(1/35-0.006, 0.5e10, "$T=35$ s")
plt.xlabel("$\omega/2\pi$ (in s$^{-1}$)")
plt.ylabel("$\\varepsilon^2 \,PSD$ (microstrain$^2$ s)")
plt.xlim(0.01,0.4)
plt.ylim(0, 40)
num_samples=5
psdRms = rolling_rms(f, psd, num_samples)

ax2 = plt.subplot(2,1,2)
plt.plot(f[num_samples::], psdRms, '--k', linewidth=1, label="Exp.")
plt.axvline(x=1/10, color='black', linewidth=0.5, linestyle='--')
plt.axvline(x=1/15, color='black', linewidth=0.5, linestyle='--')
plt.axvline(x=1/35, color='black', linewidth=0.5, linestyle='--')
ax2.text(1/10-0.006, 1e4, "$T=10$ s")
ax2.text(1/15-0.006, 1e-4, "$T=15$ s")
ax2.text(1/35-0.006, 1e4, "$T=35$ s")
plt.xlabel("$\omega/2\pi$ (in s$^{-1}$)")
plt.ylabel("$\\varepsilon_{RMS}$ (microstrain)")
plt.xlim(0.01,0.4)
plt.ylim(0, 20)

ax1.legend(loc='upper right')
ax2.legend(loc='upper right')

plt.savefig("RMSCompare.pdf", bbox_inches="tight")
plt.show()
