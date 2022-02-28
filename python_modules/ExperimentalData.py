########################################################################################################
# Python script to compute the RMS strain from experiments
# Run:
#        python3 ExperimentalData.py
#
# Data obtained from:
#      Kristensen, M., Squire, V.A., Moore, S.C., 1982. Tabular icebergs in ocean waves 297, Nature 3.
#
# Compute the RMS strain by computing the area under the curve.
########################################################################################################

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

## Frequency bin
#a=0.047
#b=0.055
a=0.094
b=0.1

D2=np.loadtxt('inc_wave.txt',delimiter=',')
fexp=abs(D2[:,0])
psd=abs(D2[:,1])
plt.figure(figsize=[5,4])
plt.plot(fexp,psd)
# Find the theoretical spectrum
al=8.1e-3
g=9.8
omegap=2*np.pi*0.093
u195=0.877*9.8/omegap
omega0=g/u195
def Spec(omega):
    return (al*g**2/omega**5)*np.exp(-0.74*(omega0/omega)**4)
SpecPierson=Spec(2*np.pi*fexp)
plt.plot(fexp,SpecPierson)
plt.xlabel('$f$ (in Hz)')
plt.ylabel('PSD (in m$^2$ s)')
plt.show()

# Strain spectrum
D2=np.loadtxt('strain.txt',delimiter=',');
f=abs(D2[:,0])
psd=abs(D2[:,1])
fig=plt.figure(figsize=[5,4])
plt.plot(f,psd)
a=0.085
b=0.1317
fq=np.linspace(a,b,2000)
psdfun=interpolate.PchipInterpolator(f,psd)
psdq=psdfun(fq)
plt.plot(fq,psdq)
amp3=np.sqrt(np.trapz(psdq,fq))*1e-6
plt.fill_between(fq,0,psdq,color=[1,0,0])
plt.xlabel('$f$ (in Hz)')
plt.ylabel('PSD (in microstrain$^2$ s)')
print("Value of strain (experimental) = ",amp3)
plt.show()
