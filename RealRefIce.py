import numpy as np
from modules.interpolateFreq import *
import matplotlib.pyplot as plt
import sys
import time

filePath=sys.argv[1];

T=np.linspace(10,100,51)
pi=np.pi
omega=2*pi/T
RC=np.zeros((21,4),dtype=complex)

npts=199
nev=16
omeganew=interpolateCoeffsFreq(2*pi/100,2*pi/10,omega,nev,filePath+"2_ModesMatrix/",npts,1)
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
V=interpolateRefCoeff(omega,omeganew,nev,filePath+"2_RefCoeff/","C")
RC=buildRMat(LAM,filePath,"C")
plt.figure(figsize=[10,4])
plt.plot(omeganew/(2*pi),np.transpose(abs(RC)),linewidth=2,label="R($\omega$)")

## Interpolating Transmission coefficients
V=interpolateRefCoeff(omega,omeganew,nev,filePath+"2_RefCoeff/","T")
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

#print("Sleep")
#time.sleep(10)

#Pearson Moscowitz Spectrum.
g=9.8
omega0=(g/15)
alpha=0.0081
beta=0.74
omeganew1=2*pi*np.linspace(0.01,0.45,200)
S=alpha*g**2/omeganew1**5*np.exp(-beta*(omega0/omeganew1)**4)+0.05;
plt.figure(figsize=[5,4])
plt.plot(omeganew1/(2*pi),S,label='Pearson Moskowitz Spectrum')
plt.xlabel('$\omega$ (in Hz)')
plt.ylabel('PSD')
S1=np.loadtxt('/Users/balaje/data1.txt',delimiter=',')
plt.plot(S1[:,0],S1[:,1],label='Kristensen et al.')
plt.legend()
plt.savefig(sys.argv[2]+"_5.pdf",bbox_inches= 'tight')
plt.show()



omeganew1=omeganew;
S=alpha*g**2/omeganew1**5*np.exp(-beta*(omega0/omeganew1)**4)+0.05;
S=omeganew**0


plt.figure(figsize=[5,4])
maxDisp=np.loadtxt("./dispVsFreq.dat")
plt.plot(omeganew/(2*pi),S*maxDisp[:,0],label="$u_x$")
plt.plot(omeganew/(2*pi),S*maxDisp[:,1],label="$u_y$")
plt.legend()
plt.xlabel('$\omega$')
#plt.xlim([1/100,0.45])
plt.title('Displacements')
plt.savefig(sys.argv[2]+"_3.pdf",bbox_inches='tight')


plt.figure(figsize=[5,4])
maxStrain=np.loadtxt("./strainVsFreq.dat")
plt.plot(omeganew/(2*pi),S*maxStrain[:,0],label="$e_{xx}$")
plt.plot(omeganew/(2*pi),S*maxStrain[:,1],label="$e_{xy}$")
plt.plot(omeganew/(2*pi),S*maxStrain[:,2],label="$e_{yy}$")
plt.legend()
plt.xlabel('$\omega$')
plt.title('Strain')
#plt.xlim([1/100,0.45])
plt.savefig(sys.argv[2]+"_4.pdf",bbox_inches='tight')
