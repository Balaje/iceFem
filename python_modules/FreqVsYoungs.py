######################################################################
# Python script for studying Real part of omega vs Young's modulus E
#
# Run:
#                  python3 FreqVsYoungs.py
###############################################

import numpy as np
from modules.interpolateFreq import *
import cplot
import matplotlib.pyplot as plt
from os import system
import sys
import math
# Import eigenvalue related function
from scipy.linalg import eig

pi = np.pi


def solve_one_ice_problem(E, omega):
    solveFreeFemProblem('iceshelf2d.edp', 3630, 500, 280, 12, omega, E, 10, '0', 'ICEBERG_COMPLEX/', 0.01, 1)
    H = read_H_matrix('0')
    return H

def computeResonanceFrequency(E, omega0, tol, ind):
    count = 1
    dw = 1e-4
    omega0s = [omega0]
    while (abs(dw)/abs(omega0) > tol) and (np.isinf(dw) == 0):
        H0 = solve_one_ice_problem(E, omega0)
        H1 = solve_one_ice_problem(E, omega0 + dw)

        dH = (H1-H0)/dw
        domega, dwvec = eig(H0, -dH)
        dw = domega[np.argmin(abs(domega))]
        detH = np.abs(np.product(domega))
        detH0 = np.abs(np.linalg.det(H0))
        condH = np.linalg.cond(H0)
        print("det(H0) = "+str("{:.4e}".format(detH0))+
              "\tProd(domega) (determinant) = "+
              str("{:.4e}".format(detH))+"\tcond(H0) = "+
              str("{:.4e}".format(condH))+
              "\tMin(domega) = "+
              str("{:.4e}".format(dw)))
        if (np.isinf(abs(dw)/abs(omega0)) == 0):
            omega0 = omega0 + dw
            omega0s.append(omega0)
    omega0s = np.array(omega0s)
    detH = np.product(domega)
    print("det(H) = "+str(detH))
    return omega0s

######################################################
# Real part of resonance freq. vs Young's Modulus
######################################################
Es = np.linspace(2.5, 5, 20)
#omega0 = 2*pi*0.012 + 0*1j
#omega0 = 2*pi*0.022 + 0*1j
#omega0 = 2*pi*0.051 + 0*1j
omega0 = 2*pi*0.1 + 0*1j
omegar = []
for m in range(0, np.size(Es)):
    omega0s = computeResonanceFrequency(Es[m]*1e9, omega0, 1e-8, 100)
    omegar.append(omega0s[-1])
    omega0 = omega0s[-1].real
    H0 = solve_one_ice_problem(Es[m]*1e9, omegar[m])
    print(str(omegar[m])+" for "+str(Es[m])+"*10^9 GPa")

# Save the output
np.savetxt("ICEBERG_COMPLEX/omegar_vs_youngs.txt", np.array(omegar),delimiter="\t", newline="\n")
omegar0 = np.genfromtxt("ICEBERG_COMPLEX/omegar_vs_youngs0.txt", dtype=complex)/(2*pi)
omegar1 = np.genfromtxt("ICEBERG_COMPLEX/omegar_vs_youngs1.txt", dtype=complex)/(2*pi)
omegar2 = np.genfromtxt("ICEBERG_COMPLEX/omegar_vs_youngs2.txt", dtype=complex)/(2*pi)
#omegar3 = np.genfromtxt("ICEBERG_COMPLEX/omegar_vs_youngs3.txt", dtype=complex)/(2*pi)
omegar3 = omegar

omega_exp3 = np.ones((30,))*(2*pi*(1/10))/(2*pi)
omega_exp2 = np.ones((30,))*(2*pi*(1/15))/(2*pi)
omega_exp1 = np.ones((30,))*(2*pi*(1/50))/(2*pi)

plt.figure(figsize = [13,6])
ax = plt.subplot(1,2,1)
# Plot the experimental and theoretical frequencies vs E
# omega1
plt.plot(Es, omegar1.real, 'r-', linewidth=2)
plt.plot(Es, omega_exp1, 'k--', linewidth=1.5)
plt.text(2.1, omega_exp1[1].real+0.003, "$T_{exp} = 50$ s", usetex=True)
# omega2
plt.plot(Es, omegar2.real, 'b-', linewidth=2)
plt.plot(Es, omega_exp2, 'k-.', linewidth=1.5)
plt.text(2.1, omega_exp2[1].real+0.003, "$T_{exp} = 15$ s", usetex=True)
# omega3
plt.plot(Es, omegar3.real, 'g-', linewidth=2)
plt.plot(Es, omega_exp3, 'k-o', linewidth=1.5)
plt.text(2.1, omega_exp3[1].real+0.003, "$T_{exp} = 10$ s", usetex=True)
# Range of omega
plt.fill_between(Es, 0, 1, where= (Es >= 1) & (Es <= 3.5),
                 color='green', alpha = 0.3, transform=ax.get_xaxis_transform())
plt.ylim(0,0.15)
plt.xlim(Es[0], Es[-1])
plt.xlabel("Young's modulus $E$ (in GPa)", usetex=True)
plt.ylabel("$\omega_r/2\pi \quad$ (in s$^{-1}$)", usetex=True)

## Plot the error in the frequency vs E
err = np.sqrt(np.abs(omegar1.real - omega_exp1)**2
              + np.abs(omegar2.real - omega_exp2)**2
              + np.abs(omegar3.real - omega_exp3)**2)
ax = plt.subplot(1,2,2)
plt.plot(Es, err, linewidth=1.5)
ind = np.where(err == np.min(err))
plt.axvline(x=Es[ind], color='black', linewidth=1, linestyle='--')
#plt.grid()
plt.xlabel("Young's modulus $E$ (in GPa)", usetex=True)
plt.ylabel("$\sqrt{\Sigma |\omega_r - \omega_{exp}|^2}\quad$ (in s$^{-1}$)", usetex=True)
Eopt = Es[ind]
plt.text(Es[ind]+0.04, 0.035, "$E_{opt} = "+str(np.round(Eopt[-1], 4))+"$ GPa")
plt.xlim(Es[0], Es[-1])
plt.savefig("omegarVsYoungs.pdf", bbox_inches='tight')
plt.show()
