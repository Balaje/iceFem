##############################################
# Function to non-dimensionalize the problem
##############################################
import numpy as np

def nonDimensionalize(L, H, h, E, omega):
    nu = 0.33
    rhoi = 922.5; rhow = 1025
    ag = 9.8; Ad = 1
    elasCons = [E/(2*(1+nu)), E*nu/((1+nu)*(1-2*nu))]
    EI = E*h**3/(12*(1-nu)**2)
    Lc = (EI/(rhow*ag))**(0.25)
    tc = np.sqrt(rhow*Lc**6/(EI*H))
    d = (rhoi/rhow)*h
    LL = L/Lc;    HH = H/Lc;    hh = h/Lc
    ndOmega = omega*tc
    alpha = HH*ndOmega**2
    Ap = ag/(1j*omega)*Ad
    return [LL, HH, hh, alpha]
