# Interpolation Routines for the frequency domain solution
# 1) def interpolateFreq(a,b,omega,Nev,filePath,npts,isSolve):
#        Implements the interpolation of reduced system on real
# frequency space.
# 2) def interpolateFreqComplex(a,b,omega,Nev,filePath,npts,isSolve):
#        Implements the interpolation of reduced system on complex
# frequency space.
# 3) def interpolateRefCoeff(a,b,omega,Nev,filePath,npts,TorC):
#        Implements the interpolation of radiation and diffraction
# reflection coefficients.



from scipy import interpolate
from numpy import loadtxt
import numpy as np

def interpolateFreq(a,b,omega,Nev,filePath,npts,isSolve):
    omegaNew=np.linspace(a,b,npts+1)
    H=np.zeros((Nev**2,len(omega)),dtype=complex)
    F=np.zeros((Nev,len(omega)),dtype=complex)
    for iter in range(0,len(omega)-1):
        HRe=loadtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=loadtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*Him
        if(isSolve):
            FRe=loadtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
            FIm=loadtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
            F[:,iter]=FRe+1j*FIm
    HNew=np.zeros((Nev**2,len(omegaNew)),dtype=complex)
    FNew=np.zeros((Nev,len(omegaNew)),dtype=complex)
    for m in range(0,Nev**2-1):
        ReH=H[m,:].real
        ImH=H[m,:].imag
        f=interpolate.interp1d(omega,ReH)
        ReHnew=f(omegaNew)
        f=interpolate.interp1d(omega,ImH)
        ImHnew=f(omegaNew)
        HNew[m,:]=ReHnew+1j*ImHnew
        if(m<Nev & isSolve):
            ReF=F[m,:].real
            ImF=F[m,:].imag
            f=interpolate.interp1d(omega,ReF)
            ReFnew=f(omegaNew)
            f=interpolate.interp1d(omega,ImF)
            ImFnew=f(omegaNew)
    np.savetxt(filePath+"Interpolated_H/ReH.dat",ReHnew,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_H/ImH.dat",ImHnew,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ReF.dat",ReFnew,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ImF.dat",ImFnew,delimiter="\t",newline="\n")
    if(isSolve):
        lambdaj=np.zeros((len(omegaNew),Nev))
        for p in range(0,len(omegaNew)-1):
            Hmat=np.reshape(HNew[:,p],(Nev,Nev),order='F')
            Fmat=FNew[:,p]
            lambdaj[p,:]=np.linalg.solve(Hmat,Fmat)
        np.savetxt(filePath+"Interpolated_L/lambdaRe.dat",lambdaj.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_L/lambdaIm.dat",lambdaj.imag,delimiter="\t",newline="\n")
    return omegaNew


def interpolateFreqComplex(omega,omegaNew,Nev,filePath):
    omegaflat=omega.flatten(order='F')
    omeganewflat=omegaNew.flatten(order='F')
    H=np.zeros((Nev**2,len(omegaflat)),dtype=complex)
    F=np.zeros((Nev**2,len(omeganewflat)),dtype=complex)
    for iter in range(0,len(omegaflat)-1):
        HRe=loadtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=loadtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*Him
        FRe=loadtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
        FIm=loadtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
        F[:,iter]=FRe+1j*FIm
    HNew=np.zeros((Nev**2,len(omegaNew)),dtype=complex)
    FNew=np.zeros((Nev,len(omegaNew)),dtype=complex)
    for m in range(0,Nev**2-1):
        ReH=np.reshape(H[m,:],(np.shape(omega)[0],np.shape(omega)[1]))
        f=interp2d(omegaflat.real, omegaflat.imag, ReH)
        HH=f(omeganewflat.real,omeganewflat.imag)
        HNew[m,:]=HH

        if(m<Nev):
            ReF=np.reshape(F[m,:],(np.shape(omega)[0],np.shape(omega)[1]))
            f=interp2d(omegaflat.real,omegaflat.imag,ReF)
            FF=f(omeganewflat.real,omeganewflat.imag)
            FNew[m,:]=FF
    np.savetxt(filePath+"Interpolated_H/ReH.dat",HNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_H/ImH.dat",HNew.imag,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ReF.dat",FNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ImF.dat",FNew.imag,delimiter="\t",newline="\n")
    lambdaj=np.zeros((len(omegaNew),Nev))
    for p in range(0,len(omegaNew)-1):
        Hmat=np.reshape(HNew[:,p],(Nev,Nev),order='F')
        Fmat=FNew[:,p]
        lambdaj[p,:]=np.linalg.solve(Hmat,Fmat)
    np.savetxt(filePath+"Interpolated_L/lambdaRe.dat",lambdaj.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_L/lambdaIm.dat",lambdaj.imag,delimiter="\t",newline="\n")

    return 0
