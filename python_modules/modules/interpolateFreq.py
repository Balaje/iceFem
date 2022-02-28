################################################################################
# Interpolation Routines for the frequency domain solution:
#
# (1) def interpolateCoeffsFreq(a,b,omega,Nev,filePath,npts,isSolve):
# Implements the interpolation of reduced system on real frequency space
# Inputs
#   a: omega(1)
#   b: omega(end)
#   omega: An array containing the original frequency space.
#   Nev: Number of in-vacuo modes used in the construction of the solutions
#   filePath: Full path to the Working Directory. Eg. `1_SIMPLE5/`
#   npts: Number of interpolation points
#   isSolve: Option (0/1) to solve the reduced system to obtain the modal coefficients
#
# 2) def interpolateFreqComplex(a,b,omega,Nev,filePath,npts,isSolve):
#        Implements the interpolation of reduced system on complex
# frequency space.
#
# 3) def interpolateRefCoeff(a,b,omega,Nev,filePath,npts,TorC):
#        Implements the interpolation of radiation and diffraction
# reflection coefficients.
#
# 4) interpolateRefCoeffComplex(a,b,c,d,npts, Nev, filePath, TorC, nptsNew):
#        Implements the interpolation of radiation and diffraction
# reflection coefficients on the complex plane
#
# 5) def buildLam(SolutionDir):
#        Build the coefficient vector using the interpolated reduced system.
#
# 6) def buildRMat(LAM,SolutionDir,TorC,p):
#        Build the reflection/transmission coefficient using the solution
# LAM = buildLam(SolutionDir):
#
# 7) def genComplexPlot(RC,omeganew):
#        Generate complex plot using the reflection coefficient and (complex) omeganew.
###################################################################################


from scipy import interpolate
import numpy as np
import cplot
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

def interpolateCoeffsFreq(a,b,omega,Nev,filePath,npts,isSolve):
    omegaNew=np.linspace(a,b,npts+1)
    H=np.zeros((Nev**2,len(omega)),dtype=complex)
    F=np.zeros((Nev,len(omega)),dtype=complex)
    for iter in range(0,len(omega)):
        HRe=np.genfromtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=np.genfromtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*HIm
        if(isSolve):
            FRe=np.genfromtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
            FIm=np.genfromtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
            F[:,iter]=FRe+1j*FIm
    HNew=np.zeros((Nev**2,len(omegaNew)),dtype=complex)
    FNew=np.zeros((Nev,len(omegaNew)),dtype=complex)

    for m in range(0,Nev**2):
        f=interpolate.interp1d(omega,H[m,:],kind='cubic')
        HNew[m,:]=f(omegaNew);
        if(m<Nev and isSolve==1):
            f=interpolate.interp1d(omega,F[m],kind='cubic')
            FNew[m,:]=f(omegaNew)

#    for m in range(0,4):
#        plt.semilogy(omega,abs(F[m,:]),marker='+',linestyle=' ',label="Original F["+str(m)+"]")
#    for m in range(0,4):
#        plt.semilogy(omegaNew,abs(FNew[m,:]),label="Interpolated F["+str(m)+"]")
#    plt.legend()
#    heading="Entries vs Frequency ("+filePath+")";
#    plt.title(heading)
#    plt.show()

    np.savetxt(filePath+"Interpolated_H/ReH.dat",HNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_H/ImH.dat",HNew.imag,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ReF.dat",FNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ImF.dat",FNew.imag,delimiter="\t",newline="\n")
    condH=np.zeros(len(omegaNew));
    if(isSolve):
        lambdaj=np.zeros((len(omegaNew),Nev),dtype=complex)
        for p in range(0,len(omegaNew)):
            Hmat=np.reshape(HNew[:,p],(Nev,Nev),order='F')
            condH[p]=abs(np.linalg.cond(Hmat))
            Fmat=FNew[:,p]
            lambdaj[p,:]=np.linalg.solve(Hmat,Fmat)
        np.savetxt(filePath+"Interpolated_L/lambdaRe.dat",lambdaj.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_L/lambdaIm.dat",lambdaj.imag,delimiter="\t",newline="\n")

    return omegaNew


def interpolateCoeffsFreqComplex(a,b,c,d,npts,Nev,filePath,nptsNew):
    x=np.linspace(a,b,npts)
    y=np.linspace(c,d,npts)
    X,Y=np.meshgrid(x,y)
    omega=X+1j*Y

    xq=np.linspace(a,b,nptsNew)
    yq=np.linspace(c,d,nptsNew)
    Xq,Yq=np.meshgrid(xq,yq)
    omegaNew=Xq+1j*Yq

    H=np.zeros((Nev**2,npts**2),dtype=complex)
    F=np.zeros((Nev,npts**2),dtype=complex)
    for iter in range(0,npts**2):
        HRe=np.genfromtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=np.genfromtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*HIm
        FRe=np.genfromtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
        FIm=np.genfromtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
        F[:,iter]=FRe+1j*FIm
    HNew=np.zeros((Nev**2,nptsNew**2),dtype=complex)
    FNew=np.zeros((Nev,nptsNew**2),dtype=complex)
    for m in range(0,Nev**2):
        f1=interpolate.interp2d(x, y, np.reshape(H[m,:].real,np.shape(omega)))
        f2=interpolate.interp2d(x, y, np.reshape(H[m,:].imag,np.shape(omega)))
        HH1=f1(xq,yq)
        HH2=f2(xq,yq)
        HNew[m,:]=HH1.flatten(order='F')+1j*HH2.flatten(order='F')
        if(m<Nev):
            f1=interpolate.interp2d(x, y, np.reshape(F[m,:].real,np.shape(omega)))
            f2=interpolate.interp2d(x, y, np.reshape(F[m,:].imag,np.shape(omega)))
            FF1=f1(xq,yq)
            FF2=f2(xq,yq)
            FNew[m,:]=FF1.flatten(order='F')+1j*FF2.flatten(order='F')

    np.savetxt(filePath+"Interpolated_H/ReH.dat",HNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_H/ImH.dat",HNew.imag,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ReF.dat",FNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_F/ImF.dat",FNew.imag,delimiter="\t",newline="\n")
    lambdaj=np.zeros((nptsNew**2,Nev),dtype=complex)
    for p in range(0,nptsNew**2):
        Hmat=np.reshape(HNew[:,p],(Nev,Nev),order='F')
        Fmat=FNew[:,p]
        lambdaj[p,:]=np.linalg.solve(Hmat,Fmat)
    np.savetxt(filePath+"Interpolated_L/lambdaRe.dat",lambdaj.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_L/lambdaIm.dat",lambdaj.imag,delimiter="\t",newline="\n")
    return 0

def interpolateRefCoeff(omega,omegaNew,Nev,filePath,NModes,TorC):
    rd=np.zeros((len(omega),NModes+1),dtype=complex)
    rr=np.zeros((len(omega),NModes+1,Nev),dtype=complex)
    for p in range(0,NModes+1):
        for m in range(0,len(omega)):
            rcDiff=np.genfromtxt(filePath+"RefCoeff_Dif/ref"+TorC+str(m+1)+".dat",delimiter="\t",dtype=None,encoding=None)
            rcRad=np.genfromtxt(filePath+"RefCoeff_Rad/ref"+TorC+str(m+1)+".dat",delimiter="\t",dtype=None,encoding=None)
            rd[m,:]=rcDiff[:,0]+1j*rcDiff[:,1]
            rr[m,p,:]=(rcRad[:,0:2:2*NModes+1]+1j*rcRad[:,1:2:2*NModes+1]).reshape(-1)

        f=interpolate.interp1d(omega,rd[:,p],kind='cubic')
        rdNew=f(omegaNew)
        rcNew=np.zeros((len(omegaNew),Nev),dtype=complex)
        for n in range(0,Nev):
            f=interpolate.interp1d(omega,rr[:,p,n],kind='cubic')
            rcNew[:,n]=f(omegaNew)

        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifRe_Mode"+str(p)+".dat",rdNew.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifIm_Mode"+str(p)+".dat",rdNew.imag,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadRe_Mode"+str(p)+".dat",rcNew.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadIm_Mode"+str(p)+".dat",rcNew.imag,delimiter="\t",newline="\n")

def interpolateRefCoeffComplex(a,b,c,d,npts, Nev, filePath, TorC, nptsNew, NModes):
    x=np.linspace(a,b,npts)
    y=np.linspace(c,d,npts)
    X,Y=np.meshgrid(x,y)
    omega=X+1j*Y

    xq=np.linspace(a,b,nptsNew)
    yq=np.linspace(c,d,nptsNew)
    Xq,Yq=np.meshgrid(xq,yq)
    omegaNew=Xq+1j*Yq

    rd=np.zeros((NModes+1,npts**2),dtype=complex)
    rr=np.zeros((NModes+1,npts**2,Nev),dtype=complex)
    for m in range(0,npts**2):
        rcDiff=np.genfromtxt(filePath+"RefCoeff_Dif/ref"+TorC+str(m+1)+".dat")
        rcRad=np.genfromtxt(filePath+"RefCoeff_Rad/ref"+TorC+str(m+1)+".dat")
        rd[:,m]=rcDiff[:,0]+1j*rcDiff[:,1]
        rr[:,m,:]=(rcRad[:,::2]+1j*rcRad[:,1::2]).transpose()

    for p in np.arange(0,NModes):
        reshapeRd=np.reshape(rd[p,:],np.shape(omega))
        fR=interpolate.interp2d(x, y, reshapeRd.real)
        fI=interpolate.interp2d(x, y, reshapeRd.imag)
        rdNew=fR(xq, yq)+1j*fI(xq, yq)
        rdNew=rdNew.flatten(order='F')

        rcNew=np.zeros((nptsNew**2,Nev),dtype=complex)
        for m in range(0,Nev):
            reshapeRc=np.reshape(rr[p,:,m],np.shape(omega))
            fR=interpolate.interp2d(x,y,reshapeRc.real)
            fI=interpolate.interp2d(x,y,reshapeRc.imag)
            rcc=fR(xq, yq)+1j*fI(xq, yq)
            rcNew[:,m]=rcc.flatten(order='F')

        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifRe_Mode"+str(p)+".dat",rdNew.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifIm_Mode"+str(p)+".dat",rdNew.imag,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadRe_Mode"+str(p)+".dat",rcNew.real,delimiter="\t",newline="\n")
        np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadIm_Mode"+str(p)+".dat",rcNew.imag,delimiter="\t",newline="\n")


def buildLam(SolutionDir):
    LAMRe=np.genfromtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
    LAMIm=np.genfromtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
    LAM=LAMRe+1j*LAMIm
    return LAM

def buildRMat(LAM,SolutionDir,TorC,p):
    RefTDifRe=np.genfromtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"DifRe_Mode"+str(p)+".dat");
    RefTDifIm=np.genfromtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"DifIm_Mode"+str(p)+".dat");
    RefTRadRe=np.genfromtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"RadRe_Mode"+str(p)+".dat");
    RefTRadIm=np.genfromtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"RadIm_Mode"+str(p)+".dat");
    RTDif=RefTDifRe+1j*RefTDifIm
    RTRad=RefTRadRe+1j*RefTRadIm
    RT=np.zeros(len(RTDif),dtype=complex)
    for indx in range(0,len(RT)):
        L=LAM[indx,:]
        RT[indx]=RTDif[indx]+np.dot(RTRad[indx,:],L)
    return RT

def genComplexPlot(RC,omeganew):
    RCNew=np.reshape(RC,np.shape(omeganew))
    RCNew=RCNew.T
    VAL=cplot.get_srgb1(RCNew,colorspace="cam16")
    VAL=np.flipud(VAL)
    a=omeganew[0,0].real; b=omeganew[-1,-1].real
    c=omeganew[0,0].imag; d=omeganew[-1,-1].imag
    plt.imshow(VAL,extent=[a,b,c,d])
    plt.xlabel("$Re(\omega)$",usetex=True)
    plt.ylabel("$Im(\omega)$",usetex=True)
