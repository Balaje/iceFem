# Interpolation Routines for the frequency domain solution
#(1) def interpolateCoeffsFreq(a,b,omega,Nev,filePath,npts,isSolve):
  # Implements the interpolation of reduced system on real frequency space
  # Inputs
  #   a: omega(1)
  #   b: omega(end)
  #   omega: An array containing the original frequency space.
  #   Nev: Number of in-vacuo modes used in the construction of the solutions
  #   filePath: Full path to the Working Directory. Eg. `1_SIMPLE5/`
  #   npts: Number of interpolation points
  #   isSolve: Option (0/1) to solve the reduced system to obtain the modal coefficients

# 2) def interpolateFreqComplex(a,b,omega,Nev,filePath,npts,isSolve):
#        Implements the interpolation of reduced system on complex
# frequency space.
# 3) def interpolateRefCoeff(a,b,omega,Nev,filePath,npts,TorC):
#        Implements the interpolation of radiation and diffraction
# reflection coefficients.



from scipy import interpolate
from numpy import loadtxt
import numpy as np
import cplot
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

def interpolateCoeffsFreq(a,b,omega,Nev,filePath,npts,isSolve):
    omegaNew=np.linspace(a,b,npts+1)
    H=np.zeros((Nev**2,len(omega)),dtype=complex)
    F=np.zeros((Nev,len(omega)),dtype=complex)
    for iter in range(0,len(omega)):
        HRe=loadtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=loadtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*HIm
        if(isSolve):
            FRe=loadtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
            FIm=loadtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
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
            condH[p]=np.linalg.cond(Hmat)
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
        HRe=loadtxt(filePath+"ReH"+str(iter+1)+".dat",delimiter="\t")
        HIm=loadtxt(filePath+"ImH"+str(iter+1)+".dat",delimiter="\t")
        H[:,iter]=HRe+1j*HIm
        FRe=loadtxt(filePath+"ReF"+str(iter+1)+".dat",delimiter="\t")
        FIm=loadtxt(filePath+"ImF"+str(iter+1)+".dat",delimiter="\t")
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

def interpolateRefCoeff(omega,omegaNew,Nev,filePath,TorC):
    rd=np.zeros(len(omega),dtype=complex)
    rr=np.zeros((len(omega),Nev),dtype=complex)
    for m in range(0,len(omega)):
        rcDiff=np.loadtxt(filePath+"RefCoeff_Dif/ref"+TorC+str(m+1)+".dat",delimiter="\t")
        rcRad=np.loadtxt(filePath+"RefCoeff_Rad/ref"+TorC+str(m+1)+".dat",delimiter="\t")
        rd[m]=rcDiff[0]+1j*rcDiff[1]
        rr[m,:]=rcRad[:,0]+1j*rcRad[:,1]
    f=interpolate.interp1d(omega,rd,kind='cubic')
    rdNew=f(omegaNew)
    rcNew=np.zeros((len(omegaNew),Nev),dtype=complex)
    for m in range(0,Nev):
        f=interpolate.interp1d(omega,rr[:,m],kind='cubic')
        rcNew[:,m]=f(omegaNew)

    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifRe.dat",rdNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifIm.dat",rdNew.imag,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadRe.dat",rcNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadIm.dat",rcNew.imag,delimiter="\t",newline="\n")

def interpolateRefCoeffComplex(a,b,c,d,npts, Nev, filePath, TorC, nptsNew):
    x=np.linspace(a,b,npts)
    y=np.linspace(c,d,npts)
    X,Y=np.meshgrid(x,y)
    omega=X+1j*Y

    xq=np.linspace(a,b,nptsNew)
    yq=np.linspace(c,d,nptsNew)
    Xq,Yq=np.meshgrid(xq,yq)
    omegaNew=Xq+1j*Yq

    rd=np.zeros((npts**2,1),dtype=complex)
    rr=np.zeros((npts**2,Nev),dtype=complex)
    for m in range(0,npts**2):
        rcDiff=np.loadtxt(filePath+"RefCoeff_Dif/ref"+TorC+str(m+1)+".dat",delimiter="\t")
        rcRad=np.loadtxt(filePath+"RefCoeff_Rad/ref"+TorC+str(m+1)+".dat",delimiter="\t")
        rd[m]=rcDiff[0]+1j*rcDiff[1]
        rr[m,:]=rcRad[:,0]+1j*rcRad[:,1]

    reshapeRd=np.reshape(rd,np.shape(omega))
    fR=interpolate.interp2d(x, y, reshapeRd.real)
    fI=interpolate.interp2d(x, y, reshapeRd.imag)
    rdNew=fR(xq, yq)+1j*fI(xq, yq)
    rdNew=rdNew.flatten(order='F')

    rcNew=np.zeros((nptsNew**2,Nev),dtype=complex)
    for m in range(0,Nev):
        reshapeRc=np.reshape(rr[:,m],np.shape(omega))
        fR=interpolate.interp2d(x,y,reshapeRc.real)
        fI=interpolate.interp2d(x,y,reshapeRc.imag)
        rcc=fR(xq, yq)+1j*fI(xq, yq)
        rcNew[:,m]=rcc.flatten(order='F')
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifRe.dat",rdNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"DifIm.dat",rdNew.imag,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadRe.dat",rcNew.real,delimiter="\t",newline="\n")
    np.savetxt(filePath+"Interpolated_R/ref"+TorC+"RadIm.dat",rcNew.imag,delimiter="\t",newline="\n")


def buildLam(SolutionDir):
    LAMRe=np.loadtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaRe.dat")
    LAMIm=np.loadtxt(SolutionDir+"2_ModesMatrix/Interpolated_L/lambdaIm.dat")
    LAM=LAMRe+1j*LAMIm
    return LAM

def buildRMat(LAM,SolutionDir,TorC):
    RefTDifRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"DifRe.dat");
    RefTDifIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"DifIm.dat");
    RefTRadRe=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"RadRe.dat");
    RefTRadIm=np.loadtxt(SolutionDir+"2_RefCoeff/Interpolated_R/ref"+TorC+"RadIm.dat");
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
    VAL=cplot.get_srgb1(RCNew,colorspace="cam16",alpha=0.5)
    VAL=np.flipud(VAL)
    a=omeganew[0,0].real; b=omeganew[-1,-1].real
    c=omeganew[0,0].imag; d=omeganew[-1,-1].imag
    plt.imshow(VAL,extent=[a,b,c,d])
    plt.xlabel("$Re(\omega)$",usetex=True)
    plt.ylabel("$Im(\omega)$",usetex=True)
