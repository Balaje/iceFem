%% Program to study the vibration of the Brunt Ice--shelf.

clc
clear
close all
format long

[~,~,~,~,E,nu,rhow,rhoi,g,~]=getProperties();
a=2*pi/200;
b=2*pi/400;
omega=linspace(a,b,51); %Solved in the HPC grid.
T=2*pi./omega;
Ad=1;
Ap=(g./(1i*omega))*Ad;


% Set the new-frequency space;
npts=500;
a1=a; b1=b;
omegaNew=linspace(a1,b1,npts+1);
Tnew=2*pi./omegaNew;
ApNew=(g./(1i*omegaNew))*Ad;
file1='3_BEDMAP2/';
rc=zeros(length(omega),1);
for m=1:length(omega)
    RC=load([file1,'2_RefCoeff/RefCoeff',num2str(m-1),'.dat']);
    rc(m)=RC(1)+1i*RC(2);
end

% Interpolate the reflection coefficient.
filePath=[file1,'2_RefCoeff'];
nev=60;
rd=zeros(length(omega),1);
rr=zeros(length(omega),nev);
lambdaj=zeros(nev,length(omega));
rrc=zeros(length(omega),1);
for m=1:length(omega)
    rcDiff = load([filePath,'/RefCoeff_Dif/refC',num2str(m-1),'.dat']);
    rcRad = load([filePath,'/RefCoeff_Rad/refC',num2str(m-1),'.dat']);    
    rd(m) = rcDiff(1)+1i*rcDiff(2);
    rr(m,:) = (1/20)*(rcRad(:,1)+1i*rcRad(:,2)).';
    
    lam = load([file1,'/2_ModesMatrix/lambdaj',num2str(m-1),'.dat']);
    lam = (lam(:,1)+1i*lam(:,2)).';
    
    lambdaj(:,m) = lam;
    rrc(m)=rd(m)+rr(m,:)*lambdaj(:,m);
end
rdNew=interp1(omega,rd,omegaNew);
rrNew=zeros(length(omegaNew),nev);
for m=1:nev
    rrNew(:,m)=interp1(omega,rr(:,m),omegaNew);
end

figure(1);
subplot(1,2,1);
plot(Tnew,real(rrNew(:,1:5)),'b',T,real(rr(:,1:5)),'r+');
subplot(1,2,2);
plot(Tnew,real(rdNew),'b',T,real(rd),'r+');



%% Interpolate the linear system.
filePath = [file1,'2_ModesMatrix/'];
[omegaNew,detH,condH] = interpolateFreq(a1,b1,omega,nev,filePath,npts,1);

LAMRe=load([filePath,'Interpolated_L/lambdaRe.dat']);
LAMIm=load([filePath,'Interpolated_L/lambdaIm.dat']);
LAM=LAMRe+1i*LAMIm;
rcNew=zeros(length(omegaNew),1);
for m=1:length(omegaNew)
    rcNew(m)=rdNew(m);
    for n=1:length(nev)
        rcNew(m)=rcNew(m)+(LAM(m,n)*rrNew(m,n));
    end
end

figure(2);
subplot(2,1,1);
semilogx(Tnew,real(rcNew),Tnew,imag(rcNew));
subplot(2,1,2);
semilogx(Tnew,condH);
