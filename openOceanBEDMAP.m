%% Routine to plot the open--ocean solution

clc
clear
close all

SolDir='1_BEDMAP2/';
nfreq=200;
a=2*pi*0.0001;
b=2*pi*0.01;
omega=linspace(a,b,200);
domega=omega(2:end)-omega(1:end-1);
domega=[domega,0];

LAMRe=load([SolDir,'2_ModesMatrix/Interpolated_L/lambdaRe.dat']);
LAMIm=load([SolDir,'2_ModesMatrix/Interpolated_L/lambdaRe.dat']);
X=load([SolDir,'2_Modes/PHI/XVAL.dat']);
PHIFreq=zeros(length(X),nfreq);

for m=1:nfreq    
    PHI0=load([SolDir,'2_Modes/PHI/phi0FS',num2str(m),'.dat']);
    PHIJ=load([SolDir,'2_Modes/PHI/phijFS',num2str(m),'.dat']);
    PHI0Re=PHI0(:,1);   PHI0Im=PHI0(:,2);
    PHIJRe=PHIJ(:,1);   PHIJIm=PHIJ(:,2);
    
    PHIJ1=reshape(PHIJRe,[4036,64])+1i*reshape(PHIJIm,[4036,64]);
    PHI01=PHI0Re+1i*PHI0Im;    
    
    LAM=LAMRe(m,:)+1i*LAMIm(m,:);
    PHIFreq(:,m)=(PHI01+PHIJ1*LAM.')*1i*omega(m)/9.8;
    
    fprintf(['Imported ',num2str(m),' ... \n']);
end
X=flipud(X);


%% Find the time-domain solution
close all
FAmp=load('FAmp.dat');
t=-5000:10:5000;
phi=zeros(length(X),1);
for m=1:length(t)
    Arr=exp(-1i*omega*t(m) - 1i*FAmp(3,:)).*domega.*FAmp(2,:);
    phi=PHIFreq*Arr.';    
    
    M=[X,real(phi)];
    filename=[SolDir,'2_Potential/phiSurf',num2str(m),'.csv'];
    csvwrite(filename,M);
    plot(X,real(phi),'b-');
    ylim([-2e-6,2e-6]);
    xlim([X(1),X(end)]);
    pause(0.01);    
end