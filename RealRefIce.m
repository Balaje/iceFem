%% Run the following command in the shell to obtain an initial set of frequency domain solutions.

% for i in $(seq 15 2.5 80); do mpirun -np 2 FreeFem++-mpi -v 0 iceberg.edp -N1 20 -N2 30 -Tr $i -L 3000 -H 2000 -h 200 -nev 8 -iter $(echo $i/2.5-5 | bc) > /dev/null; echo "Done $i"; done 

% This may take a while depending on your computer. Once the command
% finishes running, run this MATLAB script!

clc
clear

SolutionDir='TEST/';

%% From the finite element problems; (Sample Shell command to run the frequency domain solutions for iceberg.edp)
npts=27;
RC=zeros(npts,4);
for m=1:npts
    RC(m,:)=load([SolutionDir,'2_RefCoeff/rc',num2str(m),'.dat']);
end
RC1=RC(:,1)+1i*RC(:,2);
TC=RC(:,3)+1i*RC(:,4);

%% Interpolate the linear system and solve
a=15;
b=80;
T=linspace(a,b,npts);
omega=2*pi./T;
NptsNew=300;
[omegaNew,detH,conH]=interpolateFreq(2*pi/b,2*pi/a,omega,8,[SolutionDir,'2_ModesMatrix/'],NptsNew-1,1);
ReLAMBDA=load([SolutionDir,'2_ModesMatrix/Interpolated_L/lambdaRe.dat']);
ImLAMBDA=load([SolutionDir,'2_ModesMatrix/Interpolated_L/lambdaIm.dat']);
LAM=ReLAMBDA+1i*ImLAMBDA;

%% Interpolate Reflection Coeffs
V1=interpRefCoeffs(omega,omegaNew,8,[SolutionDir,'2_RefCoeff'],'C');
V2=interpRefCoeffs(omega,omegaNew,8,[SolutionDir,'2_RefCoeff'],'T');
% and load
RefCDifRe=load([SolutionDir,'2_RefCoeff/Interpolated_R/refCDifRe.dat']);
RefCDifIm=load([SolutionDir,'2_RefCoeff/Interpolated_R/refCDifIm.dat']);
RefTDifRe=load([SolutionDir,'2_RefCoeff/Interpolated_R/refTDifRe.dat']);
RefTDifIm=load([SolutionDir,'2_RefCoeff/Interpolated_R/refTDifIm.dat']);

RCDif=RefCDifRe+1i*RefCDifIm;
RTDif=RefTDifRe+1i*RefTDifIm;

RefCRadRe=load([SolutionDir,'2_RefCoeff/Interpolated_R/refCRadRe.dat']);
RefCRadIm=load([SolutionDir,'2_RefCoeff/Interpolated_R/refCRadIm.dat']);
RefTRadRe=load([SolutionDir,'2_RefCoeff/Interpolated_R/refTRadRe.dat']);
RefTRadIm=load([SolutionDir,'2_RefCoeff/Interpolated_R/refTRadIm.dat']);

RCRad=RefCRadRe+1i*RefCRadIm;
RTRad=RefTRadRe+1i*RefTRadIm;

RC=zeros(NptsNew,1);
RT=zeros(NptsNew,1);

%% Reconstruct RefCoeffs
for indx=1:NptsNew
    L=LAM(indx,:);
    RC(indx)=RCDif(indx)+dot(conj(RCRad(indx,:)),L);
    RT(indx)=RTDif(indx)+dot(conj(RTRad(indx,:)),L);
end

fig=figure(1);
set(fig,'Position',[273   123   967   582]);
subplot(3,1,1);
plot(omega,abs(RC1),'+','linewidth',2);
hold on
plot(omegaNew,abs(RC),'linewidth',2)
xlabel('$\omega$');
ylabel('$R(\omega)$');
legend('Coarse $\omega$ Space','Fine $\omega$ Space')
xlim([omegaNew(1),omegaNew(end)])

subplot(3,1,2);
plot(omegaNew,sqrt(abs(RT).^2+abs(RC).^2),'k','linewidth',2);
hold on
plot(omegaNew,abs(RT),'r','linewidth',2);
plot(omegaNew,abs(RC),'b','linewidth',2);
xlabel('$\omega$');
legend('1','$T(\omega)$','$R(\omega)$');
xlim([omegaNew(1),omegaNew(end)])

subplot(3,1,3);
for indx=1:4
    L=LAM(:,indx);
    plot(omegaNew,abs(L),'DisplayName',['$|\lambda_{',num2str(indx),'}|$'],'linewidth',2);    
    hold on
end
legend show
title('Coefficients vs $\omega$');
xlabel('$\omega$');
ylabel('$|\lambda_j|$');
xlim([omegaNew(1),omegaNew(end)])
%ylim([0,0.8])

%% Function
function V=interpRefCoeffs(omega,omegaNew,Nev,filePath,TorC)
rd = zeros(length(omega),1); %Diffraction RefCoeffs
rr = zeros(length(omega),Nev); %Radiation RefCoeffs
for m=1:length(omega)
    rcDiff = load([filePath,'/RefCoeff_Dif/ref',TorC,num2str(m),'.dat']);
    rcRad = load([filePath,'/RefCoeff_Rad/ref',TorC,num2str(m),'.dat']);
    
    rd(m) = rcDiff(1)+1i*rcDiff(2);
    rr(m,:) = (rcRad(:,1)+1i*rcRad(:,2)).';
end
rdNewRe = interp1(omega,real(rd),omegaNew,'pchip');
rdNewIm = interp1(omega,imag(rd),omegaNew,'pchip');
% Loop over each mode and interpolate the reflection coefficients
rcNew = zeros(length(omegaNew),Nev);
for m=1:Nev
    rccRe = interp1(omega,real(rr(:,m)),omegaNew,'pchip');
    rccIm = interp1(omega,imag(rr(:,m)),omegaNew,'pchip');
    rcc=rccRe+1i*rccIm;
    rcNew(:,m) = rcc;
end

%% Write the interpolated reflection coefficients to a file.
dlmwrite([filePath,'/Interpolated_R/ref',TorC,'DifRe.dat'],rdNewRe);
dlmwrite([filePath,'/Interpolated_R/ref',TorC,'DifIm.dat'],rdNewIm);
dlmwrite([filePath,'/Interpolated_R/ref',TorC,'RadRe.dat'],real(rcNew));
dlmwrite([filePath,'/Interpolated_R/ref',TorC,'RadIm.dat'],imag(rcNew));

V=0;
end
