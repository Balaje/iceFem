clc
clear

%% From the finite element problems;

npts=10;
[X,Y]=meshgrid(linspace(2*pi/80,2*pi/15,npts),linspace(-0.06,0.06,npts));
omega=X+1i*Y;
T=2*pi./omega;
RC=zeros(npts^2,4);

% count=1;
% for m=1:npts
%    for n=1:npts
%        cmd=['/usr/local/bin/mpirun -np 2 /usr/local/ff++/mpich3/bin/FreeFem++-mpi -v 0 iceberg.edp -N1 20 -N2 30 -Tr ',num2str(real(T(m,n))),' -Ti ',num2str(imag(T(m,n))), ...
%            ' -L 3000 -H 2000 -h 200 -nev 8 -iter ',num2str(count)];       
%        [a,b]=system(cmd);
%        count=count+1;
%    end
% end

for m=1:npts^2
    RC(m,:)=load(['1_ICEBERG/2_RefCoeff/rc',num2str(m),'.dat']);
end
RCOld=RC(:,1)+1i*RC(:,2);
RCOld=reshape(RCOld,[npts,npts]);
TCOld=RC(:,3)+1i*RC(:,4);
TCOld=reshape(TCOld,[npts,npts]);

NptsNew=300;
[Xq,Yq]=meshgrid(linspace(2*pi/80,2*pi/15,NptsNew),linspace(-0.06,0.06,NptsNew));
omegaNew=Xq+1i*Yq;
Tnew=2*pi./omegaNew;

V1=interpolateFreqComplex(omega,omegaNew,8,'1_ICEBERG/2_ModesMatrix');
V2=interpolateRefCoeff(omega,omegaNew,8,'1_ICEBERG/2_RefCoeff/','C');
V3=interpolateRefCoeff(omega,omegaNew,8,'1_ICEBERG/2_RefCoeff/','T');

%% Load the RefCoeffs
RefCDifRe=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refCDifRe.dat');
RefCDifIm=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refCDifIm.dat');
RefTDifRe=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refTDifRe.dat');
RefTDifIm=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refTDifIm.dat');
RCDif=RefCDifRe+1i*RefCDifIm;
RTDif=RefTDifRe+1i*RefTDifIm;

RefCRadRe=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refCRadRe.dat');
RefCRadIm=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refCRadIm.dat');
RefTRadRe=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refTRadRe.dat');
RefTRadIm=load('1_ICEBERG/2_RefCoeff/Interpolated_R/refTRadIm.dat');
RCRad=RefCRadRe+1i*RefCRadIm;
RTRad=RefTRadRe+1i*RefTRadIm;

RC=zeros(NptsNew,NptsNew);
RT=zeros(NptsNew,NptsNew);

%% Reconstruct RefCoeffs
ReLAMBDA=load('1_ICEBERG/2_ModesMatrix/Interpolated_L/lambdaRe.dat');
ImLAMBDA=load('1_ICEBERG/2_ModesMatrix/Interpolated_L/lambdaIm.dat');
LAM=ReLAMBDA+1i*ImLAMBDA;
iter=1;
for i=1:NptsNew
   for j=1:NptsNew
       rd=RTDif(iter);
       rr=RTRad(iter,:);
       xR=ReLAMBDA(iter,:);
       xI=ImLAMBDA(iter,:);
       xi=xR+1i*xI;
       RT(i,j) = (rd + sum(xi.*rr));
       iter=iter+1;
   end
end
figure(1);
subplot(3,1,1);
pltphase(omegaNew,RT,'m');
xlim([2*pi/80,2*pi/15]);
ylim([-0.06,0.06]);
title('Transmission $T(\omega)$');

%%
iter=1;
for i=1:NptsNew
   for j=1:NptsNew
       rd=RCDif(iter);
       rr=RCRad(iter,:);
       xR=ReLAMBDA(iter,:);
       xI=ImLAMBDA(iter,:);
       xi=xR+1i*xI;
       RC(i,j) = (rd + sum(xi.*rr));
       iter=iter+1;
   end
end
subplot(3,1,2)
pltphase(omegaNew,RC,'m');
xlim([2*pi/80,2*pi/15]);
ylim([-0.06,0.06]);
title('Reflection $R(\omega)$');
