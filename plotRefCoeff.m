%% Program to study the vibration of the Brunt Ice--shelf.

clc
clear
close all
format long

[~,~,~,~,E,nu,rhow,rhoi,g,~]=getProperties();
a=2*pi*0.0001;
b=2*pi*0.01;
omega=linspace(a,b,200); %Solved in the HPC grid.
T=2*pi./omega;
Ad=1;
Ap=(g./(1i*omega))*Ad;


% Set the new-frequency space;
npts=499;
a1=a; b1=b;
omegaNew=linspace(a1,b1,npts+1);
Tnew=2*pi./omegaNew;
ApNew=(g./(1i*omegaNew))*Ad;
file1='7_BEDMAP2/';
rc=zeros(length(omega),1);
for m=1:length(omega)
    RC=load([file1,'2_RefCoeff/RefCoeff',num2str(m),'.dat']);
    rc(m)=RC(1)+1i*RC(2);
end

% Interpolate the reflection coefficient.
filePath=[file1,'2_RefCoeff'];
nev=64;
rd=zeros(length(omega),1);
rr=zeros(length(omega),nev);
lambdaj=zeros(nev,length(omega));
rrc=zeros(length(omega),1);
for m=1:length(omega)
    rcDiff = load([filePath,'/RefCoeff_Dif/refC',num2str(m),'.dat']);
    rcRad = load([filePath,'/RefCoeff_Rad/refC',num2str(m),'.dat']);    
    rd(m) = rcDiff(1)+1i*rcDiff(2);
    rr(m,:) = (rcRad(:,1)+1i*rcRad(:,2)).';
    
    lam = fileread([file1,'/2_ModesMatrix/lambdaj',num2str(m),'.dat']);
    lam=str2num(lam);
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
    rcNew(m)=rcNew(m)+(LAM(m,:)*rrNew(m,:).'); 
%     
%     % Plot the lambdaj chart
%     figure(10);
%     bar(1:length(LAM(m,:)), abs(LAM(m,:)));
%     ylim([0,200]);
%     pause(0.01);    
end

%%
figure(2);
subplot(2,1,1);
semilogx(Tnew,real(rcNew),'b',Tnew,imag(rcNew),'r','Linewidth',1.5);
legend('Re parts of R($\omega$)','Im parts of R($\omega$)')
xlim([Tnew(end),Tnew(1)]);
grid on


subplot(2,1,2);
UY=load([file1,'uyabs.dat']);
semilogx(T,UY(:,1),'Linewidth',1.5);
xlabel('Wave period $T$ (in s)')
grid on
hold on
% [val,indx]=find(T>300);
% UY1=UY(indx,1);
% [pks,loc]=findpeaks(UY1(:,1));
% plot(T(loc),UY1(loc),'o');
% 
% dlmwrite('waveperiods.txt',T(loc),'delimiter','\n','precision',8);

