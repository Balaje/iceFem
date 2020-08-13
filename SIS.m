%% Program to obtain the fourier transform of the incident wave from the NASA data.

clc
clear

A=readtable('ROB002011.csv');
B=table2array(A(:,3));
Amp=B(:,1);
DT=A.Var2;
refTime=datenum(DT(1));
t=(datenum(DT)-refTime)*24*60*60; %In Minutes/Hours.
interval=t(2)-t(1);
Fs=1/interval;

Amp1=sgolayfilt(Amp,2,41);
[t,Amp,Amp1]=freqRange(1,length(t),t,Amp,Amp1);
L=length(Amp);
Amp2=Amp-Amp1;

%% Perform band-pass filter on the difference.
[Y,d1]=highpass(Amp2,0.0001,Fs);

%% Perform Fourier Transform to find the amplitude of the components.
YY=fft(Y);
P2=abs(YY/L);
% Looking at [0,L/2] due to symmetry
P1=P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f=Fs*(0:1/L:(1/2));

%% Linear regression in log scale.
f1=f(find(f>=0.0001));
P2=P1(find(f>=0.0001));
X=log10(f1); Y1=log10(P2);
bpoly=polyfit(X',Y1,1);

%% Original frequency space in the HPC
a1=2*pi*0.0001; 
b1=2*pi*0.01;
omega=2*pi*linspace(a1,b1,200); % Original frequency space in HPC
file1='1_BEDMAP2/';
% Interpolate the system, solve and write the solution.
filePath = [file1,'2_ModesMatrix/'];
nev=64;

%% New frequency space.
a=a1; 
b=b1;
npts=200;
XX=linspace(a,b,npts);
X1=log10(XX);
yhat=10.^(polyval(bpoly,X1));
% Write the frequency, amplitude and a random phase of the signal
phase=2*pi*rand(1,length(X1));
FAmp=[10.^(X1); yhat; phase];
dlmwrite('FAmp.dat',FAmp,'delimiter','\t','precision',16);

%% Interpolate the linear system.
[omegaNew,detH,condH] = interpolateFreq(a,b,omega,nev,filePath,npts-1,1);
