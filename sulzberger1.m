%% Program to obtain the fourier transform of the incident wave from the NASA data.

clc
clear
close all

A=readtable('ROB002011.csv');
B=table2array(A(:,3));
Amp=B(:,1);
DT=A.Var2;
refTime=datenum(DT(1));
t=(datenum(DT)-refTime)*24*60*60; %In Minutes/Hours.
interval=t(2)-t(1);
Fs=1/interval;
%Amp1=smoothdata(Amp,'rlowess');
[Amp11,Amp12]=envelope(Amp,30,'peak');
Amp1=0.5*(Amp11+Amp12);
[t,Amp,Amp1]=freqRange(1,length(t),t,Amp,Amp1);
L=length(Amp);

%% Plot the signal and the difference.
figure(1);
subplot(2,1,1);
plot(A.Var2,Amp,A.Var2,Amp1);
xlabel('Time in Months');
ylabel('Displacement');
legend('Original Signal','Tide Component (Smoothed)');
subplot(2,1,2);
Amp2=Amp-Amp1;
plot(A.Var2,Amp2);
xlabel('Time in Months');
ylabel('Difference');

%% Perform band-pass filter on the difference.
figure(2);
bandpass(Amp2,[0.0001,0.001],Fs);
[Y,d1]=bandpass(Amp2,[0.0001,0.001],Fs);

%% Perform Fourier Transform to find the amplitude of the components.
figure(3);
YY=fft(Y);
P2=abs(YY/L);
% Looking at [0,L/2] due to symmetry
P1=P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f=Fs*(0:1/L:(1/2));
% Plot the signal again.
subplot(2,1,1);
plot(t,Y);
xlabel('t (In seconds)');
ylabel('Displacement (m)');
xlim([t(1),t(end)]);
hold on
% Plot the FFT.
subplot(2,1,2);
plot(1000*f,P1);
xlabel('Frequency (mHz)');
ylabel('Amplitude (m)');
xlim([1000*f(1),1000*f(end)]);
hold on
