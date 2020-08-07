%% Program to interpolate the linear system and solve to obtain the new lambdaj

clc
clear 
close all

[~,~,~,~,E,nu,rhow,rhoi,g,~]=getProperties();
a=0.0001;
b=0.001;
omega=linspace(a,b,51); %Solved in the HPC grid.
T=2*pi./omega;
Ad=1;
Ap=(g./(1i*omega))*Ad;

%% Define the new frequency space.
npts=300;
a1=a; b1=b;
omegaNew=linspace(a1,b1,npts+1);
Tnew=2*pi./omegaNew;
ApNew=(g./(1i*omegaNew))*Ad;
file1='1_BEDMAP2/';

%% Interpolate the system, solve and write the solution.
filePath = [file1,'2_ModesMatrix/'];
nev=64;
[omegaNew,detH,condH] = interpolateFreq(a1,b1,omega,nev,filePath,npts-1,1);
