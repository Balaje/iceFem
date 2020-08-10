%% Program to interpolate the linear system and solve to obtain the new lambdaj

clc
clear 
close all

[~,~,~,~,E,nu,rhow,rhoi,g,~]=getProperties();
a=100;
b=1000;
T=linspace(a,b,50); %Solved in the HPC grid.
omega=2*pi./T;
Ad=1;
Ap=(g./(1i*omega))*Ad;

%% Define the new frequency space.
npts=1000;
a1=400; b1=1000;
TNew=linspace(a1,b1,npts+1);
omegaNew=2*pi./TNew;
ApNew=(g./(1i*omegaNew))*Ad;
file1='1_BEDMAP2/';

%% Interpolate the system, solve and write the solution.
filePath = [file1,'2_ModesMatrix/'];
nev=64;
[omegaNew,detH,condH] = interpolateFreq(a1,b1,omega,nev,filePath,npts-1,1);
