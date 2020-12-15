%% Program to compute the interpolated system 

clc
clear 
close all

nthickness=10;
a1=2*pi*0.0001;
b1=2*pi*0.01;
omega=linspace(a1,b1,10);
T=2*pi./omega;

filePath='1_SIMPLE5/2_ModesMatrix/th0/';
[omegaNew,detH,conH]=interpolateFreq(a1,b1,omega,16,filePath,199,1);