function startup
clc
close all

fprintf('Run:\n which FreeFem++\n in your command line to get the path and set the full path in the variable\n ff\n in startup.m');
addpath('./modules/');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);
envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

%% Should be set manually by the user.
global ff
ff='/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++';

end