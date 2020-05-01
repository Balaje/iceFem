function CreatePaths
clc
close all

fprintf('>> Run:\n\n   which FreeFem++\n\n   in your command line to get the path for FreeFem++.\n   Set the full path in the variable `ff` in CreatePaths.m\n');
addpath([pwd,'/modules/']);
addpath([pwd,'/examples/']);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);
envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

%% Should be set manually by the user.
global ff
ff='/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++';

end
