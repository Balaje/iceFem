function CreatePaths
clc
close all

addpath([pwd,'/modules/']);
addpath([pwd,'/examples/']);
set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);
envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

if exist('~/bedmap2/','dir')
 addpath('~/bedmap2/')
else
 fprintf('>> Must addpath to bedmap2\n')
end

if exist('~/Documents/MATLAB/bedmap2_toolbox_v4.6.2/','dir')
 addpath('~/Documents/MATLAB/bedmap2_toolbox_v4.6.2/')
else
 fprintf('>> Must addpath to bedmap2_toolbox*\n')
end

if exist('~/Library/Application Support/MathWorks/MATLAB Add-Ons/Toolboxes/Antarctic Mapping Tools/AntarcticMappingTools/','dir')
 addpath('~/Library/Application Support/MathWorks/MATLAB Add-Ons/Toolboxes/Antarctic Mapping Tools/AntarcticMappingTools/')
else
 fprintf('>> Must addpath to Antarctic Mapping Tools\n')
end

%% Should be set manually by the user.
global ff
if exist('/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++')
 ff='/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++';
elseif exist('/usr/local/ff++/mpich3/bin/FreeFem++')
 ff='/usr/local/ff++/mpich3/bin/FreeFem++';
else
 fprintf('>> Run:\n\n   which FreeFem++\n\n   in your command line to get the path for FreeFem++.\n   Set the full path in the variable `ff` in CreatePaths.m\n');
end

end
