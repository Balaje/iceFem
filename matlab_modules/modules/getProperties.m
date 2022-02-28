function [L, H, th, d, E, nu, rhow, rhoi, g, Ad] = getProperties()
%% Returns the default values.

L = 40000;
H = [300,350,400,450];
th = 20;
rhow = 1025;
rhoi = 922.5;
d = rhoi/rhow*th;

E = 2*10^9;
nu = 0.33;
g = 9.80665;

Ad = 0.6;
end