%% Program to plot the response and cavity geometries

close all
clear
clc

addpath('./modules/');

set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);

envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

%% Get the Parameters of the ice.
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 800;
L = 20000;
omega = 2*pi/(300); % Wave Period can be complex.
T = 2*pi/omega;
th = 200;
d = (rhoi/rhow)*th;

%% Non Dimensionalise
NModes = 5;
EI = E*th^3/(12*(1-nu^2));
Lc = (EI/(rhow*g))^(1/4);
tc = sqrt(rhow*Lc^6/(EI*H));
LL = L/Lc; HH = H/Lc; dd = d/Lc; TT = T/tc;
ndOmega = 2*pi./TT;
alpha = HH*ndOmega.^2;
beta = 1;
gamma = (rhoi/rhow)*(th/Lc);
Ad = 1;
Ap = ((g./(1i*omega))*Ad);

%% FreeFem++ bit to solve the Linear Elasticity problem.
file = 'iceSpline.edp';
ffpp=['/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++ -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(3), ' -isUniIce ',num2str(0), ' -isUniCav ',num2str(0)];
[aa,bb1]=system(cmd);
if(aa)
    error('Cannot run program. Check path of FF++ or install it');
end

%% Plot the results
[pts0,seg0,tri0] = importfilemesh('1_Forced/2_Deformation/iceMesh.msh');
[pts1,seg1,tri1] = importfilemesh('1_Forced/2_Deformation/movedIce0.msh');
[pts2,seg2,tri2] = importfilemesh('1_Forced/2_Deformation/cavityMesh.msh');
PHILE = importfiledata('1_Forced/2_Potential/potential0.bb');


fig1=figure(1);
set(fig1,'Position',[542   132   739   566]);
subplot(3,1,1);
pp1=pdeplot(pts0,seg0,tri0);
set(pp1,'Color','r');
axis equal
ylabel('$z$');
grid on
hold on
pp2=pdeplot(pts2,seg2,tri2);
set(pp2,'Color','b');
axis equal
xlabel('$x$ (in km)');
ylabel('$z/L_c$');
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
grid on
hold on

subplot(3,1,2);
pp1=pdeplot(pts0,seg0,tri0);
set(pp1,'Color','r');
axis equal
ylabel('$z/L_c$');
grid on
hold on
pp2=pdeplot(pts2,seg2,tri2);
set(pp2,'Color','b');
axis equal
xlabel('$x$ (in km)');
ylabel('$z/L_c$');
xlim([0-LL/10,0+LL/5]);
xticks([0-LL/10,0+LL/5]);
xticklabels(1/1000*[0-L/10,0+L/5]);
grid on
hold on

subplot(3,1,3);
pp1=pdeplot(pts0,seg0,tri0);
set(pp1,'Color','r');
ylabel('$z/L_c$');
grid on
hold on
pp2=pdeplot(pts2,seg2,tri2);
set(pp2,'Color','b');
axis equal
xlabel('$x$ (in km)');
ylabel('$z/L_c$');
xlim([LL-LL/5,LL]);
xticks([LL-LL/5,LL]);
xticklabels(1/1000*[L-L/5,L]);
grid on
hold on

export_fig('meshes.pdf','-pdf','-transparent',fig1);

%% Plot the Potential and the response
fig2=figure(2);
set(fig2,'Position',[1,1,654,523]);
sp2=subplot(2,1,1);
pos2=get(sp2,'Position');
pp3=pdeplot(pts1,seg1,tri1);
set(pp3,'Color','r');
xlim([0,LL]);
ylabel('$\mathbf{w}(x,z)$');
grid on
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);


sp3=subplot(2,1,2);
pdeplot(pts2,seg2,tri2,'XYData',real(PHILE)','colormap','jet');
caxis([min(real(PHILE)),max(real(PHILE))]);
xlim([0,LL]);
xlabel('$x$ (in km)');
ylabel('$z$');
yticks([-HH,-(HH+dd)/2,-dd]);
yticklabels([-H,-(H+d)/2,-d]);
grid on
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);


%% Save figure file
saveas(fig2,['femLEnon',num2str(T),'.fig'],'fig');
set(fig2, 'color', 'white');
export_fig(['femLEnon',num2str(T),'.pdf'],'-pdf','-transparent',fig2);