%% Comparison of the Euler-Bernoulli-FEM solution with the Linear Elasticity Solution

close all
clear
clc
global ff

%% Get the Parameters of the ice.
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 800;
L = 20000;
omega = 2*pi/(200); % Wave Period can be complex.
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

%% FreeFem++ Part to solve the problem using LE method
file = 'iceSpline.edp';
ffpp=[ff,' -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(5),' -isUniCav ',num2str(0)];
[aa1,bb1]=system(cmd);
if(aa1)
    error('Cannot run program. Check path of FF++ or install it');
end

RC1 = load('1_Forced/2_RefCoeff/refcoeff0.dat');
RefIV = RC1(1)+1i*RC1(2);

%% FreeFem++ Part to solve the FEM-EulerBernoulli Problem
file = 'iceshelf_submerged_moving.edp';
ffpp=[ff,' -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(5),' -isUniCav ',num2str(0)];

[aa2,bb2]=system(cmd);
if(aa2)
    error('Cannot run program. Check path of FF++ or install it');
end

%% Load the files for Linear Elasticity Solution
[pts1,seg1,tri1] = importfilemesh('1_Forced/2_Deformation/movedIce0.msh');
[pts2,seg2,tri2] = importfilemesh('1_Forced/2_Deformation/cavityMesh.msh');
PHILE = importfiledata('1_Forced/2_Potential/potential0.bb');

%% Load the files for FEM-EulerBernoulli Solution
realDisp=load('1_ThinPlate/displacement_real0.txt');
uEBfem=flipud(realDisp);
[pts3,seg3,tri3] = importfilemesh('1_ThinPlate/cavityMesh.msh');
PHIEB = importfiledata('1_ThinPlate/potential0.bb');


%% Plot and visualize the data
fig=figure(1);
set(fig,'Position',[1,1,654,523]);
sp1=subplot(3,1,1);
pos1=get(sp1,'Position');
% Modal Expansion
pp1=pdeplot(pts1,seg1,tri1);
set(pp1,'Color','m');
axis equal
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
ylabel('$z$');
grid on
hold on
% Euler Bernoulli
plot(uEBfem(1:5:end,1),uEBfem(1:5:end,2),'bo','MarkerSize',6,'LineWidth',1.5);
xlim([0,LL]);
axis equal
ylabel('$\mathbf{w}(x,z)$');
hold on
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);

%% Plot the velocity potential
sp2=subplot(3,1,2);
pdeplot(pts2,seg2,tri2,'XYData',real(PHILE)','colormap','jet');
caxis([min(real(PHILE)),max(real(PHILE))]);
ylim([-HH,-dd]);
yticks([-HH,-(HH+dd)/2,-dd]);
yticklabels([-H,-(H+d)/2,-d]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
xlim([0,LL]);
ylabel('$z$');
grid on

sp3=subplot(3,1,3);
pdeplot(pts3,seg3,tri3,'XYData',real(PHIEB)','colormap','jet');
caxis([min(real(PHILE)),max(real(PHILE))]);
ylim([-HH,-dd]);
xlabel('$x$ (in km)');
ylabel('$z$');
grid on
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
yticks([-HH,-(HH+dd)/2,-dd]);
yticklabels([-H,-(H+d)/2,-d]);

pos2=get(sp2,'Position');
set(sp1,'Position',[pos2(1),pos1(2),pos2(3),pos2(4)]);

%% Save figure file
saveas(gcf,['femEBCompareT',num2str(T),'.fig'],'fig');
set(gcf, 'color', 'white');
export_fig(['femEBCompareT',num2str(T),'.pdf'],'-pdf','-transparent',gcf);