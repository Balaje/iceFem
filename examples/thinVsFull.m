%% Comparison of the Euler-Bernoulli solution with the Linear Elasticity Solution

close all
clear
clc
global ff

my_path=which('thinVsFull.m');
my_path=my_path(1:end-13);

if ~strcmp(pwd,my_path)
 hld_path=pwd;
 cd(my_path)
end

%% Get the Parameters of the ice.
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 800;
L = 20000;
T=200;
omega = 2*pi/T; % Wave Period can be complex.
th = 200;
d = (rhoi/rhow)*th;

%% Get the potential from Euler Bernoulli solution
% Non Dimensionalise
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

% Euler Bernoulli solution
[a,b,kappa,k,~,~,~] = movingplate(LL,HH,dd,NModes,alpha,beta,gamma,Ap);
RefTP = a(1)/Ap;

%% FreeFem++ Part to solve the problem using invacuo modes
file = 'iceSpline.edp';
ffpp=[ff,' -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(5)];
[aa,bb1]=system(cmd);
if(aa)
    error('Cannot run program. Check path of FF++ or install it');
end

RC1 = load('1_Forced/2_RefCoeff/refCoeff0.dat');
RefIV = RC1(1)+1i*RC1(2);

%% FreeFem++ Part to solve the Direct FEM problem
file = 'icefem.edp';
ffpp=[ff,' -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(2)];
[aa,bb]=system(cmd);
if(aa)
    error('Cannot run program. Check path of FF++ or install it');
end

RC2 = load('1_DirectFEM/2_RefCoeff/refCoeff0.dat');
RefFEM = RC2(1)+1i*RC2(2);

%% Plot the results
[pts1,seg1,tri1] = importfilemesh('1_Forced/2_Deformation/movedIce0.msh');
[pts2,seg2,tri2] = importfilemesh('1_DirectFEM/2_Deformation/movedIce0.msh');

u_fd = @(x) -(1/(1i*omega*Lc))*( sum(transpose(b(1:NModes+3)).*exp(-kappa.*x).*...
    (-kappa).*tan(kappa*(HH-dd))) + ...
    sum(transpose(b(NModes+4:end)).*exp(kappa.*(x-LL)).*...
    (-kappa).*tan(kappa*(HH-dd))) );

xpts = linspace(0,LL,40);
Ufd = zeros(1,length(xpts));
for m=1:length(xpts)
    Ufd(m)=real(u_fd(xpts(m)));
end

%% Plot the data
fig=figure(1);
set(fig,'Position',[359   123   733   575]);

sp1=subplot(2,1,1);
pos1=get(sp1,'Position');

% ModalExpansion
pp1=pdeplot(pts1,seg1,tri1);
set(pp1,'Color','m');
axis equal
ylabel('$z$');
hold on
% EulerBernoulli
plot(xpts,Ufd,'bo','MarkerSize',6,'LineWidth',1.5);
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
axis equal
ylabel('$\mathbf{w}(x,z)$ [m]');
xlabel('$x$ [km]');
hold on
grid on

legend(sp1,'Full elasticity (modal method)','Thin plate');

% DirectFEM
sp2=subplot(2,1,2);
pos2=get(sp2,'Position');
pp2=pdeplot(pts2,seg2,tri2);
clb=[0.3010 0.7450 0.9330];
set(pp2,'Color',clb);
axis equal
xlabel('$x$ (in km)');
grid on
hold on
% EulerBernoulli
plot(xpts,Ufd,'bo','MarkerSize',6,'LineWidth',1.5);
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
axis equal
ylabel('$\mathbf{w}(x,z)$');
xlabel('$x$ (in km)');

legend(sp2,'Full elasticity (direct FEM)','Thin plate');

set(figure(1),'Name','Vertical deflection of shelf-cavity interface')

%% Compute the velocity potentials
[pts3,seg3,tri3] = importfilemesh('1_Forced/2_Deformation/cavityMesh.msh');
PHILE = importfiledata('1_Forced/2_Potential/potential0.bb');
PHIEB = zeros(1,length(pts3));
for m=1:NModes+3
    psi = @(z) cos(kappa(m)*(z+HH))/cos(kappa(m)*HH);    
    PHIEB = PHIEB + (b(m)*exp(-kappa(m)*pts3(1,:)).*psi(pts3(2,:))) + ...
        b(m+NModes+3)*exp(kappa(m)*(pts3(1,:)-LL)).*psi(pts3(2,:));
end

%% Plot the velocity potentials
fig2=figure(2);
set(fig2,'Position',[359   123   733   575]);

sp3=subplot(2,1,1);
pdeplot(pts3,seg3,tri3,'XYData',real(PHILE)','colormap','jet');
caxis([min(real(PHILE)),max(real(PHILE))]);
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
ylim([-HH,-dd]);
yticks([-HH,-(HH+dd)/2,-dd]);
yticklabels([-H,-(H+d)/2,-d]);
xlabel('$x$ [km]');
ylabel('$z$ [m]');
grid on

title(sp3,'Velocity potential Re($\phi(x,z)$) (modal method)')

sp4=subplot(2,1,2);
pdeplot(pts3,seg3,tri3,'XYData',real(PHIEB)','colormap','jet');
caxis([min(real(PHILE)),max(real(PHILE))]);
xlim([0,LL]);
xticks([0,LL/4,LL/2,3*LL/4,LL]);
xticklabels([0,L/4000,L/2000,3*L/4000,L/1000]);
ylim([-HH,-dd]);
yticks([-HH,-(HH+dd)/2,-dd]);
yticklabels([-H,-(H+d)/2,-d]);
ylabel('$z$');
xlabel('$x$ (in km)');
grid on

title(sp4,'Velocity potential Re($\phi(x,z)$) (thin plate)')

set(figure(2),'Name','Velocity potential in sub-shelf water cavity')

%% Save figure file (Uncomment to save)
% saveas(fig,['Compare1T',num2str(T),'.fig'],'fig');
% saveas(fig2,['Compare2T',num2str(T),'.fig'],'fig');
% filename1=['Compare1T',num2str(T),'.pdf'];
% filename2=['Compare2T',num2str(T),'.pdf'];
% 
% set(fig,'color','white');
% export_fig(filename1,'-pdf','-transparent',fig);
% set(fig2,'color','white');
% export_fig(filename2,'-pdf','-transparent',fig2);


%% Print the Uniform square mesh.
[ptsU,segU,triU] = importfilemesh('1_DirectFEM/2_Deformation/iceMesh.msh');
[ptsC,segC,triC] = importfilemesh('1_DirectFEM/2_Deformation/cavityMesh.msh');

fig10=figure(10);
set(fig10,'Position',[542   132   739   566]);
subplot(3,1,1);
pp1=pdeplot(ptsU,segU,triU);
set(pp1,'Color','r');
axis equal
ylabel('$z$');
grid on
hold on
pp2=pdeplot(ptsC,segC,triC);
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
pp1=pdeplot(ptsU,segU,triU);
set(pp1,'Color','r');
axis equal
ylabel('$z/L_c$');
grid on
hold on
pp2=pdeplot(ptsC,segC,triC);
set(pp2,'Color','b');
axis equal
xlabel('$x$ [km]');
ylabel('$z/L_c$');
xlim([0-LL/10,0+LL/5]);
xticks([0-LL/10,0+LL/5]);
xticklabels(1/1000*[0-L/10,0+L/5]);
grid on
hold on

subplot(3,1,3);
pp1=pdeplot(ptsU,segU,triU);
set(pp1,'Color','r');
ylabel('$z/L_c$');
grid on
hold on
pp2=pdeplot(ptsC,segC,triC);
set(pp2,'Color','b');
axis equal
xlabel('$x$ [km]');
ylabel('$z/L_c$');
xlim([LL-LL/5,LL]);
xticks([LL-LL/5,LL]);
xticklabels(1/1000*[L-L/5,L]);
grid on
hold on

set(fig10,'Name','Meshes')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('hld_path','var')
 cd(hld_path)
 clear hld_path
end

% Uncomment to write to a file.
%export_fig('meshesUni.pdf','-pdf','-transparent',fig10);