%% Example program to illustrate running FreeFem++ and plotting with MATLAB

clc
clear
close all
global ff

%% Parameters
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 800;
L = 20000;
omega = 2*pi/(300); % Wave Period can be complex.
T = 2*pi/omega;
th = 200;
d = (rhoi/rhow)*th;

%% Run the code;
file = 'simple1.edp';
ffpp=[ff,' -nw -ne ', file];
cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(5), ' -isUniIce ',num2str(0), ' -isUniCav ',num2str(0)];
[aa,bb1]=system(cmd);
if(aa)
    error('Cannot run program. Check path of FF++ or install it');
end

%% From the FreeFem to MATLAB converter.
[pts1,seg1,tri1] = importfilemesh('meshIce0.msh');
[pts2,seg2,tri2] = importfilemesh('meshCav0.msh');
UX = importfiledata('xDisp0.bb');
UY = importfiledata('yDisp0.bb');
PHILE = importfiledata('potentialCav0.bb');

% Define new pts for ice.
pts3=zeros(size(pts1));
pts3(1,:)=pts1(1,:)+UX;
pts3(2,:)=pts1(2,:)+UY;

% Plot using PDEPLOT
figure(1)
subplot(2,2,[1,2])
pdeplot(pts1,seg1,tri1)
axis equal
grid on
xlabel('$x$');
ylabel('$z$');
subplot(2,2,3);
pdeplot(pts1,seg1,tri1)
axis equal
grid on
xlim([-0.5,1]);
xlabel('$x$');
ylabel('$z$');
subplot(2,2,4);
pdeplot(pts1,seg1,tri1)
axis equal
grid on
xlim([19,21])
xlabel('$x$');
ylabel('$z$');

figure(2)
subplot(2,1,1);
pdeplot(pts3,seg1,tri1);
axis equal
grid on
xlabel('$x$')
ylabel('$z$')

subplot(2,1,2);
pdeplot(pts2,seg2,tri2,'XYData',real(PHILE)','colormap','jet');
grid on
xlabel('$x$')
ylabel('$z$')