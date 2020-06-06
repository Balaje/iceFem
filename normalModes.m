%% Interpolate the matrices;

clc
clear
global ff
close all

%% Parameters
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 500;
L = 20000;
omega = linspace(2*pi/300, 2*pi/50, 20);
T = 2*pi./omega;
th = 200;
d = (rhoi/rhow)*th;


%% Run the code to solve the problem.
file = 'simple2.edp';
ffpp=[ff,' -nw -ne ', file];

for m=1:length(omega)
    cmd=[ffpp,' -iter ',num2str(m),' -Tr ',num2str(real(T(m))),' -Ti ',num2str(imag(T(m))),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
        ,num2str(th),' -N ',num2str(5), ' -isUniIce ',num2str(1), ' -isUniCav ',num2str(1)];
    [aa,bb1]=system(cmd);
    if(aa)
        error('Cannot run program. Check path of FF++ or install it');
    end
    fprintf('Finish m = %d\n',m);
end

%% Compute the determinant of H
filePath='1_EigenModes/2_ModesMatrix/';
[omegaNew,detH,condH] = interpolateFreq(omega(1),omega(end),omega,20,filePath,1000,0);

open('refCoefQ.fig');
hold on
subplot(2,1,2);
semilogy(omegaNew,condH);
xlim([2*pi/300, 2*pi/50])

%% Find the resonance frequency by computing the peaks in the condition number,
[pks,loc] = findpeaks(condH);
hold on
scatter(omegaNew(loc),condH(loc));

ReHNew = load([filePath, 'Interpolated_H/ReH.dat']);
ImHNew = load([filePath, 'Interpolated_H/ImH.dat']);
Hmats = ReHNew(:,loc) + 1i*ImHNew(:,loc);

%% Plot the mode shapes.
[pts,seg,tri] = importfilemesh('1_EigenModes/2_ModesMatrix/iceMesh.msh');
VX = zeros(20,length(pts));
VY = zeros(20,length(pts));
for m=1:20
    % Load the m-th mode shape and store in an array.
    VX(m,:)=importfiledata(['1_EigenModes/2_Modes/modex',num2str(m-1),'.bb']);
    VY(m,:)=importfiledata(['1_EigenModes/2_Modes/modey',num2str(m-1),'.bb']);        
end
for n=1:length(loc)
    H=reshape(Hmats(:,n),20,20);
    [V,D]=eig(H,'vector');
    lam=V(:,(abs(D)==min(abs(D))));
    MODESHAPE=zeros(2,length(pts));    
    for m=1:20
        MODESHAPE(1,:) = MODESHAPE(1,:) + lam(m)*VX(m,:);        
        MODESHAPE(2,:) = MODESHAPE(2,:) + lam(m)*VY(m,:);
    end
    scale=100;
    DISP=pts+scale*MODESHAPE;
    figure(2)
    subplot(3,2,n);
    pdeplot(abs(DISP),seg,tri);
    grid on
    pause();    
end