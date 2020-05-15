%% Program to work on generating the analytic extension plot.

clear
close all
global ff

%% Get the properties of the shelf.
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 500;
L = 20000;
th = 200;

a = 2*pi/300;
b = 2*pi/50;
[X,Y] = meshgrid(linspace(a,b,20),linspace(-0.02,0.02,20));
omega = X+1i*Y;

T = 2*pi./omega;
d = (rhoi/rhow)*th;
Ad = 1;
Ap = ((g./(1i*omega))*Ad);
EI = E*th^3/(12*(1-nu^2));
Lc = (EI/(rhow*g))^(1/4);
tc = sqrt(rhow*Lc^6/(EI*H));
LL = L/Lc; HH = H/Lc; dd = d/Lc;

%% Set the New Frequency Space
npts=300;
a1 = a;
b1 = b;
[X,Y] = meshgrid(linspace(a1,b1,npts+1),linspace(-0.02,0.02,npts+1));
omegaNew = X+1i*Y;
ApNew = (g./(1i*omegaNew))*Ad;

%% Define Reflection Coefficient Matrix.
rc = zeros(size(omega,1),size(omega,2));

%% Run the FreeFem++ code and obtain the reflection coefficients.
file = 'iceSpline.edp';
file1 = '1_Forced/';
ffpp=[ff,' -nw -ne ', file];
iter=1;
for m=1:size(omega,1)
    for n=1:size(omega,2)
%         cmd=[ffpp,' -Tr ',num2str(real(T(m,n))),' -Ti ',num2str(imag(T(m,n))),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
%             ,num2str(th),' -N ',num2str(4), ' -iter ', num2str(iter),' -isUni ',num2str(1)];
%         [aa,bb]=system(cmd);
%         if(aa)
%             error('Cannot run program. Check path of FF++ or install it');
%         end
        
        fprintf('Finish (m,n)= (%d,%d)\n',m,n);
        RC = load([file1, '2_RefCoeff/refCoeff',num2str(iter),'.dat']);
        rc(m,n) = RC(1)+1i*RC(2);
        iter = iter+1;
    end
end

%% Get the Reflection coefficients in the new frequency space.
filePath = [file1,'2_RefCoeff'];
nev=20;
flag1 = interpolateRefCoeff(omega,omegaNew,nev,filePath);

RefDifRe = load([filePath,'/Interpolated_R/refCDifRe.dat']);
RefDifIm = load([filePath,'/Interpolated_R/refCDifIm.dat']);
RefRadRe = load([filePath,'/Interpolated_R/refCRadRe.dat']);
RefRadIm = load([filePath,'/Interpolated_R/refCRadIm.dat']);

disp(flag1);

%% Interpolate the mode contributions
filePath = [file1,'2_ModesMatrix'];
flag2 = interpolateFreqComplex(omega,omegaNew,nev,filePath);

disp(flag2);

%% Loop over omegaNew
iter = 1;
xiR = load([filePath,'/Interpolated_L/lambdaRe.dat']);
xiI = load([filePath,'/Interpolated_L/lambdaIm.dat']);
rNew = zeros(size(omegaNew));
for m=1:size(omegaNew,1)
    for n=1:size(omegaNew,1)
        rd = RefDifRe(iter) + 1i*RefDifIm(iter);
        rr = RefRadRe(iter,:) + 1i*RefRadIm(iter,:);
        xR = xiR(iter,:); xI = xiI(iter,:);
        xi = xR + 1i*xI;
        rNew(m,n) = (rd + sum(xi.*rr))/ApNew(m,n);
        iter = iter+1;
    end
end

%% Plot the reflection coeffcients on the complex plane
fig1=figure(1);
set(fig1,'Position',[360,72,798,626]);
subplot(3,1,1);
pltphase(omega,rc);
pts1=num2str(length(omega));
title(['Before Interpolation for ',pts1,'x',pts1,' Grid']);
ylabel('$\Im(\omega)$','Interpreter','latex')
%xlabel('$\Re(\omega)$','Interpreter','latex')
xlim([a,b])


subplot(3,1,2);
pltphase(omegaNew,rNew);
pts2=num2str(length(omegaNew));
title(['After Interpolation for ',pts2,'x',pts2,' Grid']);
ylabel('$\Im(\omega)$','Interpreter','latex')
%xlabel('$\Re(\omega)$','Interpreter','latex')
xlim([a,b])
hold on

pause(0.01)


%% Compute the resonance frequencies in the complex plane.
guesses = [0.02,0.03,0.04,0.05,0.06,0.08,0.1];
roots = zeros(length(guesses),1);
for m=1:length(guesses)
    guess=guesses(m);
    roots(m) = findResonanceCplx(L,H,th,guess,file,'1_Forced/2_ModesMatrix/');
end

%% Mark the root on the complex plane
subplot(3,1,2);
scatter(real(roots),imag(roots),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
hold on
scatter(real(roots),-imag(roots),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);

%% %%%%%%%%%%%%% Run the Shallow water code to check the result.
RefSW = zeros(size(omegaNew));
T = 2*pi./omegaNew;

% EI = E*th^3/(12*(1-nu^2));
% Lc = (EI/(rhow*g))^(1/4);
% tc = sqrt(rhow*Lc^6/(EI*H));
% LL = L/Lc; HH = H/Lc; dd = d/Lc;

for m=1:size(omegaNew,1)
    for n=1:size(omegaNew,2)
        TT = T(m,n)/tc;
        ndOmega = 2*pi./TT;
        alpha = HH*ndOmega.^2;
        beta = 1;
        gamma = (rhoi/rhow)*(th/Lc);
        
        % Shallow water solution.
        [c,R,kSw,~,A] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,ApNew(m,n));
        RefSW(m,n) = R/ApNew(m,n);
    end
end
subplot(3,1,3);
pltphase(omegaNew,RefSW);
title('Shallow Water Solution');
ylabel('$\Im(\omega)$','Interpreter','latex')
xlabel('$\Re(\omega)$','Interpreter','latex')
xlim([a,b])
hold on

guesses = [0.02,0.03,0.04,0.05,0.06,0.08,0.1];
rootsSW = zeros(length(guesses),1);
for m=1:length(guesses)
    rootsSW(m) = eigenFreqSW(L,H,th,Ad,guesses(m));
end

subplot(3,1,3);
scatter(real(rootsSW),imag(rootsSW),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
hold on
scatter(real(rootsSW),-imag(rootsSW),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);

%-Find the maximum displacement with frequency
omegaSW = linspace(a1,b1,npts+1);

TSW=2*pi./omegaSW;
ApSW=(g./(1i*omegaSW))*Ad;
maxdisp_sh = zeros(length(TSW),1);
x=linspace(0,LL,100);
for m=1:length(TSW)
    TTSW = TSW(m)/tc;
    ndOmegaSW = 2*pi/TTSW;
    alpha = HH*ndOmegaSW^2;
    [c,R,~,kSw,A] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,ApSW(m));
    u_sh = @(x) ((H-d)*1/(1i*omegaSW(m)*Lc^2))*( c(1)*(kSw(1))^2*exp(kSw(1)*(x-LL)) + ...
        c(2)*(kSw(2))^2*exp(kSw(2)*(x-LL)) + ...
        sum(transpose(c(3:end)).*(kSw(3:6)).^2.*exp(kSw(3:6)*x)) );
    
    U_fd = zeros(length(x),1);
    U_sh = zeros(length(x),1);
    for p=1:length(x)
        U_sh(p) = u_sh(x(p));
    end
    
    maxdisp_sh(m) = max(abs(U_sh/Ad));
end


%% Plot the analytic extension against the resonance modes.
%{ 
(To test only)
lmdm = load('lmdmLambda.dat');
maxDisp = load('maxDispLambda.dat');
maxlocs = load('maxlocs.dat');
omegaNewReal = load('omegaNewLambda.dat');

fig2=figure(2,'Visible','off');
set(fig2,'Position',[360,72,798,626]);
s1=subplot(3,1,1);
pltphase(omegaNew,rNew);
pts2=num2str(length(omegaNew));
title(['After Interpolation for ',pts2,'x',pts2,' Grid']);
ylabel('$\Im(\omega)$','Interpreter','latex')
xlabel('$\Re(\omega)$','Interpreter','latex')
hold on
scatter(real(roots),imag(roots),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
hold on
scatter(real(roots),-imag(roots),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
xlim([a,b]);
s2=subplot(3,1,2);
plot(omegaNewReal,maxDisp,'linewidth',1);
[ps2,ls2,w2,p2] = findpeaks(maxDisp);
hold on
scatter(omegaNewReal(ls2),maxDisp(ls2),'o');
grid on
xlabel('$\omega$')
ylabel('$\max|u|$')
xlim([a,b]);
s3=subplot(3,1,3);
plot(omegaNewReal,lmdm,'linewidth',1);
[ps3,ls3,w3,p3] = findpeaks(lmdm);
hold on
scatter(omegaNewReal(ls3),lmdm(ls3),'o');
grid on
xlabel('$\omega$');
ylabel('$|\lambda_{14}|$');
xlim([a,b]);

%--> Same for Shallow water solution
fig3=figure(3,'Visible','off');
title('Resonance for Shallow Water Solution');
set(fig3,'Position',[360,72,798,626]);
subplot(2,1,1);
pltphase(omegaNew,RefSW);
title('Shallow Water Solution');
ylabel('$\Im(\omega)$','Interpreter','latex')
xlabel('$\Re(\omega)$','Interpreter','latex')
xlim([a,b])
hold on
scatter(real(rootsSW),imag(rootsSW),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
hold on
scatter(real(rootsSW),-imag(rootsSW),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[1 1 1],...
    'LineWidth',1.5);
subplot(2,1,2);
plot(omegaSW,maxdisp_sh,'linewidth',1);
xlim([omegaSW(1),omegaSW(end)]);
[ps3,ls3,w3,p3] = findpeaks(maxdisp_sh);
hold on
scatter(omegaSW(ls3),maxdisp_sh(ls3),'o');
xlim([omegaSW(1),omegaSW(end)]);
grid on

%}
%% Save all the figures. (Uncomment to save)
% saveas(fig1,'out1.png','png');
% saveas(fig2,'out2.png','png');
% saveas(fig3,'out3,png','png');
