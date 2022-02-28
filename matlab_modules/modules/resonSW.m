%% Phase Plot - Shallow Water.

clear
clc
close all

addpath('modules/');

set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',16);

%% Ice - Shelf parameters
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
L = 20000;
H = 500;
th = 200;
d = (rhoi/rhow)*th;
a = 2*pi/300;
b = 2*pi/30;
Tr = linspace(a,b,200);
Ti = linspace(-0.02,0.02,200);
[X,Y] = meshgrid(Tr,Ti);

omega = X+1i*Y;

T = 2*pi./omega;
Ad = 1;
Ap = (g./(1i*omega))*Ad;

RefSW = zeros(length(X),length(Y));

%% Compute the Reflection Coefficient
for m=1:size(T,1)
    for n=1:size(T,2)
        EI = E*th^3/(12*(1-nu^2));
        Lc = (EI/(rhow*g))^(1/4);
        tc = sqrt(rhow*Lc^6/(EI*H));
        
        LL = L/Lc; HH = H/Lc; dd = d/Lc; TT = T(m,n)/tc;
        ndOmega = 2*pi./TT;
        alpha = HH*ndOmega.^2;
        beta = 1;
        gamma = (rhoi/rhow)*(th/Lc);
        
        [c,R,kSw,~,A] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,Ap(m,n));
        RefSW(m,n) = R/Ap(m,n);
    end
end

%% Plot the reflection coefficients
figure(1);
subplot(2,1,1);
PP=pltphase(omega,RefSW);
ylabel('$\Im(\omega)$','Interpreter','latex')
xlabel('$\Re(\omega)$','Interpreter','latex')
hold on

%% Find the Eigenfrequency with a guess.
omegap=[0.02,0.03,0.035,0.04,0.05,0.06,0.08,0.1,0.14, 0.18];
omegap2=zeros(length(omegap),1);
for m=1:length(omegap)
    omegap2(m)=eigenFreqSW(L,H,th,Ad,omegap(m));
end

%% Plot the points
scatter(real(omegap2),imag(omegap2),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',1.5)
scatter(real(omegap2),-imag(omegap2),'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',1.5)
xlim([Tr(1),Tr(end)]);

%% Find the maximum displacement for a set of real frequencies.
x = linspace(0,LL,500);
Tr = linspace(2*pi/b,2*pi/a,1000);
Ush = zeros(length(x),1);
maxDisp = zeros(length(Tr),1);
omega = 2*pi./Tr;
Ap = (g./(1i*omega))*Ad;
for m=1:length(Tr)
    EI = E*th^3/(12*(1-nu^2));
    Lc = (EI/(rhow*g))^(1/4);
    tc = sqrt(rhow*Lc^6/(EI*H));
    
    LL = L/Lc; HH = H/Lc; dd = d/Lc; TT = Tr(m)/tc;
    ndOmega = 2*pi/TT;
    alpha = HH*ndOmega.^2;
    beta = 1;
    gamma = (rhoi/rhow)*(th/Lc);
    
    [c,R,k,kappa_s,A] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,Ap(m));
    
    u_sh = @(x) ((H-d)*1/(1i*omega(m)*Lc^2))*( c(1)*(kappa_s(1))^2*exp(kappa_s(1)*(x-LL)) + ...
        c(2)*(kappa_s(2))^2*exp(kappa_s(2)*(x-LL)) + ...
        sum(transpose(c(3:end)).*(kappa_s(3:6)).^2.*exp(kappa_s(3:6)*x)) );
    for j=1:length(x)        
        Ush(j) = u_sh(x(j));
    end
    maxDisp(m) = max(abs(Ush));
end

%% Plot the maximum displacement with frequency
subplot(2,1,2);
plot(omega, maxDisp);
xlim([a,b]);