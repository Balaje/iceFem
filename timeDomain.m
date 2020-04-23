%% Program to generate the time domain solution.

clc
clear
close all

addpath('modules');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',16);

envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

%% Get the properties of the shelf.
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();
H = 800;
L = 20000;
omega = linspace(2*pi/(2000),2*pi/(20),201);
T = 2*pi./omega;
th = 200;
d = (rhoi/rhow)*th;
Ad = 1;
Ap = ((g./(1i*omega))*Ad);
EI = E*th^3/(12*(1-nu^2));
Lc = (EI/(rhow*g))^(1/4);
tc = sqrt(rhow*Lc^6/(EI*H));
LL = L/Lc; HH = H/Lc; dd = d/Lc;

%% Run the FreeFem++ code.
file = 'iceSpline.edp';
ffpp=['/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++ -nw -ne ', file];
for m=1:length(omega)
    cmd=[ffpp,' -Tr ',num2str(real(T)),' -Ti ',num2str(imag(T)),' -H ',num2str(H), ' -L ',num2str(L),' -h '...
    ,num2str(th),' -N ',num2str(3), ' -isUniIce ',num2str(0), ' -isUniCav ',num2str(0)];
    [aa,bb]=system(cmd);
    if(aa)
        error('Cannot run program. Check path of FF++ or install it');
    end
    fprintf('Finish m= %d\n',m);
end

%% All the wavenumbers
TT = T/tc;
ndOmega = 2*pi./TT;
alpha = HH*ndOmega.^2;

K = zeros(length(alpha),1);
for i=1:length(alpha)
    K(i) = fzero(@(x) tanh(H/Lc*x) - alpha(i)/x,[0.0001,100*alpha(i)]);
end
K_0 = 1;
fhat = @(K) (2/pi)*exp(-1*(K - K_0).^2);


%% Perform the Time-Domain conversion.
[pts,seg,tri] = importfilemesh('1_Forced/2_Deformation/iceMesh.msh');
xpts=pts(1,:);
ypts=pts(2,:);
ndofs=length(xpts);
UX = @(t) (0*zeros(1,ndofs)*t);
UY = @(t) (0*zeros(1,ndofs)*t);
w = @(x,t) 0*x*t;

for m=2:length(omega)
    % ----- IceDomain
    ReUX = importfiledata(['1_Forced/2_Deformation/ReUX',num2str(m),'.bb']);
    ImUX = importfiledata(['1_Forced/2_Deformation/ImUX',num2str(m),'.bb']);
    ReUY = importfiledata(['1_Forced/2_Deformation/ReUY',num2str(m),'.bb']);
    ImUY = importfiledata(['1_Forced/2_Deformation/ImUY',num2str(m),'.bb']);
  
    ux = ReUX+1i*ImUX;
    uy = ReUY+1i*ImUY;
    
    sumx = @(t) (1/pi*fhat(K(m))*(ux)*exp(-1i*omega(m)*t)*(K(m)-K(m-1)));
    sumy = @(t) (1/pi*fhat(K(m))*(uy)*exp(-1i*omega(m)*t)*(K(m)-K(m-1)));
    UX = @(t) (UX(t) + sumx(t));
    UY = @(t) (UY(t) + sumy(t));
    
    % ----- WaterDomain
    cr = load(['1_Forced/2_ModesMatrix/reC',num2str(m),'.dat']);
    ci = load(['1_Forced/2_ModesMatrix/imC',num2str(m),'.dat']);
    c = cr+1i*ci;
    NModes = length(c)-1;
    
    beta = 1;
    gamma = (rhoi/rhow)*(th/Lc);
    [a,b,kappa,k,~,~,~] = movingplate(LL,HH,dd,NModes,alpha(m),beta,gamma,Ap(m));
    k(1) = -k(1);
    
    xi = @(x) (1/(Lc*omega(m)*1i))*...
        ( Ap(m)*exp(-k(1)*x)*k(1)*tan(k(1)*HH) + ...
        sum(transpose(c).*exp(k*x).*k.*tan(k*HH)) );
    w = @(x,t) ( w(x,t) + 1/pi*fhat(K(m))*xi(x)*exp(-1i*omega(m)*t)*(K(m)-K(m-1)) );
    
    fprintf('End m=%d\n',m);
end

%% Make the movie
v = VideoWriter('timeDomain.avi');
v.FrameRate = 10;
open(v);

TAX = linspace(-200,2000,201);
X = linspace(-L/Lc,0,200);
W = zeros(length(X),1);
for m=1:length(TAX)
    NEWPTS = [xpts+real(UX(TAX(m))); ypts+real(UY(TAX(m)))];
                
    for n=1:length(X)
        W(n) = w(X(n),TAX(m));        
    end

    figure(1)
    pdeplot(NEWPTS,seg,tri);    
    hold on    
    plot(X,real(W),'LineWidth',1);
    hold off
    %title(['Time = ',num2str(TAX(m))]);    
    xlim([min(X),L/Lc])
    %xlim([0,L/Lc])
    ylim([-2,2])
    axis square
    
    if(TAX(m)==0)
       pause() 
    end
      
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
