clear
clc

[~,~,~,~,E,nu,rhow,rhoi,g,~]=getProperties();
H=800;
L=10000;
th=200;
T=400;
omega=2*pi/T;
d=(rhoi/rhow)*th;
NModes = 10;
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
[a,b,kappa,k,~,~,~] = movingplate(LL,HH,dd,NModes,alpha,beta,gamma,Ap);
RefTP = a(1)/Ap;

u_fd = @(x) -(1/(1i*omega*Lc))*( sum(transpose(b(1:NModes+3)).*exp(-kappa.*x).*...
    (-kappa).*tan(kappa*(HH-dd))) + ...
    sum(transpose(b(NModes+4:end)).*exp(kappa.*(x-LL)).*...
    (-kappa).*tan(kappa*(HH-dd))) );
xpts = linspace(0,LL,100);
Ufd = zeros(1,length(xpts));
for m=1:length(xpts)
    Ufd(m)=real(u_fd(xpts(m)));
end

figure(1)
plot(xpts,Ufd,'linewidth',3);
xlim([xpts(1),xpts(end)]);
hold on
grid on