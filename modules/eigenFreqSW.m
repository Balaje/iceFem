function omegap = eigenFreqSW(L,H,th,Ad,guess)
%% Finding the Resonance frequencies.

%% Get all the properties
[~,~,~,~,E,nu,rhow,rhoi,g,~] = getProperties();

d = (rhoi/rhow)*th;

omegap = guess;
dw = 1e-3;
tol = 1e-7;

EI = E*th^3/(12*(1-nu^2));
Lc = (EI/(rhow*g))^(1/4);
tc = sqrt(rhow*Lc^6/(EI*H));

beta = 1;
gamma = (rhoi/rhow)*(th/Lc);

while (abs(dw) > tol)
    %% First Omega.
    T = 2*pi/omegap;
    Ap = (g/(1i*omegap))*Ad;
    LL = L/Lc; HH = H/Lc; dd = d/Lc; TT = T/tc;
    ndOmega = 2*pi/TT;
    alpha = HH*ndOmega^2;
    [~,~,~,~,Ax] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,Ap);
    
    %% Second Omega.
    T = 2*pi/(omegap+dw);
    Ap = (g/(1i*(omegap+dw)))*Ad;
    LL = L/Lc; HH = H/Lc; dd = d/Lc; TT = T/tc;
    ndOmega = 2*pi/TT;
    alpha = HH*ndOmega^2;
    [~,~,~,~,Axp] = shallowmovingplate(LL,HH,dd,alpha,beta,gamma,Ap);    
    
    %% Update Step
    AD = (Axp-Ax)/dw;    
    dwnew=eig(Ax,-AD);
    dw=min(dwnew); 
    omegap=omegap+dw;
end
end