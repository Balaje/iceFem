function [a,b,kappa,k,chi,LHS,D1] = movingplate(L,h,D,NModes,alpha,beta,gamma,Ap)
% Program to solve the moving ice shelf problem with draft
%   alpha = (omega^2/g,  beta = D/(rho_w g),  gamma = (rho_i h)/(rho_w)
%   INPUT:  L, H, d, NModes
%   OUTPUT: a,b,kappa,k
%tic
M = NModes;
H = h; LL = L; d = D;

%% Get the roots of the dispersion equation
k = dispersion_free_surface(alpha,M,H);
kappa = dispersion_elastic_surface(alpha,beta,gamma,M+2,H-d);
chi = (0:M)*(pi/(H-d));
lambda = k; lambda(1) = -k(1);

% %% Get the matrices
% D1 = innerproduct(chi,kappa,H,d,[H-d,H-d]);
% D2 = innerproduct(chi,k,H,d,[H-d,H]);
% D3 = innerproduct(k,kappa,H,d,[H,H-d]);
D1 = zeros(M+1,M+3);
D2 = zeros(M+1,M+1);
D3 = zeros(M+1,M+3);

for m=1:M+1
   for n=1:M+3
      D1(m,n) = innerproduct(chi(m),kappa(n),H,d)/...
                (cos(chi(m)*(H-d))*cos(kappa(n)*(H-d))); 
      if(n<=M+1)
        D2(m,n) = innerproduct(chi(m),k(n),H,d)/...
                (cos(chi(m)*(H-d))*cos(k(n)*H));
      end
      D3(m,n) = innerproduct(k(m),kappa(n),H,d)/...
                (cos(k(m)*H)*cos(kappa(n)*(H-d)));
   end
end
%%% Amplitude depends on \alpha = H/Lc
A = diag(0.5*(cos(k*H).*sin(k*H) + k*H)./(k.*(cos(k*H)).^2)); %H=1;
B1 = [-D1, -D1.*repmat(exp(-kappa*LL),[M+1,1])];
B2 = [D3.*repmat(kappa,[M+1,1]), D3.*repmat(-(kappa).*exp(-kappa*LL),[M+1,1])];
B3 = [D1.*repmat((kappa).*exp(-kappa*LL),[M+1,1]), -D1.*repmat(kappa,[M+1,1])];

B4 = [-kappa.*exp(-kappa*LL).*tan(kappa*(H-d)); ...
       (kappa.^2).*exp(-kappa*LL).*tan(kappa*(H-d)); ...
      -(kappa.^3).*tan(kappa*(H-d)); ...
       (kappa.^4).*tan(kappa*(H-d))];

B5 = [-kappa.*tan(kappa*(H-d)); ...
       -(kappa.^2).*tan(kappa*(H-d)); ...
       -(kappa.^3).*tan(kappa*(H-d)).*exp(-kappa*LL); ...
       -(kappa.^4).*tan(kappa*(H-d)).*exp(-kappa*LL)]; 
   
f1 = -Ap*D2(:,1); %Amplitude depends on h
f2 = zeros(M+1,1);  f2(1) = -k(1)*Ap*A(1,1); %Amplitude depends on h

LHS = [D2,B1; ...       
       diag(lambda)*A,B2; ...       
       zeros(M+1),B3; ...
       zeros(4,M+1),B4,B5];
RHS = [f1; f2; zeros(M+1,1); zeros(4,1)];

% lA = diag(lambda)*A;
% LHS = [D2, B1; ...
%        zeros(M+1), D2\B1 - (lA\eye(M+1))*B2;...
%        zeros(M+1), B3; ...
%        zeros(4,M+1), B4, B5];
% RHS = [f1; D2\f1 - (lA\eye(M+1))*f2; zeros(M+5,1)];

sol = LHS\RHS;
a = sol(1:M+1);
b = sol(M+2:end);

%%% Plot routine - Uncomment to check%%%
% phi = @(z) cos(k*(z+H))./cos(k*H);
% psi = @(z) cos(kappa*(z+H))./cos(kappa*(H-d));
% 
% PHI1 = @(x,z) Ap*exp(k(1)*x)*cos(k(1)*(z+H))/cos(k(1)*H) + ...
%     sum(transpose(a).*exp(lambda*(x)).*phi(z));
% PHI2 = @(x,z) sum(transpose(b(1:M+3)).*exp(-kappa*(x)).*psi(z)) + ...
%     sum(transpose(b(M+4:end)).*exp(kappa*(x-LL)).*psi(z));
% % Begin plotting
% z_0 = linspace(-h,-D,30);
% phi1 = zeros(length(z_0),1);
% phi2 = zeros(length(z_0),1);
% for n=1:length(z_0)
%     phi1(n) = PHI1(0,z_0(n));
%     phi2(n) = PHI2(0,z_0(n));
% end
% figure(3)
% plot(z_0,abs(phi1),'r-+',z_0,abs(phi2),'b-')
% pause()
% toc
end

function V = innerproduct(k, kappa, H, d)
% if k and kappa are not close
if(abs(k - kappa) >= 1e-7)
    V = ( kappa*sin(kappa*(H-d))*cos(k*(H-d)) - ...
          k*cos(kappa*(H-d))*sin(k*(H-d)) )/(kappa^2-k^2);
else
    V = (H-d)/2 + sin(2*k*(H-d))/(4*k);
end

% [K, KAPPA] = meshgrid(k, kappa);
% V = (K.*sin((H-d)*K).*cos((H-d)*KAPPA) - ...
%         KAPPA.*cos((H-d)*K).*sin((H-d)*KAPPA))./...
%         ((K.^2 - KAPPA.^2).*cos(K*(factor(1))).*cos(KAPPA*factor(2)));
% V = transpose(V);
end