function [c,R,k,m,A] = shallowmovingplate(L,h,D,alpha,beta,gamma,Ap)
%%% Program to solve the shallow water model
%%% Input:  Length (L), depth of ocean (h), depth of cavity (hs)
%%%         submerged depth (D)?, alpha, beta, gamma


hs = h - D;
H = h;   Hs = hs;   LL = L;

k = sqrt(alpha/H);

poly = [beta, 0, 0, 0, (1-gamma*alpha), 0, alpha/Hs];
mm = roots(poly);

[~,index] = sort(real(mm)>1e-7,'descend');
m = mm(index);

m = transpose(m);
A = zeros(7,7);

A(1,1:2) = m(1:2);  
A(1,3:6) = m(3:6).*exp(m(3:6)*LL);

A(2,1:2) = (m(1:2)).^2; 
A(2,3:6) = (m(3:6)).^2.*exp(m(3:6)*LL);

A(3,1:2) = (m(1:2)).^3; 
A(3,3:6) = (m(3:6)).^3.*exp(m(3:6)*LL);

A(4,1:2) = (m(1:2)).^4.*exp(-m(1:2)*LL);
A(4,3:6) = (m(3:6)).^4;

A(5,1:2) = (m(1:2)).^5.*exp(-m(1:2)*LL);
A(5,3:6) = (m(3:6)).^5;

A(6,1:2) = exp(-m(1:2)*LL);
A(6,3:6) = 1;

A(7,1:2) = m(1:2).*exp(-m(1:2)*LL);
A(7,3:6) = m(3:6);

A(6,7) = -1;
A(7,7) = 1i*k*(H/Hs);

b = zeros(7,1);
b(6) = Ap;   b(7) = Ap*1i*k*(H/Hs);

sol = A\b;
c = sol(1:6);
R = sol(7);
end