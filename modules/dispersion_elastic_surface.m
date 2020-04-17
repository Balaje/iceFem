function mroots= dispersion_elastic_surface(alpha,beta,gamma,N,h)
%mroots = dispersion_ice_standard(alpha,beta,gamma,N)
%calculates the real and N imaginary (and complex) 
% solutions of alpha/(beta*z^4 + (1-alpha*gamma)) = - z*tan(z)
%. 
% this program is the same dispersion_ice except that the representation of the 
% roots is the standard one. The travelling root is negative imaginary and 
% the damped roots are positive real. 
%
% program using the programs written by tim

%rescaling 
alpha=alpha*h;
beta=beta/h^4;
gamma=gamma/h;

del = (1 - alpha*gamma)/alpha * (beta/alpha)^(-1/5);
H = (beta/alpha)^(-1/5);

K = RTS_ice_roots(del,H,N+1).';
mroots = -i*(beta/alpha)^(-1/5)*K(1:N+1);
mroots = mroots./h;

v=mroots(3);
mroots(3)=-mroots(1);
mroots(1)=v;


function K = RTS_ice_roots(del,H,N,varargin)
%% CALL: K = RTS_ice_roots(del,H,N); Z={del,H}
%% calculates the real root, 1 complex one, and N imaginary solutions of 
%% the dispersion relation 1/(K^4 + del) = K*tanh(K*H).
%%
%% DEPENDENCIES: RTS_imag_roots_ice.c (.mexglx), RTS_imag_root_ice.c (.mexglx).
do_test=0;

%% GET REAL ROOT:
Kgs=roots([1 0 0 0 del -1]);
guess0=(Kgs>0 & imag(Kgs)==0)'*Kgs;
K0=gen_root_ice(del,H,guess0); K=K0;

if N>0%%%%
   K=zeros(N+3,1); K(1)=K0; %miK=0*K;
   Del=pi/H; delc0=-Del^4; tol=1e-8;
   nc=(del/delc0)^.25;
   if del>=0
	nc=1;
   elseif nc~=round(nc)
	nc=round(nc+.5);%% round down because del<=nc^4*delc0
			%% iff (del/delc0)^.25>=nc.
   end%,nc

   if nc<=N%% GET "EASY" IMAGINARY ROOTS:
     w=RTS_imag_roots_ice(del,H,nc,N);%w/pi,length(w), find(w==0)
     j_notcvg=find( w(2:length(w))==0 );
     for r=1:length(j_notcvg)
	%j=j_notcvg(r)+1; i1=j*pi; i0=i1-pi/2;
	j=j_notcvg(r)+nc; i1=j*pi; i0=i1-pi/2;
	w(j+1-nc)=RTS_imag_root_ice(del,H,i0,i1);
     end%,w/pi
     if w(1)==0%%do nc-th interval specially
	i1=nc*pi; i0=i1-pi;
	w(1)=RTS_imag_root_ice(del,H,i0-tol,i1+tol);
     end
     jj=nc:N;
     K(jj+3)=i*w/H; %w/pi
   else
     i1=nc*pi; i0=i1-pi;
     K(nc+3)=i/H*RTS_imag_root_ice(del,H,i0-tol,i1+tol);
   end
   %% GET MORE DIFFICULT IMAGINARY ROOTS:
   for j=1:min(N,nc-1)
	i0=(j-1)*pi; i1=i0+pi/2; 
	w=RTS_imag_root_ice(del,H,i0,i1); K(j+3)=i*w/H;
   end

   %% ATTEMPT TO FIND COMPLEX ROOTS:
   if nargin==4
	guess1=varargin{:};
   else
	guess1=(Kgs>0 & imag(Kgs)>0)'*Kgs;
   end
   k=gen_root_ice(del,H,guess1); tol=1e-8;
   [k_cx,k_im]=is_complex(k,{K(nc+3),nc*pi,H});
   if k_cx==0 & k_im==0%% no complex or imaginary root found.
	K=complex_roots(del,H,K,nc);
   elseif k_cx==0%% found another imaginary root.
	K(3)=k_im; K=complex_roots(del,H,K,nc);%disp('2nd imag root')
   else
	K(2)=k_cx; K(3)=-K(2)';
   end
   K=K(1:N+3);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if do_test%% do some tests
	subplot(1,2,1), plot(real(K),imag(K)/Del,'x'), hold on;	
	plot(-real(K),-imag(K)/Del,'x');
	Xlim=(K(1)+1)*[-1 1]; Ylim=5*[-1 1];
	plot(Xlim,0*Xlim,':r'), plot(0*Xlim,Ylim,':r'), hold off;
	xlim(Xlim), ylim(Ylim);
	%%
	n=5000; wf=(0:n)'*5*Del/n;
	pf=(wf.^4+del).*wf.*sin(wf*H)+cos(wf*H);
	scale=abs(wf).^5+1;
	subplot(1,2,2), plot(wf/Del,[pf./scale,0*wf]); hold on;
	%%
	jj=find(abs(real(K)) < tol);
	plot(real(-i*K(jj))/Del,0*K(jj),'.r'); hold off;
	xlim([0 5]);
   end
end%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gen_root_ice: given an initial starting guess it uses Newton-Rhapson to find
%%	the nearest root of the dispersion relation
%%	f=(k^4+del)*k*sinh(k*H)-cosh(k*H)=0.
function [k,df]=gen_root_ice(del,H,guess)

tol=1e-8; max_reps=50; k0=guess;

[f,df]=f_df(k0,del,H); 
dk=f/df; k=k0-dk; reps=0;
while (abs(dk) > tol) & (reps < max_reps)
	reps=reps+1; k0=k;
	[f,df]=f_df(k0,del,H);
	dk=f/df; k=k0-dk;
end

if reps==max_reps
	k=0;
end

function [f,df]=f_df(K,del,H)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% f/df, where f has the same zeros as of the dispersion function, 
%% is the correction term in Newton's method for finding zeros in f

Lam=K^4+del; Lampr=5*K^4+del; x=7.5;
if real(K*H)<=x
	f=Lam*K*sinh(K*H)-cosh(K*H);
	df=Lam*K*H*cosh(K*H)+(Lampr-H)*sinh(K*H);
else
	f=Lam*K*tanh(K*H)-1;
	df=Lam*K*H+(Lampr-H)*tanh(K*H);
end
%% END OF MAIN PROGRAM---REQUIRED SUBROUTINES FOLLOW.%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% is_complex: checks whether the root k is complex or not
%   	---if it is y_cx=k, otherwise y_cx=0;
%%  	if a 2nd argument K, containing all the other roots is passed,
%%  	this function also compares k to all of them
%%  	---if is imaginary and is different to all the other roots,
%%  	y_im=k, otherwise y_im=0.
function [y_cx,y_im]=is_complex(k,varargin)

tol=1e-8; y_cx=0; y_im=0;
k=abs(real(k))+i*abs(imag(k));
if real(k)>tol
    if imag(k)>tol%%has converged to a complex root.
	y_cx=k;
    %%else has converged to a real root.
    end
elseif nargin==2%%has converged to an imaginary root;
    %%check if it is different from K(nc) (if it is provided).
    	k=i*sign(imag(k))*imag(k); 
	Z=varargin{1}; Kc=Z{1}; i1=Z{2}; H=Z{3}; i0=i1-pi; w=imag(k)*H;
    if abs(k-Kc)>tol/H & (w-i0+tol)*(i1+tol-w)>0
	y_im=k;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finds the complex roots/extra imaginary roots if they haven't been found yet.
function K=complex_roots(del,H,K,nc)

%Tst=p_fxn(-i*K*H,del,H)
tol=1e-8; Del=pi/H; N=length(K)-3; req23=1;
I1=nc*pi; I0=I1-pi;
k=K(nc+3); %k_=K([4:nc+2,nc+4:N+3]);
w=imag(k)*H; [p,dp]=p_fxn(w,del,H);

%% IN CERTAIN CASES, CAN IMMEDIATELY TELL THAT THERE ARE >1 IMAG ROOTS
%% IN INTERVAL #nc:
if K(3)~=0
	k_1=K(3); w_1=imag(k_1)*H; [p,dp_1]=p_fxn(w_1,del,H);
	%tst_roots=p_fxn([w w_1],del,H)
	req23=0;%%have 2 imag roots in i/val #nc => must be another.
    if w > w_1
	W = w; dP = dp; W_1 = w_1; dP_1 = dp_1;
    else
	W = w_1; dP = dp_1; W_1 = w; dP_1 = dp;
    end
    if dp*dp_1>0%%the 3rd root is between the other 2
	i0=W_1+tol; i1=W-tol;
	w_2=RTS_imag_root_ice(del,H,i0,i1); K(2)=i*w_2/H;
    elseif dP*(-1)^nc<0%%the 3rd root is between W & I1;
	i0=W+tol; i1=I1;
	w_2=RTS_imag_root_ice(del,H,i0,i1); K(2)=i*w_2/H;
    else%%the 3rd root is between I_0 & W_1;
	i0=I0; i1=W_1-tol; %p_fxn([i0 i1],del,H)
	w_2=RTS_imag_root_ice(del,H,i0,i1); K(2)=i*w_2/H;
	%K(2)=imag_root_iceB(del,H,[i0,i1])%%
    end
elseif dp*(-1)^nc<0%%only have 1 root in i/val,
		%% but there are definitely 3 roots in the interval.
	i0=I0; i1=w-tol; req23=0;
	w_=RTS_imag_root_ice(del,H,i0,i1); K(3)=i*w_/H;
	i0=w+tol; i1=I1;
	w_=RTS_imag_root_ice(del,H,i0,i1); K(2)=i*w_/H;
else%%check whether k is a triple root
	[H3,del3,w3]=triple_root(nc);
    if abs(H-H3)<tol & abs(del-del3)<tol
	K(2:3)=k; req23=0;
    end
end
%% NB neither N-R nor bisection will cvg to a dbl root;
%% but bisection will locate a triple root.

%% DETERMINE WHETHER THERE ARE OTHER IMAG ROOTS IN I/VAL #nc; IF SO, FIND THEM.
if req23==1
	%%first find where the double roots are.
	[w20,del20,w21,del21]=double_roots(H,nc,{H3,del3,w3});
    if abs(del-del20)<tol%% the 2 roots are on the 1st dbl root.
	K(2:3)=i*w20/H; req23=0;
    elseif abs(del-del21)<tol%% the 2 roots are on the 2nd dbl root.
	K(2:3)=i*w21/H; req23=0;
    elseif del>del20+tol & del<del21-tol
	%% the 2 roots are between the 2 dbl roots.
	nn=10; req23=0;
	Dw=(w21-w20)/nn; ww=w20:Dw:w21; pp=p_fxn(ww,del,H);
      if w>w21
	j=find(pp*dp >= 0);
	while isempty(j)==1
		Dw=Dw/10; ww=w20:Dw:w21; pp=p_fxn(ww,del,H); j=find(pp*dp >= 0);
	end
      else%%w<w20
	  j=find(pp*dp <= 0);
	  while isempty(j)==1
		Dw=Dw/10; ww=w20:Dw:w21; pp=p_fxn(ww,del,H); j=find(pp*dp <= 0);
	  end
      end
      if length(j)==1
	  if pp(j)==0
		K(2)=i*ww(j)/H; i0=ww(j)+tol; i1=ww(j+1);
		K(3)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
	  else
		i0=ww(j-1); i1=ww(j); K(2)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
		i0=i1; i1=ww(j+1); K(3)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
	  end
      else
	  nj=length(j);
	  if pp(j(1))==0 & pp(j(nj))==0
		K(2)=i*ww(j(1))/H; K(3)=i*ww(j(nj))/H;
	  elseif pp(j(1))==0
		K(2)=i*ww(j(1))/H; i0=ww(j(nj)); i1=ww(j(nj)+1);
		K(3)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
	  elseif pp(j(nj))==0
		K(3)=i*ww(j(nj))/H; i0=ww(j(1)-1); i1=ww(j(1));
		K(2)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
	  else
		i0=ww(j(1)-1); i1=ww(j(1));
		K(2)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
		i0=ww(j(nj)); i1=ww(j(nj)+1);
		K(3)=i*RTS_imag_root_ice(del,H,i0,i1)/H;
	  end
      end
    end
end

%% IF MORE IMAGINARY ROOTS HAVE BEEN FOUND WE CAN FINISH
%% BY SORTING THE IMAGINARY ROOTS IN ORDER OF INCREASING MAGNITUDE.
if req23==0
	j=2:nc+3; w=sort(imag(K(j))); K(j)=i*w;
end

%% OTHERWISE THERE ARE DEFINITELY 2 AS YET UNFOUND COMPLEX ROOTS.

%% 1ST: TRY USING THE SHALLOW WATER ROOTS AS A GUESS
if req23==1 & nc==1 & del<del20 & del>-1/3/(2*H)^(2/3);
	%%try using shallow water roots as a guess
	w=roots([1 0 del*H^(2/3) -1]); k0=sqrt(w*H^(1/3));
	k0=k0(find(imag(k0)>0));
	k=gen_root_ice(del,H,k0); k=is_complex(k);
    if k~=0
	K(2)=k; K(3)=-k'; req23=0;
    end	
end

%% 2ND: IF THAT DOESN'T WORK, WE KNOW FOR A GIVEN VALUE OF del,
%% IF ONE WATER DEPTH HAS COMPLEX ROOTS, THERE WILL ALSO BE COMPLEX ROOTS
%% FOR ALL HIGHER WATER DEPTHS. HENCE WE KEEP INCREASING THE DEPTH TILL
%% THEY ARE FOUND SUCCESFULLY; WHEN WE DO FIND A COMPLEX ROOT,
%% WE USE THAT RESULT TO APPROXIMATE THE RESULT FOR A SLIGHTLY SMALLER DEPTH,
%% AND KEEP DECREASING IT UNTIL WE GET BACK TO THE ORIGINAL DEPTH.
if req23==1
	%disp('looking with "complex_seed.m"')
	H_=H+.1;
	r=roots([1 0 0 0 del -1]); gs=find(r>0 & imag(r)>0);
	[k_,df]=gen_root_ice(del,H_,gs); k_=is_complex(k_);
    while k_==0 %& H_<2.5
	H_=H_+.1;
	r=roots([1 0 0 0 del -1]); gs=r(find(r>0 & imag(r)>0));
	[k_,df]=gen_root_ice(del,H_,gs); k_=is_complex(k_);
    end
	Z={k_,df,H_}; %dH=-min(.05,H);
	k=complex_seed(del,H,Z);
	K(2)=k; K(3)=-k';
end

function [p,dp]=p_fxn(w,del,H)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=(w.^4+del*H^4).*w.*sin(w)+H^5*cos(w);
dp=(5*w.^4+H^4*(del-H)).*sin(w)+(w.^4+del*H^4).*w.*cos(w);

function y=complex_seed(del,H,Z)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=Z{1}; df_k=Z{2}; H1=Z{3}; dH=-.05; tol=1e-8;
if length(Z)==4
	dH=Z{4};
end
while H1>H & k~=0 & abs(H1-H)>=tol
%% sts get rounding errors which make it seem like H>H1,
%% when they are actually equal
	H1=H1+dH;
	Lam=k^4+del; df_H=Lam*k^2*cosh(k*H)-k*sinh(k*H);
	gs=k-df_H*dH/df_k; [k,df_k]=gen_root_ice(del,H,gs);
	k=is_complex(k);
    if k~=0
	Z={k,df_k,H1};
    end
end
if abs(H1-H)<tol & k~=0
	y=k;
elseif k==0;%%find last successful result, decrease the step-size & start over.
	Z{4}=dH/2; y=complex_seed(del,H,Z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% double_roots: finds the double roots in an interval, if there are any.
%% 	if p=(w^4+del*H^4)*w*sin(w)+H^5*cos(w), the roots are found by setting
%%	p=p'=0 and eliminating del*H^4 to give q(w,H)=0 (cf. q_fxn below).
%%	q=0 is solved by using matlab's fzero routine.
function [w20,del20,w21,del21]=double_roots(H,nc,varargin)

if nargin==3
	Z=varargin{:}; H3=Z{1}; del3=Z{2}; w3=Z{3}; %[H3,H]
else
	[H3,del3,w3]=triple_root(nc);
end

I1=nc*pi; I0=I1-pi; tol=1e-8; I0=max([I0 tol]);
jc=2*nc-1; Hc=( (jc*pi)^4/4 )^.2; wc=jc*pi/2;
if H>H3%%no roots.
	w20=Inf; w21=Inf;
elseif H==H3%%get a triple root.
	w20=w3; w21=w3; 
else
	%%know w21 is always in (w3,I1)
	i0=w3; i1=I1; w21=double_root(H,i0,i1);
    if H==Hc %%know w20=wc.
	w20=wc;
    elseif H<Hc %%know w20 is to the left of wc
	i0=I0; i1=wc; w20=double_root(H,i0,i1);
    else %%w20 is in (I0,w3)
	i0=wc; i1=w3; w20=double_root(H,i0,i1);
    end
end
del20=-((w20/H)^4+H*cot(w20)/w20);
del21=-((w21/H)^4+H*cot(w21)/w21);

function w=double_root(H,i0,i1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=1e-8; w=fzero(@q_fxn,[i0, i1],optimset('TolX',tol),H);

function q=q_fxn(w,H)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c2=cos(2*w); s2=sin(2*w);
q=0*w-2*H^5; j=find(w~=0);
q(j)=2*w(j).^4.*(1-c2(j))-H^5*(1+s2(j)/2./w(j));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% triple_roots: returns the triple roots of the dispersion relation in the 
%% 	nth interval. If q(w)=2*w^4*(1-c2)-H^5*(1+s2/2/w), 
%% 	(c2=cos(2*w), s2=sin(2*w)),
%%	the roots are found by setting q=q'=0 and eliminating H^5.
function [H,del,w]=triple_root(n)

I=[(2*n-1)*pi/2 n*pi]; w=fzero(@r_fxn,I,optimset('TolX',1e-8));
num=2-2*cos(2*w)+w*sin(2*w);
den=2*w*cos(2*w)-sin(2*w); H5=8*w^5*num/den;
H=H5^.2; del=-( (w/H)^4+H*cot(w)/w );

function y=r_fxn(w)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%function to be "zeroed"
c2=cos(2*w); s2=sin(2*w);
lhs=(4*w+2*s2).*(2-2*c2+w.*s2);
rhs=(1-c2).*(2*w.*c2-s2); y=rhs-lhs;

function x=RTS_imag_roots_ice(del,H,N0,N1)

EPS=1e-9;
MAXIT=50;
x=zeros(N1+1-N0,1);

for j=N0:N1
  %% initial guess
  w=j*pi;
  i1=w;
  i0=i1-pi/2;
  if j==N0
    i0=i1-pi;
  end
  for its=1:MAXIT
    w4=w*w*w*w;
    H4=H*H*H*H;
    p=(w4+del*H4)*w*sin(w)+H4*H*cos(w);
    dp=( 5*w4+H4*(del-H) )*sin(w)+(w4+del*H4)*w*cos(w);
    w1=w;
    w=w1-p/dp;
    if abs(w-w1) <= EPS
      break;
    end
  end
  %%check root has converged & is in correct interval:
  if ( (its < MAXIT) & ((w-i0+EPS)*(i1+EPS-w)>=0) )
    x(j+1-N0)=w;
  else
    x(j+1-N0)=0;
  end
end

function [x,exitflag]=RTS_imag_root_ice(del,H,i0,i1)

%% convergence parameters:
EPS=1e-9;
MAXIT=100;

H4=H*H*H*H;
x4=i0*i0*i0*i0;
Lam=x4+del*H4;
fl=Lam*i0*sin(i0)+H*H4*cos(i0);
%%
x4=i1*i1*i1*i1;
Lam=x4+del*H4;
fh=Lam*i1*sin(i1)+H4*H*cos(i1);

%% check there is a root in interval
if (fl > 0 & fh > 0) | (fl < 0 & fh < 0 )
  x=0;
  exitflag=0;
  return;
end

%% if either endpoint is "close enough", just use that for root
if (fl == 0)
  x = i0;
  exitflag=1;
  return;
end
if (fh == 0)
  x = i1;
  exitflag=2;
  return;
end

%% define search direction to be in dirn of increasing p
if (fl < 0.0)
  xl=i0;
  xh=i1;
else
  xl=i1;
  xh=i0;
end

%% initial guess for NR:
x=0.5*(i0+i1);
dxold=i1-i0;
dx=dxold;
x4=x*x*x*x;
Lam=x4+del*H4;
p=Lam*x*sin(x)+H4*H*cos(x);
dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);

%% find a root using N-R:*/
for its=1:MAXIT
  crit1=( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0 );
  crit2=( abs(2*p) > abs(dxold*dp) );
  if crit1 | crit2
    %% bisect if NR jumps out of range,
    %% or if not decreasing fast enough.*/
    dxold=dx;
    dx=0.5*(xh-xl);
    x=xl+dx;
    if (x == xl)
      %% change in root is negligible.
      exitflag=3;
      return;
    end
  else
    %% use NR
    dxold=dx;
    dx=p/dp;
    x0=x;
    x=x0-dx;
    if x == x0
      exitflag=3;
      return;
      %% change in root is negligible.
    end
  end
  if (abs(dx) < EPS)
    exitflag=3;
    return;
  end
  x4=x*x*x*x;
  Lam=x4+del*H4;
  p=Lam*x*sin(x)+H4*H*cos(x);
  dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);
  if (p < 0.0 ) %% maintain the bracket on the root
    xl=x;
  else
    xh=x;
  end
end
x=0;
exitflag=4;
return;
%% if program reaches here it
%% has done MAXIT runs but hasn't converged yet