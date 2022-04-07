function mroots = dispersion_free_surface(alpha,N,h)
% mroots = dispersion_free_surface(alpha,N,h)
% calculates the positive imaginary and N first positive real solutions of -alpha = k*tan(k h)
% for complex alpha. It uses three methods - homotopy for
% starting with alpha =1, guess with linear expansion
% and a linear expansion.
% The first roots is positive imaginary and the next are the first N positive real
% ordered from smallest. 
%
% If the value for h is not given the default value is h = 1;

% it would be easy to write a much faster program for only real alpha

if nargin == 2
    h = 1;
else
    alpha = h*alpha; % scaling for the solution 
end

%first of all we calculate the root for alpha = 1;
mroots = zeros(1,N+1);

if N ==0; % the N = 0 case does not involve any of the special methods and is treated separately.
    count = 0;
    mroots(count+1) = homotopy(alpha,count);
else
    
count = 0;
mroots(count+1) = homotopy(alpha,count);
  count = count + 1;
while 0 <= 1
  %first of all we calculate the root for alpha = 1;
  
  mroots(count+1) = homotopy(alpha,count);
 if  abs(mroots(count + 1) - (i*count*pi + alpha/(i*count*pi))) < 0.01
   %sprintf('Now we can use the close guess')
     while 0 <=1
        mroots(count + 1) = oneroot(alpha,i*count*pi + alpha/(i*count*pi));
        %abs(mroots(count + 1) - (i*count*pi + alpha/(i*count*pi)))
        if abs(mroots(count + 1) - (i*count*pi + alpha/(i*count*pi))) < 1e-8
          % sprintf('Now we can work the rest out easily')
           mroots(count+1:N+1) = i*(count:N)*pi + alpha./(i*(count:N)*pi);
           count = N;
           break
        end
        if count ==N
           break
        end
        count = count + 1;
     end
  end
  
  if count == N
     break
  end
  count = count +1;

end

end

mroots = -i/h*mroots;mroots(1) = -mroots(1);

function mroot = homotopy(alpha,N)
%calculates the Nth root using the homotopy method

if N == 0;
   mroot = oneroot(1,1);
else
   mroot = oneroot(1,i*N*pi);
end

step =0.043;
if abs(alpha) < 1
    alphastep = ([1:-step:abs(alpha),abs(alpha)]);
else
    alphastep = ([1:step:abs(alpha),abs(alpha)]);
end

for k=2:length(alphastep)
        mroot = oneroot(alphastep(k),mroot);
end


if angle(alpha) > 0
    alphastep = abs(alpha)*exp(i*[0:pi/30:angle(alpha),angle(alpha)]);
else
    alphastep = abs(alpha)*exp(i*[0:-pi/30:angle(alpha),angle(alpha)]);
end


for k=2:length(alphastep)
   mroot = oneroot(alphastep(k),mroot);
end
  

function out = oneroot(alpha,guess)
%calculates the root nearest the root guess.
ans1 = guess+1;
out = guess;
while abs(ans1 - out) > 1e-9
    ans1 = out;
    out = ans1 - f(ans1,alpha)/difff(ans1);
end


function out = f(z,alpha)
out = z*tanh(z) - alpha;

function out = difff(z)
out = tanh(z) + z*sech(z).^2;
