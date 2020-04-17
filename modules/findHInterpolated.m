function [Hmat,Fmat] = findHInterpolated(omegaSpace,H,F,omegaq)
Nev = size(F,2);
%% Loop over the modes of all frequencies and find the interpolated values.
HH = zeros(Nev^2,1);
FF= zeros(Nev^2,1);
for m=1:Nev^2
    Hm = reshape(H(:,m),[length(omegaSpace),length(omegaSpace)]);  
    HH(m) = interp2(real(omegaSpace),imag(omegaSpace),Hm,real(omegaq),imag(omegaq),'makima');
    
    if(m<=Nev)
        Fm = reshape(F(:,m), [length(omegaSpace),length(omegaSpace)]);
        FF(m) = interp2(real(omegaSpace),imag(omegaSpace),Fm,real(omegaq),imag(omegaq),'makima');        
    end
end

Hmat = reshape(HH,[Nev,Nev]);
Fmat = FF;


end