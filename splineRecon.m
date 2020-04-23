function V=splineRecon(t,coeffs,breaks)
npieces=length(coeffs);
V=0;
for i=1:npieces
    a=coeffs(i,1);
    b=coeffs(i,2);
    c=coeffs(i,3);
    d=coeffs(i,4);
    
    if(t>=breaks(i) && t<breaks(i+1))
        V=a*(t-breaks(i))^3+b*(t-breaks(i))^2+c*(t-breaks(i))+d;
    end
end

end