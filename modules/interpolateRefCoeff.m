function V = interpolateRefCoeff(omega,omegaNew,Nev,filePath)
%% Interpolate the Diffraction Potential over the new-frequency space.
omegalist=omega(:);
omegaNewlist = omegaNew(:);

rd = zeros(length(omegalist),1); %Diffraction RefCoeffs
rr = zeros(length(omegalist),Nev); %Radiation RefCoeffs
for m=1:length(omegalist)
    rcDiff = load([filePath,'/RefCoeff_Dif/refC',num2str(m),'.dat']);
    rcRad = load([filePath,'/RefCoeff_Rad/refC_',num2str(m),'.dat']);
    
    rd(m) = rcDiff(1)+1i*rcDiff(2);
    rr(m,:) = (rcRad(:,1)+1i*rcRad(:,2)).';
end

% Interpolate the diffraction reflection Coefficients
ReshapedRd = reshape(rd,[length(omega),length(omega)]);
rdNew = interp2(real(omega),imag(omega),ReshapedRd,real(omegaNew),imag(omegaNew),'cubic');
rdNew = rdNew(:);

% Loop over each mode and interpolate the reflection coefficients
rcNew = zeros(length(omegaNewlist),Nev);
for m=1:Nev
    ReshapedRc = reshape(rr(:,m),[length(omega),length(omega)]);
    rcc = interp2(real(omega),imag(omega),ReshapedRc,real(omegaNew),imag(omegaNew),'cubic');
    rcNew(:,m) = rcc(:);
end

%% Write the interpolated reflection coefficients to a file.
dlmwrite([filePath,'/Interpolated_R/refCDifRe.dat'],real(rdNew));
dlmwrite([filePath,'/Interpolated_R/refCDifIm.dat'],imag(rdNew));
dlmwrite([filePath,'/Interpolated_R/refCRadRe.dat'],real(rcNew));
dlmwrite([filePath,'/Interpolated_R/refCRadIm.dat'],imag(rcNew));


V = 0;
end
