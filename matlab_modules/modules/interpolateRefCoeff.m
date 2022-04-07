function V = interpolateRefCoeff(omega,omegaNew,Nev,filePath,TorC,NModes)
%% Interpolate the Diffraction Potential over the new-frequency space.
omegalist=omega(:);
omegaNewlist = omegaNew(:);

rd = zeros(1,NModes+1,length(omegalist)); %Diffraction RefCoeffs
rr = zeros(Nev,NModes+1,length(omegalist)); %Radiation RefCoeffs
for m=1:length(omegalist)
    rcDiff = load([filePath,'/RefCoeff_Dif/ref',TorC,num2str(m),'.dat']);
    rcRad = load([filePath,'/RefCoeff_Rad/ref',TorC,num2str(m),'.dat']);
    
    rd(1,:,m) = rcDiff(:,1)+1i*rcDiff(:,2);
    rr(:,:,m) = (rcRad(:,1:2:2*NModes+2)+1i*rcRad(:,2:2:2*NModes+2));
end

%% Saves more than one file now
for p=1:NModes+1
    % Interpolate the diffraction reflection Coefficients
    ReshapedRd = reshape(rd(1,p,:),[length(omega),length(omega)]);
    rdNew = interp2(real(omega),imag(omega),ReshapedRd,real(omegaNew),imag(omegaNew),'cubic');
    rdNew = rdNew(:);
    
    % Loop over each mode and interpolate the reflection coefficients
    rcNew = zeros(length(omegaNewlist),Nev);
    for m=1:Nev
        ReshapedRc = reshape(rr(m,p,:),[length(omega),length(omega)]);
        rcc = interp2(real(omega),imag(omega),ReshapedRc,real(omegaNew),imag(omegaNew),'cubic');
        rcNew(:,m) = rcc(:);
    end
    
    %% Write the interpolated reflection coefficients to a file.
    dlmwrite([filePath,'/Interpolated_R/ref',TorC,'DifRe_Mode',num2str(p-1),'.dat'],real(rdNew));
    dlmwrite([filePath,'/Interpolated_R/ref',TorC,'DifIm_Mode',num2str(p-1),'.dat'],imag(rdNew));
    dlmwrite([filePath,'/Interpolated_R/ref',TorC,'RadRe_Mode',num2str(p-1),'.dat'],real(rcNew));
    dlmwrite([filePath,'/Interpolated_R/ref',TorC,'RadIm_Mode',num2str(p-1),'.dat'],imag(rcNew));

end

V = 0;
end
