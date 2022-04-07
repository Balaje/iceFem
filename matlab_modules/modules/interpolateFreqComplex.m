function V = interpolateFreqComplex(omega,omegaNew,Nev,filePath)
%% Interpolate the matrix to a larger frequency space.
H = zeros(Nev^2,length(omega(:)));
F = zeros(Nev,length(omega(:)));

for iter=1:length(omega(:))
    HRe = load([filePath,'/ReH',num2str(iter),'.dat']);
    HIm = load([filePath,'/ImH',num2str(iter),'.dat']);
    HMat = HRe + 1i*HIm;  
    
    FRe = load([filePath,'/ReF',num2str(iter),'.dat']);
    FIm = load([filePath,'/ImF',num2str(iter),'.dat']);
    FMat = FRe + 1i*FIm;
    
    % Store all the Entries - (Frequency-wise)
    H(:,iter) = HMat(:);
    F(:,iter) = FMat;
end

%% Use Matlab's interp2 function to interpolate the values in the matrix H

HNew = zeros(Nev^2,length(omegaNew(:)));
FNew = zeros(Nev,length(omegaNew(:)));

for m=1:Nev^2
    ReH = reshape(H(m,:), [size(omega,1),size(omega,2)]);    
    HH = interp2(real(omega),imag(omega), ReH, real(omegaNew),imag(omegaNew),'cubic');    
    HNew(m,:) = HH(:);
    
    if(m <= Nev)
        ReF = reshape(F(m,:), [size(omega,1),size(omega,2)]);        
        FF = interp2(real(omega),imag(omega), ReF, real(omegaNew),imag(omegaNew),'cubic');
        FNew(m,:) = FF(:);
    end
    
%         figure(1);
%         subplot(1,2,1);
%         surf(real(omega),imag(omega),reshape(abs(H(m,:)),[length(omega),length(omega)]));
%         subplot(1,2,2);
%         surf(real(omegaNew),imag(omegaNew),reshape(abs(HNew(m,:)),[length(omegaNew),length(omegaNew)]));
%     
%         pause();
end

%% Write the interpolated H-Matrix and F-Vector to a file.
dlmwrite([filePath,'/Interpolated_H/ReH.dat'],real(HNew)','delimiter',' ');
dlmwrite([filePath,'/Interpolated_H/ImH.dat'],imag(HNew)','delimiter',' ');
dlmwrite([filePath,'/Interpolated_F/ReF.dat'],real(FNew)','delimiter',' ');
dlmwrite([filePath,'/Interpolated_F/ImF.dat'],imag(FNew)','delimiter',' ');

%% Solve for lambdaj for each frequency and dump the solution onto a file.
lambdaj = zeros(Nev,length(omegaNew(:)));
for p=1:length(omegaNew(:))
    Hmat = reshape(HNew(:,p),[Nev,Nev]);
    Fmat = FNew(:,p);
    lambdajNew = Hmat\Fmat;
    lambdaj(:,p) = lambdajNew;
end

dlmwrite([filePath,'/Interpolated_L/lambdaRe.dat'],real(lambdaj)','delimiter',' ');
dlmwrite([filePath,'/Interpolated_L/lambdaIm.dat'],imag(lambdaj)','delimiter',' ');

V = 0;
end
