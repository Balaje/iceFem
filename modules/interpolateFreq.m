function [omegaNew,detH,condH] = interpolateFreq(a,b,omega,Nev,filePath,npts,isSolve)
%% Interpolate the matrix to a larger frequency space.
omegaNew= linspace(a,b,npts+1);

H = zeros(Nev^2,length(omega(:)));
F = zeros(Nev,length(omega(:)));

for iter=1:length(omega(:))
    HRe = load([filePath,'ReH',num2str(iter),'.dat']);
    HIm = load([filePath,'ImH',num2str(iter),'.dat']);
    HMat = HRe + 1i*HIm;
    % Store all the Entries - (Frequency-wise)
    H(:,iter) = HMat(:);
    
    if(isSolve)
        FRe = load([filePath,'ReF',num2str(iter),'.dat']);
        FIm = load([filePath,'ImF',num2str(iter),'.dat']);
        FMat = FRe + 1i*FIm;
        F(:,iter) = FMat;
    end
    
end

%% Use Matlab's interp2 function to interpolate the values in the matrix H

HNew = zeros(Nev^2,length(omegaNew(:)));
FNew = zeros(Nev,length(omegaNew(:)));

for m=1:Nev^2
    ReH = real(H(m,:));
    ImH = imag(H(m,:));
    
    ReHnew = interp1(omega, ReH, omegaNew, 'pchip');
    ImHnew = interp1(omega, ImH, omegaNew, 'pchip');
    HNew(m,:) = ReHnew(:) + 1i*ImHnew(:);
    
    if(m <= Nev && isSolve)
        ReF = real(F(m,:));
        ImF = imag(F(m,:));
        ReFnew = interp1(omega, ReF, omegaNew, 'pchip');
        ImFnew = interp1(omega, ImF, omegaNew, 'pchip');
        FNew(m,:) = ReFnew(:) + 1i*ImFnew(:);
    end
    
%     figure(1);
%     subplot(1,2,1);
%     plot(omega,abs(F(m,:)),'+-');
%     hold off
%     subplot(1,2,2);
%     plot(omegaNew,abs(FNew(m,:)),'+-');    
%     hold on
%     pause();
%     hold off
end

detH=zeros(length(omegaNew),1);
condH=zeros(length(omegaNew),1);
for m=1:length(omegaNew)
    AA=reshape(HNew(:,m),Nev,Nev);    
    detH(m)=det(AA);
    condH(m)=cond(AA);
end

%% Dump all the Interpolated matrices into a folder.
% [[Careful with the directory path]]
dlmwrite([filePath,'Interpolated_H/ReH.dat']...
    ,real(HNew),'delimiter',' ');
dlmwrite([filePath,'Interpolated_H/ImH.dat']...
    ,imag(HNew),'delimiter',' ');
dlmwrite([filePath,'Interpolated_F/ReF.dat']...
    ,real(FNew),'delimiter',' ');
dlmwrite([filePath,'Interpolated_F/ImF.dat']...
    ,imag(FNew),'delimiter',' ');

%% Solve for lambdaj for each frequency and dump the solution onto a file.
if(isSolve)
    lambdaj = zeros(Nev,length(omegaNew(:)));
    for p=1:length(omegaNew(:))
        Hmat = reshape(HNew(:,p),[Nev,Nev]);        
        Fmat = FNew(:,p);        
%         spparms('bandden',1);
%         spparms('piv_tol',1e-6);
%         spparms('sym_tol',1e-6);
        lambdajNew = Hmat.'\Fmat;        
        lambdaj(:,p) = lambdajNew;
    end
    
    dlmwrite([filePath,'Interpolated_L/lambdaRe.dat'],real(lambdaj)','delimiter',' ','Precision',16);
    dlmwrite([filePath,'Interpolated_L/lambdaIm.dat'],imag(lambdaj)','delimiter',' ','Precision',16);
end
end
