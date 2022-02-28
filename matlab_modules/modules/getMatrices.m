function HNew = getMatrices(omega0,omega,Nev,filePath)
H = zeros(Nev^2,length(omega(:)));

%% Get all the matrices for the frequencies
for iter=1:length(omega(:))
    HRe = load([filePath,'/ReH',num2str(iter),'.dat']);
    HIm = load([filePath,'/ImH',num2str(iter),'.dat']);
    HMat = HRe + 1i*HIm;        
    
    % Store all the Entries - (Frequency-wise)
    H(:,iter) = HMat(:);
end

%% Interpolate to get the new-matrix
HNew = zeros(Nev^2,1);

for iter=1:Nev^2
    ReH = reshape(H(iter,:), [length(omega),length(omega)]);   
    HH = interp2(real(omega),imag(omega), ReH, real(omega0),imag(omega0),'makima');    
    HNew(iter,1) = HH(:);
end

HNew = reshape(HNew,[Nev,Nev]);

end