%% Program to generate the water and ice-shelf meshes.

clc
clear
close all
global ff

file='solveBEDMAP2.edp';

isMesh = true; % To run the mesher.

%% Load the map of the antarctica.
if(isMesh)
    antmap
    load coast
    patchm(lat,long, [0.5,0.5,0.5]);
    bedmap2 patchshelves
    [hice,hbed,hwater]=bedmap2_profile();
    
    %% Extract the points and construct the spline
    iceCoord=hice.Vertices;
    bedCoord=hbed.Vertices;
    iceCoord=[iceCoord(:,1)*1000, iceCoord(:,2)];
    bedCoord=[bedCoord(:,1)*1000, bedCoord(:,2)];
    
    % Find the top surface of the ice.
    idx1=find(gradient(iceCoord(:,1))>0);
    iceTop=iceCoord(idx1,:);
    iceTop=iceTop(find(abs(iceTop(:,2))>0),:);
    % Find the bottom surface of the ice.
    idx2=find(gradient(iceCoord(:,1))<0);
    iceCavInt=iceCoord(idx2,:);
    iceCavInt=iceCavInt(find(abs(iceCavInt(:,2))>0),:);    
    
    
    %% Extract the ocean-bed and find the grounding line.
    ppIceCavInt=spline(iceCavInt(:,1),iceCavInt(:,2));
    ppBedFull=spline(bedCoord(1:end-2,1),bedCoord(1:end-2,2));
    x0=min(iceCoord(:,1)); xf=max(iceCoord(:,1));
    x=linspace(x0,xf,2000);
    % Compute the cavity interfaces.
    fBed=fnval(ppBedFull,x);
    fCav=fnval(ppIceCavInt,x);
    idx=find(abs(fBed-fCav)>50);
    bedPts=[x(idx)', fBed(idx)'];
    
    %% Perform coordinate transformation (if applicable)
    if(~ismember([1,1],ismember(iceCoord,[0,0]),'rows'))
        xf=max(abs(iceCoord(:,1)));
        iceCavInt(:,1)=xf-iceCavInt(:,1);
        iceTop(:,1)=xf-iceTop(:,1);
        bedPts(:,1)=xf-bedPts(:,1);
    end
    
    subplot(2,1,1);
    plot(iceTop(:,1),iceTop(:,2),'r.-');
    hold on
    plot(iceCavInt(:,1),iceCavInt(:,2),'g.-');
    plot(bedPts(:,1),bedPts(:,2),'m.-');
    
    subplot(2,1,2);
    plot(iceTop(:,1),iceTop(:,2),'r.-');
    hold on
    plot(iceCavInt(:,1),iceCavInt(:,2),'g.-');
    plot(bedPts(:,1),bedPts(:,2),'m.-');
    
    pause();
    
    %% Construct the spline and write the data to a file.
    ppIceCavInt=spline(iceCavInt(:,1),iceCavInt(:,2));
    ppCav=spline(bedPts(:,1),bedPts(:,2));
    ppIceTop=spline(iceTop(:,1),iceTop(:,2));
    
    % Write the coefficients to a file.
    dlmwrite('./Meshes/BEDMAP2/iceCavInt_Coeffs.dat',ppIceCavInt.coefs,'precision',16,'delimiter','\t');
    dlmwrite('./Meshes/BEDMAP2/iceCavInt_Breaks.dat',ppIceCavInt.breaks','precision',16,'delimiter','\t');
    dlmwrite('./Meshes/BEDMAP2/iceTop_Coeffs.dat',ppIceTop.coefs,'precision',16,'delimiter','\t');
    dlmwrite('./Meshes/BEDMAP2/iceTop_Breaks.dat',ppIceTop.breaks','precision',16,'delimiter','\t');
    dlmwrite('./Meshes/BEDMAP2/cavBed_Coeffs.dat',ppCav.coefs,'precision',16,'delimiter','\t');
    dlmwrite('./Meshes/BEDMAP2/cavBed_Breaks.dat',ppCav.breaks','precision',16,'delimiter','\t');
    
    arrIce=[min(iceCavInt(:,1)); % Bottom left Ice
        max(iceCavInt(:,1)); % Bottom right Ice
        ppIceCavInt.pieces; % Number of pieces
        min(bedPts(:,1)); % Bottom left Cavity
        max(bedPts(:,1)); % Bottom right Cavity
        ppCav.pieces;
        fnval(ppIceCavInt,max(bedPts(:,1)))];
    
    dlmwrite('./Meshes/BEDMAP2/iceDat.dat',arrIce,'precision',16,'delimiter','\t');
    
end

%% Run the FreeFem++ code
Tr=200;
fprintf('Running Program ....\n');
ffpp=[ff,' -nw -ne ',file,' -isMesh ',num2str(isMesh)...
    ,' -Tr ',num2str(Tr),' -Ti ',num2str(0),' -isOpen ',num2str(1)];
[aa,bb]=system(ffpp);
if(aa)
    error('Cannot run mesher. Check output\n');
end
fprintf('FreeFem++ message: \n');
disp(bb);

%% Load the mesh
fprintf('FEM Problem solved. Importing Meshes ....\n');
[pts,seg,tri]=importfilemesh('Meshes/iceMeshBEDMAP.msh');
[pts1,seg1,tri1]=importfilemesh('Meshes/cavMeshBEDMAP.msh');

%% Visualize the solution.
UX=importfiledata('xDisp0.bb');
UY=importfiledata('yDisp0.bb');
IUX=importfiledata('xDispI0.bb');
IUY=importfiledata('yDispI0.bb');
POT=importfiledata('potentialCav0.bb');

figure(1);
subplot(2,1,1);
DISP=[pts(1,:)/20+(UX+1i*IUX); pts(2,:)+(UY+1i*IUY)];
pdeplot(real(DISP),seg,tri);
grid on
xlabel('$0.05\times x/L_c$')
ylabel('$z/L_c$')
title(['Incident Wave Period $T=',num2str(Tr),'$\,s'])
subplot(2,1,2);
pdeplot([pts1(1,:)/20; pts1(2,:)],seg1,tri1,'XYData',POT,'Colormap','jet');
xlabel('$0.05\times x/L_c$')
ylabel('$z/L_c$')
hold on
pdeplot([pts(1,:)/20; pts(2,:)],seg,tri);
xlabel('$0.05\times x/L_c$')
grid on