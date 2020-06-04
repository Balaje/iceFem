%% This is a program to generate the cavity and shelf profile using the BEDMAP2 dataset. 
% 1) The cavity is assumed to be restricted to the length of the ice.

clc
clear
close all
global ff

%% Generate the map of Antarctica and the ice--shelves
antmap
load coast
patchm(lat,long, [0.5,0.5,0.5]);
bedmap2 patchshelves
[hice,hbed,hwater]=bedmap2_profile();

%% Extract the profile.
iceCoord=hice.Vertices;
bedCoord=hbed.Vertices;

[rowIce,~]=find(iceCoord(:,2)~=0);
[freeSurf,~]=find(iceCoord(:,2)==0);
icePts=[iceCoord(rowIce,1), iceCoord(rowIce,2)];
bedCoord1 = bedCoord(find(bedCoord(1:end-2,1)>=min(icePts(1:end-2,1))),:);
bedCoord1 = bedCoord1(find(bedCoord1(1:end-2,1)<=max(icePts(1:end-2,1))),:);
[rowCav,~]=find(bedCoord1(:,2)~=0);


%% Perform a coordinate transformation and write the data onto a file.
ice=icePts;
if(iceCoord(end,1)==iceCoord(end,2))
    frontloc=iceCoord(freeSurf(end/2),1);
    a=1;
else
    frontloc=iceCoord(freeSurf(end),1);
    a=2;
end
xf=frontloc;
% Extract the points after transformation
ice1=[(3-2*a)*(ice(:,1)-xf)*1000, ice(:,2)]; % Change to km
ice1=[ice1(:,1)-min(ice1(:,1)),ice1(:,2)];
bedPts=[(3-2*a)*(bedCoord1(rowCav,1)-xf)*1000, bedCoord1(rowCav,2)]; % Change to km
bedPts=[bedPts(:,1)-min(bedPts(:,1)), bedPts(:,2)];
iceCavInt=ice1(find(ice1(:,2)<0),:);
iceTop=ice1(find(ice1(:,2)>0),:);

%% Plot the data to verify
plot(iceCavInt(:,1)/1000,iceCavInt(:,2),'m.','LineWidth',2);
hold on
plot(iceTop(:,1)/1000,iceTop(:,2),'m.','LineWidth',2);
plot(bedPts(:,1)/1000,bedPts(:,2),'m.','LineWidth',2);
pause(1);

%%
ppIceCavInt=spline(iceCavInt(:,1),iceCavInt(:,2));
ppIceTop=spline(iceTop(:,1),iceTop(:,2));
ppCav=spline(bedPts(:,1),bedPts(:,2));

% Write the cubic spline descriptions
dlmwrite('./Meshes/BEDMAP2/iceCavInt_Coeffs.dat',ppIceCavInt.coefs,'precision',16,'delimiter','\t');
dlmwrite('./Meshes/BEDMAP2/iceCavInt_Breaks.dat',ppIceCavInt.breaks','precision',16,'delimiter','\t');
dlmwrite('./Meshes/BEDMAP2/iceTop_Coeffs.dat',ppIceTop.coefs,'precision',16,'delimiter','\t');
dlmwrite('./Meshes/BEDMAP2/iceTop_Breaks.dat',ppIceTop.breaks','precision',16,'delimiter','\t');
dlmwrite('./Meshes/BEDMAP2/cavBed_Coeffs.dat',ppCav.coefs,'precision',16,'delimiter','\t');
dlmwrite('./Meshes/BEDMAP2/cavBed_Breaks.dat',ppCav.breaks','precision',16,'delimiter','\t');

% Write the ice-data onto a file
arrIce=[min(iceCavInt(:,1)); % Bottom left Ice
    max(iceCavInt(:,1)); % Bottom right Ice
    ppIceCavInt.pieces; % Number of pieces
    min(bedPts(:,1)); % Bottom left Cavity
    max(bedPts(:,1)); % Bottom right Cavity
    ppCav.pieces;
    fnval(ppIceCavInt,max(bedPts(:,1)))];

dlmwrite('./Meshes/BEDMAP2/iceDat.dat',arrIce,'precision',16,'delimiter','\t');

%% Run the meshing code and plot the results
% fprintf('Running Mesher ....\n');
% file='iceShelfBEDMAP2.edp';
% ffpp=[ff,' -nw -ne ', file];
% [aa,bb]=system(ffpp);
% if(aa)
%    error('Cannot run mesher. Check output\n'); 
% end

Tr=200;
file='solveBEDMAP2.edp';
fprintf('Running Program ....\n');
ffpp=[ff,' -nw -ne ',file,' -isMesh ',num2str(1)...
    ,' -Tr ',num2str(Tr),' -Ti ',num2str(0),' -nborders ',num2str(4)];
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

%% Plot the results
UX=importfiledata('xDisp0.bb');
UY=importfiledata('yDisp0.bb');
IUX=importfiledata('xDispI0.bb');
IUY=importfiledata('yDispI0.bb');
POT=importfiledata('potentialCav0.bb');

figure(4);
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