clc
clear
close all

addpath('./modules/');

set(0,'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0,'defaultaxesfontsize',20);

envvar = [pwd,'/include'];
setenv('FF_INCLUDEPATH',envvar);

%% Generate the map of Antarctica and the ice--shelves
antmap
load coast
patchm(lat,long, [0.5,0.5,0.5]);
bedmap2 patchshelves
[hice,hbed,hwater]=bedmap2_profile();

%%
iceCoord=hice.Vertices;
bedCoord=hbed.Vertices;

[rowIce,~]=find(iceCoord(:,2)~=0);
icePts=[iceCoord(rowIce,1), iceCoord(rowIce,2)];
bedCoord1 = bedCoord(find(bedCoord(1:end-2,1)>=min(icePts(1:end-2,1))),:);
bedCoord1 = bedCoord1(find(bedCoord1(1:end-2,1)<=max(icePts(1:end-2,1))),:);
[rowCav,~]=find(bedCoord1(:,2)~=0);
bedPts=[bedCoord1(rowCav,1)*1000, bedCoord1(rowCav,2)]; % Change to km

%% Write the data onto a file.
ice=icePts;
ice=[ice(:,1)*1000, ice(:,2)]; % Change to km
iceCavInt=ice(find(ice(:,2)<0),:);
iceTop=ice(find(ice(:,2)>0),:);

%% Plot the data to verify
plot(iceCavInt(:,1)/1000,iceCavInt(:,2),'m.','LineWidth',2);
hold on
plot(iceTop(:,1)/1000,iceTop(:,2),'m.','LineWidth',2);
plot(bedPts(:,1)/1000,bedPts(:,2),'m.','LineWidth',2);
pause();

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
    ppCav.pieces];

dlmwrite('./Meshes/BEDMAP2/iceDat.dat',arrIce,'precision',16,'delimiter','\t');

%% Run the meshing code and plot the results
fprintf('Running Mesher ....\n');
[aa,bb]=system('/usr/local/ff++/openmpi-2.1/3.61-1/bin/FreeFem++ -ne -nw iceShelfBEDMAP2.edp');
if(aa)
   error('Cannot run mesher. Check output\n'); 
end

%% Load the mesh 
[pts,seg,tri]=importfilemesh('iceMeshBEDMAP.msh');
[pts1,seg1,tri1]=importfilemesh('cavMeshBEDMAP.msh');

%% Plot the results
fig=figure(10);
set(fig,'Position',[284    23   952   634]);
pp1=pdeplot(pts,seg,tri); hold on
pp2=pdeplot(pts1,seg1,tri1);
set(pp2,'Color','b')
%axis equal
axis square tight %To view the profiles
grid on

fig=figure(11);
set(fig,'Position',[284    23   952   634]);
pp1=pdeplot(pts,seg,tri); hold on
pp2=pdeplot(pts1,seg1,tri1);
set(pp2,'Color','b')
xlim([0.45*10^5,0.55*10^5]);
grid on

% %% Generate the cubic spline curve (This is to test. Real code in FF++)
% xq=linspace(min(iceCavInt(:,1)),max(iceCavInt(:,1)),1000);
% yq=zeros(length(xq),1);
% yq1=zeros(length(xq),1);
% for m=1:length(xq)
%     yq(m)=splineRecon(xq(m),ppIceCavInt.coefs,ppIceCavInt.breaks);
%     yq1(m)=splineRecon(xq(m),ppCav.coefs,ppCav.breaks);
% end