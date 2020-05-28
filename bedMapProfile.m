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
    ppCav.pieces];

dlmwrite('./Meshes/BEDMAP2/iceDat.dat',arrIce,'precision',16,'delimiter','\t');

%% Run the meshing code and plot the results
fprintf('Running Mesher ....\n');
file='iceShelfBEDMAP2.edp';
ffpp=[ff,' -nw -ne ', file];
[aa,bb]=system(ffpp);
if(aa)
   error('Cannot run mesher. Check output\n'); 
end

%% Load the mesh 
[pts,seg,tri]=importfilemesh('Meshes/iceMeshBEDMAP.msh');
[pts1,seg1,tri1]=importfilemesh('Meshes/cavMeshBEDMAP.msh');

%% Plot the results
fig=figure(10);
subplot(2,2,[1,2]);
%set(fig,'Position',[284    23   952   634]);
set(fig,'Position',[0 0 2880 1800]);
ptsScaled=[pts(1,:)/1; pts(2,:)];
pts1Scaled=[pts1(1,:)/1; pts1(2,:)];
pp1=pdeplot(ptsScaled,seg,tri); hold on
pp2=pdeplot(pts1Scaled,seg1,tri1);
set(pp2,'Color','b')
title('Scaled $x-$axis $(0.1 \times)$')
axis equal
xlabel('$x$')
ylabel('$z$')
grid on

subplot(2,2,3)
pp1=pdeplot(pts,seg,tri); hold on
pp2=pdeplot(pts1,seg1,tri1);
set(pp2,'Color','b')
axis equal
xlim([0.15*10^5,0.25*10^5]);
title('Real $x-$axis $(1 \times)$')
grid on
xlabel('$x$')
ylabel('$z$')

subplot(2,2,4)
pp1=pdeplot(ptsScaled,seg,tri); hold on
pp2=pdeplot(pts1Scaled,seg1,tri1);
set(pp2,'Color','b')
axis equal
xlim([1.8*10^4,2*10^4]);
title('Scaled $x-$axis $(0.1 \times)$')
grid on
xlabel('$x$')
ylabel('$z$')
% %% Generate the cubic spline curve (This is to test. Real code in FF++)
% xq=linspace(min(iceCavInt(:,1)),max(iceCavInt(:,1)),1000);
% yq=zeros(length(xq),1);
% yq1=zeros(length(xq),1);
% for m=1:length(xq)
%     yq(m)=splineRecon(xq(m),ppIceCavInt.coefs,ppIceCavInt.breaks);
%     yq1(m)=splineRecon(xq(m),ppCav.coefs,ppCav.breaks);
% end