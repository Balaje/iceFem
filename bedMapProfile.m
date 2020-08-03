%% This is a program to generate the cavity and shelf profile using the BEDMAP2 dataset.
% 1) The cavity is assumed to be restricted to the length of the ice.

clc
clear
close all


%% Zoom in on to the Sulzberger Ice Shelf and generate the profiles.
mapzoom('Vollmer Island','mapwidth',100,'inset','se')
bedmap2 gl
bedmap2('patchcoast');
bedmap2('patchshelves','frame','on');
scalebar
scarclick
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
figure(3);
subplot(2,1,1);
plot(iceTop(:,1),iceTop(:,2),'r.-');
hold on
plot(iceCavInt(:,1),iceCavInt(:,2),'g.-');
plot(bedPts(:,1),bedPts(:,2),'m.-');

if(~ismember([1,1],ismember(iceCoord,[0,0]),'rows'))
    xf=max(abs(iceCavInt(:,1)));
    iceCavInt(:,1)=xf-iceCavInt(:,1);
    iceTop(:,1)=xf-iceTop(:,1);
    bedPts(:,1)=xf-bedPts(:,1);
    
    % Find only the cavity part underneath the ice-shelf. i.e. x>0
    bedPts=bedPts(find(bedPts(:,1)>0),:);
    ppBed=spline(bedPts(:,1),bedPts(:,2));
    bedPts=[bedPts(:,1), bedPts(:,2);
        0, fnval(ppBed,0);];
else
    xf=min(abs(iceCavInt(:,1)));
    iceCavInt(:,1)=iceCavInt(:,1)-xf;
    iceTop(:,1)=iceTop(:,1)-xf;
    bedPts(:,1)=bedPts(:,1)-xf;
    
    % Find only the cavity part underneath the ice-shelf. i.e. x>0
    bedPts=bedPts(find(bedPts(:,1)>0),:);
    ppBed=spline(bedPts(:,1),bedPts(:,2));
    bedPts=[0, fnval(ppBed,0);
        bedPts(:,1), bedPts(:,2)];
end

% Sort by ascending x
%     iceCavInt=sort(iceCavInt,1);
%     bedPts=sort(bedPts,1);
%     iceTop=sort(iceTop,1);
%

subplot(2,1,2);
plot(iceTop(:,1),iceTop(:,2),'r.-');
hold on
plot(iceCavInt(:,1),iceCavInt(:,2),'g.-');
plot(bedPts(:,1),bedPts(:,2),'m.-');

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