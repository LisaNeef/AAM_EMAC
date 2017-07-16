% plot_spatialAMcorr.m----------------
%
%  Contour plots of the covariance between spatial components of AM and the global value
%  22 April 2010
%  Mods:
%	11 May 2010 - make into a function for easier handling.
%
%  INPUTS:
%	comp = component (1 for equatorial, 2 for axial)
%	iT: select timescale 
%		1 = subseasonal
%		2 = terannual
%		3 = semiannual
%		4 = annual
%		5 = interanual
%	term = 'm' for mass or 'w' for wind
%-------------------------------------------------------------------------
function plot_spatialAMcorr(comp,iT,term)


%---paths and file settings
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/m_map/');

%---select the components and the timescales that we want to plot.

TT = zeros(5,2);
TT(1,:) = [30,90];      % subseasonal.
TT(2,:) = [90,150];     % terannual
TT(3,:) = [150,230];    % semiannual
TT(4,:) = [230,450];    % annual
TT(5,:) = [730,1825];   % interannual
nT = size(TT,1);


% define plotting properties
TITLE = {'Covariance','Correlation'};
if iT == 1, Tbig = {'\chi_{eq} Subseasonal','\Delta LOD Subseasonal'}; end
if iT == 2, Tbig = {'\chi_{eq} Terannual','\Delta LOD Terannual'}; end
if iT == 3, Tbig = {'\chi_{eq} Semiannual','\Delta LOD Semiannual'}; end
if iT == 4, Tbig = {'\chi_{eq} Annual','\Delta LOD Annual'}; end
if iT == 5, Tbig = {'\chi_{eq} Interannual','\Delta LOD Interannual'}; end
   

XT = [-180 -120 -60 0 60 120 180];
YT = [-90 -60 -30 0 30 60 90];


% cycle the components and timescales and gather correaltions
%nlat = 46;
nlat = 64;

R = zeros(2,128,nlat);
  [R(1,:,:),L,lon,lat] = spatialcorr_AM(term,comp,TT(iT,:),1);	% covariance
  [R(2,:,:),L,lon,lat] = spatialcorr_AM(term,comp,TT(iT,:),0);	% correlation



% --cycle through components and timescales and make tons of plots!

figure(1),clf
subplot1(1,2,'Gap',[.03 .01]);
for iterm = 1:2
    subplot1(iterm);
    m_proj('miller');
    m_pcolor(lon,lat,squeeze(abs(R))');
    %m_pcolor(lon,lat,squeeze(abs(R(iterm,:,:)))');
    shading interp
    m_coast('line','Color',zeros(3,1),'LineWidth',2);
    if iterm == 2, set(gca,'Clim',[0 1]), end
    %title(TITLE(iterm))
    m_grid('xtick',XT,'tickdir','out','ytick',YT,'yaxislocation','right','box','fancy','fontsize',7);
    colorbar('Location','SouthOutside')
end
%text(-1.7*pi,-1.0*pi,Tbig(comp))


