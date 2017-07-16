function plot_spatialAMcorr(comp,iT,variable,do_covar)
% plot_spatialcorr_variables.m----------------
%
%  Contour plots of the covariance between model variables and global AM
%  26 April 2010
%  ModificationsL
%	6 July 2010 - make into a function for easier handling
%
%  INPUT:
%	comp - component (1 or 2 for eq or axial)
%	iT - timescale to show
%	variable - which variable to plot
%	do_covar - set to 1 to do covariance, 0 to do correlation (not implemented yet!!)


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

% other plot settings

if do_covar == 1
  TITLE = 'Covariance';
else
  TITLE = 'Correlation';
end

% cycle the components and timescales and gather correaltions
R = zeros(128,64);

[R(:,:),L,lon,lat] = spatialcorr_variables(comp,TT(iT,:),variable,do_covar);    % covariance


% define plotting properties
XT = [-180 -120 -60 0 60 120 180];
YT = [-90 -60 -30 0 30 60 90];

ifig = do_covar+1;

  figure(ifig),clf
    m_proj('miller');
    %m_contourf(lon,lat,squeeze(real(R(comp,iT,:,:)))')
    if do_covar == 1
      m_pcolor(lon,lat,R');
    else
      m_pcolor(lon,lat,abs(R)');
    end
    shading interp
    m_coast('line','Color',zeros(3,1),'LineWidth',2);
    if do_covar == 0; set(gca,'Clim',[0 1]); end
    title(TITLE)
    m_grid('xtick',XT,'tickdir','out','ytick',YT,'yaxislocation','right','box','fancy','fontsize',7);
    colorbar('Location','WestOutside')








