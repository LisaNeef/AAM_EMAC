 function plot_transfer_functions(comp,v)
%
% Make a contour plot of the transfer functions for the angular momentum integrals.
% Lisa Neef, started 16 May 2011
% 
% Inputs:
%	comp: vector component of AAM (1, 2, or 3)
%	v: which variable (u, v, m)
%
% Mods:
% 12 Jul 2011: cosmetic changes, so that you don't go blind looking at
% these.
%---------------------------------------------------


% generate a lat / lon grid
lat = -90:90;
lon = 0:360;

d2r = pi/180;
rlat = lat*d2r;
rlon = lon*d2r;

[LAT,LON] = meshgrid(rlat,rlon);

% compute transfer functions according to variable and component.

switch comp
  case 1
    switch v 
      case 'u'
        TF = sin(LAT).*cos(LAT).*cos(LON);
      case 'v'
        TF = -cos(LAT).*sin(LON);
      case 'm'
        TF = sin(LAT).*cos(LAT).*cos(LAT).*cos(LON);
    end 
  case 2
    switch v 
      case 'u'
        TF = sin(LAT).*cos(LAT).*sin(LON);
      case 'v'
        TF = cos(LAT).*cos(LON);
      case 'm'
        TF = sin(LAT).*cos(LAT).*cos(LAT).*sin(LON);
    end 
  case 3
    switch v 
      case 'u'
        TF = cos(LAT).*cos(LAT);
      case 'v'
        disp('Variable v does not affect the axial excitation function at all!')
        return
      case 'm'
        TF = cos(LAT).*cos(LAT).*cos(LAT);
    end 
end

% other plot settings 

nlevels = 10;

gray = 0.7*ones(1,3);
col_neg = [0.2622    0.6028    0.7112]; % dark teal
col_pos = [0.9961    0.0782    0.4427]; % fuschia
col_zero = ones(1,3);

col_map_corr  = makeColorMap(col_neg,col_zero,col_pos,nlevels); 


%----and now make a plot!


figure(1),clf

  %-----------------------
  
  ax = axesm('MapProjection','eqdcylin','grid','on',...
      'MeridianLabel','on','ParallelLabel','on',...
      'PLabelLocation',30,'MLabelLocation',90);
   
  % NOTE: sometgimes the plot disappears when we add a grid.
  % this prob seems to only be present in matlab2010b (not a)
      
  contourfm(lat,lon,TF',nlevels,'LineStyle','none');
  caxis([-1 1]);
  colormap(col_map_corr);

  setm(ax,'Fontsize',16);
  setm(ax,'Fontcolor','Black');
  c = load('coast');
  plotm(c.lat,c.long,'Color','black','LineWidth',1);

  %  For ps, mask out the sea.  This is done with a weird trick...
  if strcmp(v,'m')
    geoshow(flipud(c.lat), flipud(c.long),'DisplayType', 'polygon', 'FaceColor', gray)
  end
  
  % add a title
  title(['\chi_',num2str(comp),' TF for  ',v])
  
  % clean up a little bit
  box off
  axis off
  grid on
  h = colorbar; % Make colorbar and save handle.
  
  
  
%----export plots

plot_dir = '/home/ig/neef/Documents/Plots/AAM/';

%--- figure settings
LW = 3;
ph = 8;        % paper height
pw = 10;        % paper width
FS = 20;        % fontsize


fig_name = [plot_dir,'TF_slice',slice_type,'_',v,'_chi',num2str(comp),'.png'];
exportfig(1,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');





