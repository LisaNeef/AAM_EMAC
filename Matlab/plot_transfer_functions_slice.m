%% plot_TF_lat
%
% A simple plot showing the latitude slice of one of the EAM transfer
% functions.
%
% Started 25 Nov 2011
%
%-------------------------------------------------------------


%% user inputs

comp    = 3;
v     = 'u';
slice_type = 'lat';


%% compute the appropriate function


% generate a lat / lon grid
lat = -90:90;
lon = 0:360;

d2r = pi/180;
rlat = lat*d2r;
rlon = lon*d2r;

[LAT,LON] = meshgrid(rlat,rlon);

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
        T = 'cos^2 \phi';
      case 'v'
        disp('Variable v does not affect the axial excitation function at all!')
        return
      case 'm'
        TF = cos(LAT).*cos(LAT).*cos(LAT);
    end 
end


nlat = length(lat);
nlon = length(lon);
mlon = floor(nlon/2);
mlat = floor(nlat/2);

switch slice_type
    case 'lat'
        xx = squeeze(TF(mlon,:));
    case 'lon'
        xx = squeeze(TF(:,mlat));
end        


%% plot settings

LW = 2;
col = zeros(1,3);

%% plot


figure(1),clf

switch slice_type
    case 'lat'
        plot(xx,lat,'LineWidth',LW,'Color',col);     
        ylabel('latitude') 
        title(T)
    case 'lon'
        disp('still need to write out the code for lon profiles!')
end

grid on



%% export



plot_dir = '/home/ig/neef/Documents/Plots/';

%--- figure settings
LW = 1;
ph = 10;        % paper height
pw = 8;        % paper width
FS = 20;        % fontsize


fig_name = [plot_dir,'TF_',v,'_chi',num2str(comp),'.png'];
exportfig(1,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');




