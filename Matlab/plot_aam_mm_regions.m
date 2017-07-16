% plot_aam_contrib_mm--------------
%
% Examine the regional contributions to AAM in EMAC output, 
% by making monthly-mean maps of:
%	zonal wind at 300 hPa
%	surface pressure
%	and the local components of Xm and Xw
% 25 Aug 2010 
% Mods:
%  21 Oct 2010: make more flexible for other datasets
%  	and make it so that the constants are loaded from a function.
%
%-----------------------------------------------------------------
clear all;

% choose what to plot:
%pp  = 1;	disp('plotting monthly mean surface pressure') 
pp  = 2; 	disp('plotting monthly mean pressure term of X1')
%pp  = 3;	disp('plotting monthly mean  pressure term of X2')
%pp  = 4; 	disp('plotting monthly mean  pressure term of X3')
%pp  = 5;	disp('plotting monthly mean  300 hPa u')
%pp  = 6;	disp('plotting monthly mean  vertically-integrated wind term of X1')
%pp  = 7;	disp('plotting monthly mean  vertically-integrated wind term of X2')
%pp  = 8; 	disp('plotting monthly mean  vertically-integrated wind term of X3')


% choose whether to plot anomaly or absolute field
anom = 1;	disp('plotting anomalies')
%anom = 0; 	disp('plotting absolute field')

%---paths...

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/bipolar_colormap');


% load the data 
CCMval = 0;
%runid = 'ref2_T31L39'
runid = 't7_T42L39'
sdir = ['/dsk/nathan/lisa/ex/',runid,'/mat/'];
ff =  dir([sdir,'*comp*.mat']);
fname = [sdir,ff.name]
load(fname)

nm = length(months);

% reorder the longitudes
dum = find(lon > 180);
lon2=lon;
lon2(dum) = -(360 - lon(dum));
[a,b] = sort(lon2);
lon = a;

% here are lat and lon in radians
d2r = pi/180.;         %  degrees to radian
rlat = lat*d2r;
rlon = lon*d2r;


% if it's a CCMVal run, also reorder latitudes
%if CCMval == 1
%  [c,d]=sort(lat);
%else
  nlat = size(lat);
  d = 1:nlat;
%end

%---for X plots, we also need prefactors - load from predefined.
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)

% prefactors (Gross 09)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));
pm3  = 0.997*0.750*Re^4/(g*Cm);
pw3  = 0.997*Re^3/(Q*g*Cm);
pm = [pm12; pm12; pm3];
pw = [pw12; pw12; pw3];

% this is the older form (Barnes 83)
%pm12 = -Re^4/(g*CminusA);
%pw12 = -1.43*Re^3/(Q*g*CminusA);
%pm3  = 0.7*Re^4/(g*C);
%pw3  = Re^3/(Q*g*C);


% select which option to plot
if pp == 1, R = (1e-2)*ps_mean(:,d,b); end				% units hPa
if pp == 2, R = rad2mas*pm12*squeeze(Xm_mean(1,:,:,b)); end		% units mas
if pp == 3, R = rad2mas*pm12*squeeze(Xm_mean(2,:,:,b)); end		% units mas
if pp == 4
  Xm_mean_anom = zeros(nm,64,128);
  for im = 1:nm, Xm_mean_anom(im,:,:) = squeeze(Xm_mean(3,im,:,b)-mean(Xm_mean(3,:,:,b),2)); end
  R = LOD0_ms*pm3*Xm_mean_anom; %units microsec
end
if pp == 5, R = u300_mean(:,:,b); end					% units m/s

if pp == 6, R = rad2mas*pw12*squeeze(Xw_mean(1,:,:,b)); end		% units mas
if pp == 7, R = rad2mas*pw12*squeeze(Xw_mean(2,:,:,b)); end		% units mas
if pp == 8
  Xw_mean_anom = zeros(nm,64,128);
  for im = 1:nm, Xw_mean_anom(im,:,:) = squeeze(Xw_mean(3,im,:,b)-mean(Xw_mean(3,:,:,b),2)); end
  R = LOD0_ms*pw3*Xw_mean_anom; %units microsec
end


% if anom option selected, compute anomalies
if anom == 1
  D = R*0;
  Rmean = mean(R(:,:,:),1);
  [nm,nlat,nlon] = size(R);
  for im = 1:nm
    D(im,:,:) = R(im,:,:)-Rmean;
  end
else
  D = R;
end


% plot settings
anom
if pp == 1
  if anom == 1, cax = 10*[-1,1]; else cax = [500,1100]; end
  T = 'IB-Corrected Monthly-Mean Ps (hPa): ';
end
if pp == 2
  cax = 30*[-1,1];
  T = 'IB-Corrected Monthly-Mean \chi _1^m anomaly (mas): ';
end
if pp == 3
  cax = 30*[-1,1];
  T = 'IB-Corrected Monthly-Mean \chi_2^m anomaly (mas): ';
end
if pp == 4
  cax = 200*[-1,1];
  T = 'IB-Corrected Monthly-Mean \Delta LOD anomaly (\mu s): ';
end
if pp == 5
  cax = 25*[-1,1];
  T = 'Monthly-Mean Zonal Wind at 300 hPa (m/s): ';
end
if pp == 6
  cax = 400*[-1,1];
  T = 'Monthly-Mean \chi _1^w (mas): ';
end
if pp == 7
  cax = 400*[-1,1];
  T = 'Monthly-Mean \chi _2^w (mas): ';
end
if pp == 8
  cax = 400*[-1,1];
  T = 'Monthly-Mean \chi _3^w in terms of LOD (\mu s): ';
end

figure(1),clf
iplot=1;
%for im = 1:1:nm
for im = 2:2:nm
  figure(iplot),clf
  axesm robinson; framem;gridm
  c = load('coast');
  plotm(c.lat,c.long,'Color',0*ones(3,1))
  %contourfm(lat,lon,squeeze(D(:,:)),30);
  contourfm(lat,lon,squeeze(D(im,:,:)),30);
  caxis(cax)
  if anom == 1
    title([T 'Anomaly, Month ' num2str(months(im))])
    colormap(bluewhitered(256)), colorbar('SouthOutside')
  else
    title([T 'Month ' num2str(months(im))])
    if pp == 1, colormap(jet(256)),  end
    if pp == 5|6, colormap(bluewhitered(256)),  end
    colorbar('SouthOutside')
  end
  iplot = iplot + 1;
end


% as a check, for AAM functions, integrate globally to get the apprioximate monthly anomaly
R_dlon = trapz(rlon, R,3);
R_dlon_dlat = -trapz(rlat, R_dlon,2);
if pp == 2|3
  disp(['Global integrated Xm for Feb (mas):  ',num2str(R_dlon_dlat(1))])
end

save temp_X1monthly_anomaly_t7_T42L39.mat D pp lat lon months 
