function D = compute_aam_regions(runid,vector_comp,term)
% compute_aam_regions--------------
%
% Function that returns the monthly-mean AAM contribution per 
%  model gridcell, in EAF function units, by loading an existing mat file of AAM 
%  monthly-mean contributions.
% This function can, e.g., be used by plot_aam_mm_latprofiles_hres.m 
%
% INPUTS
%	runid: the name of the run we want
%	vector_comp: which EAF?  X1, X2, or X2
%	term: choose 1 for wind and 2 for pressure
% 26 October 2010
%
%-----------------------------------------------------------------


% ---load the data 
sdir = ['/dsk/nathan/lisa/ex/',runid,'/mat/'];
ff =  dir([sdir,'*comp*.mat']);
fname = char([sdir,ff.name])
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

%---load EAF prefacotors
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)
% prefactors (Gross 09)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));
pm3  = 0.997*0.750*Re^4/(g*Cm);
pw3  = 0.997*Re^3/(Q*g*Cm);

% here is a matrix ordered by vector_comp (3 possibilities) and term (2 possibilities)
P = zeros(3,2);
P(1,1) = pw12;	P(1,2) = pm12;
P(2,1) = pw12;	P(2,2) = pm12;
P(3,1) = pw3;	P(3,2) = pm3;


% now from the monthly means, subtract out the long-term mean
[nt,nm,nlat,nlon] = size(Xm_mean);
XX = zeros(2,nt,nm,nlat,nlon)+NaN;
XX(1,:,:,:,:) = Xw_mean;
XX(2,:,:,:,:) = Xm_mean;

% ----compute EAFs based on the plotting option
PICK UP HERE!!!  -- OR RATHER: CHANGE THE CODE THAT COMPUTES MONTHLY MEANS TO ALSO TAKE OUT THE LT MEAN.


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
