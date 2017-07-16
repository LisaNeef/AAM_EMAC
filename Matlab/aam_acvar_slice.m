function [lat,lon,varX,varU300,varPS] = aam_acvar_slice(comp,runid)
% aam_acvar_slice--------------
%
% Take existing AAM mat file for an EMAC run and return 
% the annual variance for each point in lat and lon, of the 
% AAM gridpoint contributions (integrated vertically), and u300 and ps
% gridpoint contributions.
% 17 Jan 2011  -- see notes vol. 4, p. 69.
%
% INPUTS:
%  comp: the desired AAM vector component.
%  runid: ID of the model run to investigate.
%
% OUTPUTS:
%  lat, lon arrays
%  varX: variance in individual aam contributions (wind, mass, total), 
%	for vector component comp
%  varU300, varPS: variances for these variable components.

% load the data 

if strcmp(runid,'CCMVal'), runid = 'CCMval'; end

sdir = ['/dsk/nathan/lisa/EMAC_ex/',runid,'/mat/'];
if strcmp(runid,'ERAint'), sdir = '/dsk/nathan/lisa/ERAinterim/mat/'; end
ff =  dir([sdir,'*comp*allyears*.mat']);
fname = [sdir,ff(1).name];

load(fname)
[~,nm,ny,nlat,nlon] = size(Xm_mean);

% reorder the longitudes
if strcmp(runid,'ERA') == 0
  dum = find(lon > 180);
  lon2=lon;
  lon2(dum) = -(360 - lon(dum));
  [a,b] = sort(lon2);
  lon = a;
end

% here are lat and lon in radians
d2r = pi/180.;         %  degrees to radian
rlat = lat*d2r;
rlon = lon*d2r;

%---load the right prefactors
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)

% prefactors (Gross 09)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));
pm3  = 0.997*0.750*Re^4/(g*Cm);
pw3  = 0.997*Re^3/(Q*g*Cm);
P = zeros(2,3);
P(2,:) = [pm12; pm12; pm3];
P(1,:) = [pw12; pw12; pw3];

% select the appropriate prefactor and unit conversion
pw = P(1,comp);
pm = P(2,comp);
if comp == 3
  uc = LOD0_ms;		% conversion to milliseconds for axial term (LOD)
else
  uc = rad2mas;		% conversion to mas for equatorial terms
end

% subtract out the LT mean and simultaneously convert units
dX = zeros(3,nm,ny,nlat,nlon);
dU300 = zeros(3,nm,ny,nlat,nlon);	% 2 dummy rows in the 1st dimension to make arrays similar
dPS = zeros(3,nm,ny,nlat,nlon);

Xw_LT = squeeze(nanmean(nanmean(Xw_mean(comp,:,:,:,:),2),3));
Xm_LT = squeeze(nanmean(nanmean(Xm_mean(comp,:,:,:,:),2),3));
u300_LT = squeeze(nanmean(nanmean(u300_mean(:,:,:,:),1),2));
ps_LT = squeeze(nanmean(nanmean(ps_mean(:,:,:,:),1),2));

for im = 1:nm
  for iy = 1:ny
    dX(2,im,iy,:,:) = pw*uc*(squeeze(Xw_mean(comp,im,iy,:,:)) - Xw_LT); 
    dX(3,im,iy,:,:) = pm*uc*(squeeze(Xm_mean(comp,im,iy,:,:)) - Xm_LT); 
    dU300(1,im,iy,:,:) = squeeze(u300_mean(im,iy,:,:))-u300_LT;
    dPS(1,im,iy,:,:) = squeeze(ps_mean(im,iy,:,:))-ps_LT;
  end
end
dX(1,:,:,:,:) = dX(2,:,:,:,:) +  dX(3,:,:,:,:);


% compute the annual variances in the terms
% first compute variance over year, for all the years
% --> then plot average annual variance over run. (see notes vol. 4, p. 69)



varX = squeeze(nanmean((dX.^2),2));
varU300 = squeeze(nanmean((dU300.^2),2));
varPS = squeeze(nanmean((dPS.^2),2));

