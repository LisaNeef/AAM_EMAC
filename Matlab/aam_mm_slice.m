function [lat,lon,X,u300,ps] = aam_mm_slice(comp,runid,twenty_years)
% aam_mm_slice--------------
%
% Take existing AAM mat file for an EMAC run and return monthly-mean
% AAM contribution fields for each month.
% Here monthly values are computed by averaging over the entire run.
% SEE NOTES VOL 4 P 2
% Mods:
%  8 Nov 2010 - make it  so that we can also get fields for ERA
% 18 Nov 2010 - option to only use last 20 years to make it
%	comparable with ERA

% load the data 
sdir = ['/dsk/nathan/lisa/ex/',runid,'/mat/'];
if runid(1:3) == 'ERA', sdir = ['/dsk/nathan/lisa/ERAinterim/mat/']; end
ff =  dir([sdir,'*comp*.mat']);
if twenty_years, ff =  dir([sdir,'*comp*last20only.mat']); end
fname = [sdir,ff(1).name]
load(fname)
[dum,nm,nlat,nlon] = size(Xm_mean);

% reorder the longitudes
if runid(1:3) ~= 'ERA'
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
X = zeros(2, nm,nlat,nlon);
for im = 1:nm
  X(1,im,:,:) = pw*uc*squeeze(Xw_mean(comp,im,:,:) - nanmean(Xw_mean(comp,:,:,:),2));
  X(2,im,:,:) = pm*uc*squeeze(Xm_mean(comp,im,:,:) - nanmean(Xm_mean(comp,:,:,:),2));
end

u300 = u300_mean;
ps = ps_mean;

