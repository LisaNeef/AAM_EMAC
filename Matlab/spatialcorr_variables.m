% spatialcorr_variables.m----------------
%
%  spatialcorr_variables.m(var,component,TT)
%  Compute spatial correlations between  emac model fields and AAM filtered for a certain timescale.
%  started 26 April 2010  
%  Mods:
%	7 july 2010 - TO DO: add option to do covariance or correlation
%
%  INPUTS:
%	variable: choose u or v for respective 300hPa winds, or p for surface pressure
%	TT: vector or min and max periods on which to focus
%	comp: component (1 for equatorial, 2 for axial)
%	do_covar: set to 1 for covariance, 0 for correlation
%---------------------------------------------------------
function [Rout,LAG,lon,lat] = spatialcorr_variables(comp,TT,variable,do_covar)

%---paths and file settings
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/m_map/');

%---compute AAM prefactors (see Barnes et al., 1983, p. 48)

d2r = pi/180.;         %  degrees to radian
LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate
r2m = (180/pi)*60*60*1000;    % conversion of radians to mas

%--load data: chi terms in space and time.
%load spatialAM_1995_2000.mat
load emac_uvpfields_ccmval.mat
[nlon,nlat,nt] = size(U300);

%---retrieve the global AAM and filter out the desired timescale

[t,M,O,R,G,H,MJD2] = aam_compare_tscales(TT(1),TT(2),'IERS','ERA');


%---slap everything on the same timeaxis.  Call these new functions Y
k0 = find(MJD2 == min(MJD)+0.5);
if max(MJD2) <= max(MJD)
  kf = length(MJD2);
else
  kf = find(MJD2 == max(MJD)+0.5);
end
nt2 = kf-k0+1;
Y = zeros(2,nt2);
for kk = k0:kf
  Y(2,kk-k0+1) = M(3,kk);
  Y(1,kk-k0+1) = M(1,kk)+i*M(2,kk);
end

if variable == 'u', y = U300(:,:,1:nt2); end
if variable == 'v', y = V300(:,:,1:nt2); end
if variable == 'p', y = Ps(:,:,1:nt2);  end

%---cycle through points and compute correlations.
R = zeros(nlon,nlat);
LAG = zeros(nlon,nlat);


for ilat = 1:nlat
  for ilon = 1:nlon
    yt = squeeze(y(ilon,ilat,:));
    yt2 = yt-mean(yt);
    if do_covar == 1 
      [c1,lags1] = xcov(yt2,Y(comp,:),'none');
    else
      [c1,lags1] = xcorr(yt2,Y(comp,:),'coeff');
    end
    if comp == 1; c = sqrt(imag(c1).^2+real(c1).^2); else c = c1; end
    [a,b] = max(c);
    R(ilon,ilat) = a;
    LAG(ilon,ilat) = lags1(b);
  end
end

% also straighten out longitude and sort.
dum = find(lon > 180);
lon1 = lon;
lon2 = lon;
lon2(dum) = lon(dum)-360;
[a,b] = sort(lon2);

Rout = R(b,:);
lon = lon2(b);
