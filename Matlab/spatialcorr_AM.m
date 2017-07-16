% spatialAMcorr.m----------------
%
%  spatialAMcorr(term,component,TT)
%  Compute spatial correlations between different components of AAM and the total on a certain timescale.
%  started 20 April 2010  
%  mods:
%	 convert to function and make customizable, 22 April 2010
%	 add the option to compute variances instead of correlations, 11 May 2010
%
%  INPUTS:
%	term: choose m for mass, w for motion
%	comp: choose 1 for equatorial or 2 for axial
%	TT: vector or min and max periods on which to focus
%	covar: set to 1 to compute covariances, 0 to compute correlations
%---------------------------------------------------------
function [Rout,LAG,lon,lat] = spatialcorr_AM(term,comp,TT,covar)

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
%load spatialAM_1978_1989.mat
load spatialAM_1989_2000.mat
[nv,nlon,nlat,nt] = size(xxm);
%lat2 = mean(lband,2);

% combine equatorial EFs to two complex functions
% and convert to equatorial terms to mas while converting axial term to LOD changes
% while also removing the mean in each case
dum = squeeze(xxm(3,1,1,:));
g = find(isfinite(dum) == 1);
meanm = mean(xxm(:,:,:,g),4);
meanw = mean(xxw(:,:,:,g),4);
XXw = zeros(nlon,nlat,nt);
XXm = zeros(nlon,nlat,nt);
for it = 1:nt
  if comp == 2
    XXm(:,:,it) = LOD0*(xxm(3,:,:,it) - meanm(3,:,:));
    XXw(:,:,it) = LOD0*(xxw(3,:,:,it) - meanw(3,:,:));
  else
    XXm(:,:,it) = r2m*squeeze(xxm(1,:,:,it)-meanm(1,:,:)+i*xxm(2,:,:,it)-i*meanm(2,:,:));
    XXw(:,:,it) = r2m*squeeze(xxw(1,:,:,it)-meanw(1,:,:)+i*xxw(2,:,:,it)-i*meanw(2,:,:));
  end
end


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
if term == 'm'
  y = XXm(:,:,1:nt2);
else
  y = XXw(:,:,1:nt2);
end
Y = zeros(1,nt2);

if comp == 2
  for kk = k0:kf, Y(kk-k0+1) = M(3,kk); end
else
  for kk = k0:kf, Y(kk-k0+1) = M(1,kk)+i*M(2,kk); end
end

% also compute the global contribution in time via simple addition - this removes the area weight
dum = XXw+XXm;
dxyp = gridcellarea(nlat,nlon);
%dum2 = dum;
%for ilon = 1:128
%  for it = 1:2191
%    dum2(ilon,:,it) = squeeze(dum(ilon,:,it))./dxyp';
%  end
%end
Y2 = squeeze(sum(sum((dum(:,:,1:nt2)),1),2));

%---cycle through points and compute covariances.

R = zeros(nlon,nlat);
R2 = zeros(nlon,nlat);
LAG = zeros(nlon,nlat);
dxyp = gridcellarea(nlat,nlon);



% slap everything on the same timeaxis.  Call these new functions Y
k0 = find(MJD2 == min(MJD)+0.5);
kf = length(MJD2);
nt2 = kf-k0+1;
for ilat = 1:nlat
  for ilon = 1:nlon
    % divide out gridcell area to account for this effect in the integral
    yt = squeeze(y(ilon,ilat,:));
    %yt2 = squeeze(y(ilon,ilat,:))/dxyp(ilat);
    g = find(isfinite(yt) == 1);
    if covar == 1 
      [c,lags] = xcov(yt,Y,'none');
    else
      [c,lags] = xcorr(yt,Y,'coeff');
    end
    cn = sqrt(imag(c).^2+real(c).^2);
    [a,b] = max(cn);
    %R(ilon,ilat) = a;
    R(ilon,ilat) = c(find(lags == 0));
    LAG(ilon,ilat) = lags(b);
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
%lat = lat2;
