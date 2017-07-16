function [X_mm,Y] = aam_geo_mm
%------aam_era_mm.m-----------------------------------------------------
%  Read in the geodeitcally-implied excitation functions
%  and compute monthly mean values in order to examine the annual
%  cycle.
%  See notes vol. 4, p. era 17
%  Started 12 Nov 2010.
%
%  Mods:
%   7 Jul 2011: updated call to read_eops for modified code.
%-----------------------------------------------------------------------


%---paths
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');


%---read in the geodetic data
[X1,X2,dlod,MJD] = read_eops;

[y, m, d, h] = mjd2date(MJD);
y0 = min(y);
yf = max(y);
nyears = yf-y0+1;
nm = 12;
year = y0:yf;
month = 1:12;

%--pack components together and detrend.
aam_constants_gross

X = [X1';X2';dlod'];
X_dt = X*0;

for ii=1:3, X_dt(ii,:) = detrend(X(ii,:),'constant'); end

%---cycle years and for each year, compute monthly mean

X_mm = zeros(3,nm, nyears);
Y = zeros(nyears);

for iy = 1:nyears
  Y(iy) = year(iy);
  for im = 1:nm
    tt = find((y == year(iy)) & m == month(im));
    X_mm(:,im,iy) = mean(X_dt(:,tt),2);
  end % loop over months
end   % loop over years





