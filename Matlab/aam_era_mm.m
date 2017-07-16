function [Xw_mm,Xm_mm,Y] = aam_era_mm(AM,era_type)
%------aam_era_mm.m-----------------------------------------------------
%  Read in the AAM excitation functions from ERAInterim (or ERA-40)
%  and compute monthly mean values in order to examine the annual
%  cycle.
%  See notes vol. 3 p. 124
%	and vol. 4, p. 117
%  Inputs:
%	AM type (OAM, HAM, AAM)
%	era_type: (1) for ERA-40, (2) for Interim
%
%  Mods: 
%    12 Nov 2010 - add possibility of reading out HAM and OAM as well
%     2 Mar 2011 - add choice of reading in ERA-40 or Interim.
%-----------------------------------------------------------------------


%---paths
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');


%---read in the era data
switch era_type
  case 1
    [xw,xm,MJD]= read_EFs(AM,'ERAinterim',0);
  case 2
    [xw,xm,MJD]= read_EFs(AM,'ERA40',0);
end

[y, m, d, h] = mjd2date(MJD);
y0 = min(y);
yf = max(y);
nyears = yf-y0+1;
nm = 12;
year = y0:yf;
month = 1:12;

%--next, detrend each row and convert units to mas (X1 and X2) or ms (X3)
aam_constants_gross

for ii = 1:3
  xm_dt(ii,:) = detrend(xm(ii,:),'constant');
  xw_dt(ii,:) = detrend(xw(ii,:),'constant');
end

xm_scaled = xm*0;
xw_scaled = xm*0;

xw_scaled(1:2,:) = rad2mas*xw_dt(1:2,:);
xw_scaled(3,:) = LOD0_ms*xw_dt(3,:);
xm_scaled(1:2,:) = rad2mas*xm_dt(1:2,:);
xm_scaled(3,:) = LOD0_ms*xm_dt(3,:);

%---cycle years and for each year, compute monthly mean

Xm_mm = zeros(3,nm, nyears);
Xw_mm = zeros(3,nm, nyears);
Y = zeros(nyears);

for iy = 1:nyears
  Y(iy) = year(iy);
  for im = 1:nm
    tt = find((y == year(iy)) & m == month(im));
    Xm_mm(:,im,iy) = mean(xm_scaled(:,tt),2);
    Xw_mm(:,im,iy) = mean(xw_scaled(:,tt),2);
  end % loop over months
end   % loop over years





