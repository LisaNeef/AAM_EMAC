% ----compare_sst_aam.m--------------------
% Read in Hadley SSTs (Emac input), compute variability in time, and compare to 
%  the variability of AAM computed from EMAC (as in Rosen & Salstein 2000, Fig 3)
% started 15 Feb 2010
% see notes vol. 2 p. 19
% modifications to this - notes vol. 2 p. 34

clear all;

%--- paths and other settings

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');

dir = '/dsk/nathan/lisa/Data/SCOUT/MESSy_input/SCOUT/HADLEY/T42/';
pref = 'T42_amip2sst_';
suff = '.nc';

%---initialize arrays.

year = 1950:1:2006;
nf = length(year);
SSTm = zeros(nf,12);
dum_MJD = zeros(nf,12);

%---SSTs: load first file for lat / lon arrays

y = '1950';
ff = [pref,y,suff];
[data name] = r_netcdf(dir,ff);
n = size(data,2);
for i = 1:n
  eval([name{i},'=data{i};',]);
end

% ENSO3 region definition
enso3_lat = find(abs(lat) <= 5);
enso3_lon = find(lon >= 210 & lon <= 270);

%---cycle through SST files and read out SSTs, average, and store

year = 1950:1:2006;

for ifile = 1:nf
  y = num2str(year(ifile));
  ff = [pref,y,suff];
  [data name] = r_netcdf(dir,ff);
  n = size(data,2);
  for i = 1:n
    eval([name{i},'=data{i};',]);
  end
  % for each month, average the SST over the ENSO 3 regions
  dum = sst(enso3_lon,enso3_lat,:);
  SSTm(ifile,:) = (squeeze(mean(mean(dum))))';
  dum_MJD(ifile,:) = date2mjd(year(ifile),1:12,0,0,0,0);
end

%---create a time series of SST anomalies
% ...as well as a time array in MJD
SSTclim = mean(SSTm);
SSTanom = SSTm - ones(nf,1)*SSTclim;
SSTA = zeros(1,nf*12)+NaN;
MJD_SST = SSTA*0;
for iyear = 1:nf
  k1 = 12*iyear-11;
  if iyear == 1, k1 = 1; end
  k2 = 12*iyear ;
  SSTA(k1:k2)=SSTanom(iyear,:);
  MJD_SST(k1:k2)=dum_MJD(iyear,:);
end

%--- cycle through anomalies and compute running standard deviation
k0 = find(year == 1960);
kf = find(year == 1996);
nmonths = length(SSTA);
w = (21/2)*12;		% half of a 21-year window in terms of months
k0 = 1+w;
kf = nmonths - w;
S = SSTA*0+NaN;
for k = k0:kf
  k1 = k-w;
  k2 = k+10;
  S(k) = var(SSTA(k1:k2));
end

% following Rosen & Salstein 2001, divide this series by its own STD and subtract its mean
g1 = find(isfinite(S) == 1);
NS = S/std(S(g1));

%---restore AAM functions computed from EMAC

load aam_emac_1960_2000.mat
g = find(isfinite(X(3,:)) == 1);
X3 = X(3,g);
nt = length(X3);
MJD_emac = MJD(g);

wd = round((21/2)*365);		% half of a 21-year window in terms of days
k0 = 1+wd;
kf = nt - wd;

SX3 = zeros(nt,1)+NaN;

% also compute STD in time.
for k = k0:kf
  k1 = k-wd;
  k2 = k+wd;
  SX3(k) = var(X3(k1:k2));
end

g2 = find(isfinite(SX3) == 1);
NSX3 = SX3/std(SX3(g2));


% ---make a plot comparing normalized standard deviations in both tseries
figure(1),clf
  sst = plot(MJD_SST,NS-mean(NS(g1)),'LineWidth',2,'Color',rand(3,1));
  hold on
  emac = plot(MJD_emac,NSX3-mean(NSX3(g2)),'k','LineWidth',2);
  xlabel('time')
  title('Standardized Running Variance NINO-3 Regional Mean SST Anomaly')
  ylabel('Normalized STD')
  lhandle=[sst(1) emac(1)];
  legend(lhandle,'SST','EMAC X_3')



