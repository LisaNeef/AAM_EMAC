%---emac_uvpfields.m-------------------------:
% cycle through emac output files and compile arrays of u, v, p, and mslp.
% see notes vol. 2 p. 81
% started 22 April 2010 - see notes vol 2 p 84
%  mod 10 may 2010 - also pull out mslp field
%
% INPUTS:
%	y0: start year
%	yf: stop year

% choose desired years:
y0 = 1995;
yf = 2000;

%---paths and file settings
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/m_map/');

dir = '/dsk/nathan/lisa/Data/CCMval/';
pref_uv = 'scout02____dm_';
suff_uv = '_uv.nc';
pref_aps = 'O3_';
suff_aps = '_aps.nc';
pref_mslp = 'scout02____dm_mslp_';
suff_mslp = '.nc';

%---start and stop times and time axis
mm = 1:1:12;
day0 = date2mjd(y0,1,1,0,0,0);
dayf = date2mjd(yf,12,31,0,0,0);
ndays = dayf - day0;

% (these are the reference dates from which each file starts counting)
refday_uv = date2mjd(1958,1,1,0,0,0);
refday_aps = date2mjd(1901,1,15,0,0,0);

%---initialize arrays 
nlon = 128;
nlat = 64;
U300 = zeros(nlon,nlat,ndays)+NaN;
V300 = zeros(nlon,nlat,ndays)+NaN;
Ps = zeros(nlon,nlat,ndays)+NaN;
MSLP = zeros(nlon,nlat,ndays)+NaN;
MJD = zeros(1,ndays)+NaN;

%---cycle through netcdf files by year and month and store the fields
kold = 0;
for y = y0:1:yf
  for m = 1:12 
    mstr = num2str(m);

    if size(mstr,2)==1, mstr=['0',mstr]; end
    datestr = [num2str(y),mstr];
    ff_uv = [pref_uv,num2str(y),mstr,suff_uv];
    ff_ps = [pref_aps,num2str(y),mstr,suff_aps];

    if exist([dir,ff_uv]) + exist([dir,ff_ps]) == 4 

      % load wind data
      % variables 3, 4, 5, and 6 are, respectively, level, time, u, v
      [data name] = r_netcdf(dir,ff_uv);
      nt = length(data{4});
      k1 = kold+1;
      k2 = k1+nt-1;
      lon = data{1};
      lat = data{2};
      lev = data{3};
      MJD(k1:k2) = data{4}+refday_uv;
      u = data{5};
      v = data{6};

      % store winds at 300 hPa.
      U300(:,:,k1:k2) = squeeze(u(:,:,find(lev == 300),:));
      V300(:,:,k1:k2) = squeeze(v(:,:,find(lev == 300),:));

      % load surface pressure data 
      % here variables 3 and 4 are time and surface pressure
      [data name] = r_netcdf(dir,ff_ps);
      nt_aps = length(data{3});
        MJD(k1:k2) = data{3}+refday_aps;	% the reference date is different for aps data
      ps = data{4};
      Ps(:,:,k1:k2) = squeeze(ps(:,:,:));

      kold = k2;
      disp(['Finished   ',num2str(y),'-',mstr])

    end % file exist if loop

  end	% month loop
end	% year loop




% save the big state arrays
  save 'emac_uvpfields.mat' U300 V300 Ps y0 yf  MJD lon lat

