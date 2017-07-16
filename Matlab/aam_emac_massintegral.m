%function aam_emac_massintegral(runid,decade)
%---aam_emac.m-------------------------:
% cycle through EMAC output files and compute atmospheric angular momentum
% in this version the integral is computed over mass elements, which allows us to use
% the model output on hybrid levels.
% Lisa Neef
% started 19 Jul 2010
% see notes vol. 3 p. 26, 40
% mods:

%---temporary inputs-----------
clear all;
runid = 't7_T42L39';
dec = 1980;
savecomp = 1;

%---paths and file settings

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');

decade = num2str(dec);
datadir = ['/dsk/nathan/lisa/ex/',runid,'/uvm/'];
flist = dir([datadir,'*',decade(1:3),'*messy_uvm.nc']);
nf = size(flist,1);

%---define some constants

Re = 6.371e6;    	% radius of earth in meters
Q = 7.292115e-5; 	% rotation rate of the earth in rad/s
g = 9.81;        	% grav constant in m/s2
C = 7.1236e37;  	% principal moment of inertia of the mantle in kgm2
CminusA = 2.34e35;       % C minus eq. MOI, same units 
d2r = pi/180.;         %  degrees to radian

%---compute AAM prefactors (see Barnes et al., 1983, p. 48)
p12 = 1/(Q*CminusA);
p3  = 1/(Q*C);
prefac_I = Q*([p12; p12; 0.698*p3]);
prefac_h = [1.43*p12; 1.43*p12; p3];

%---load landmask data (for IB correction)
ff_lm = '/home/ig/neef/Data/sample_landmask.nc';
slm = nc_varget(ff_lm,'slm');
lat = nc_varget(ff_lm,'lat');
lon = nc_varget(ff_lm,'lon');
rlat = lat*d2r;
rlon = lon*d2r;

ref_date_str = '19580101';
yr = str2num(ref_date_str(1:4));
mr = str2num(ref_date_str(5:6));
dr = str2num(ref_date_str(7:8));
ref_day = date2mjd(yr,mr,dr,0,0,0);		% initial date in avail. files

%---set up grids and time axis

f0 = [datadir,flist(1).name];
ff = [datadir,flist(nf).name];

dm = nc_varget(ff,'grmass');
[ntf,nlev,nlat,nlon] = size(dm);

% specify initial date of run
time0 = nc_varget(f0,'time');
day0 = time0(1)+ref_day;

% final date of the times considered here.
timef = nc_varget(ff,'time');
dayf = timef(ntf)+ref_day;

ndays = dayf - day0;				% number of days in avail. files
opd = round(ntf/30);				% number of obs per day

% matricies of lats and lons that are the same grid as the u and v fields
% (make this array bigger than 31 days to leave from for any NaNs in the output)
RLAT31 = zeros(33*opd,nlev,nlat,nlon);
RLON31 = RLAT31;
for ilev=1:nlev
  for it = 1:33*opd
    RLAT31(it,ilev,:,:) = rlat*ones(1,128);
    RLON31(it,ilev,:,:) = ones(64,1)*rlon';
  end
end

%---initialize arrays for AM and LOD

ww = zeros(3,ndays*opd,nlat,nlon)+NaN;
mm = zeros(3,ndays*opd,nlat,nlon)+NaN;
Xw = zeros(3,ndays*opd)+NaN;
Xm  = zeros(3,ndays*opd)+NaN;
MJD = zeros(1,ndays*opd)+NaN;


%---cycle through netcdf files by year and month and perform AAM integrals
kold = 0;

for ifile = 1:nf
  f = [datadir,flist(ifile).name];

  if exist(f) == 2  
    % compute AAM only if file exists and u has a time dimension
    % (i.e. is not screwed up somehow)
    ucoslat = nc_varget(f,'um1');
    vcoslat = nc_varget(f,'vm1');
    u = ucoslat*0;  v = vcoslat*0;

    % if the time-dimension isn't screwed up here, do the integral
    if size(size(u),2)==4

      % note that um1 and vm1 are multiplied by cos(lat)
      for ilat = 1:nlat  
        u(:,:,ilat,:) = ucoslat(:,:,ilat,:)/cos(rlat(ilat));
        v(:,:,ilat,:) = vcoslat(:,:,ilat,:)/cos(rlat(ilat));
      end

      dm = nc_varget(f,'grmass');
      time = nc_varget(f,'time');
      nt = size(dm,1);
      k1 = kold+1;
      k2 = k1+nt-1;

      no_data = find(time == 0);
        u(no_data) = NaN;
        v(no_data) = NaN;
        mjd(no_data) = NaN;
        dm(no_data) = NaN;
        time(no_data) = NaN;
      MJD(k1:k2) = time+ref_day;

      RLAT = RLAT31(1:nt,:,:,:);
      RLON = RLON31(1:nt,:,:,:);

      % ***insert inverse barometer correction here?

      h_field = zeros(3,nt,nlev,nlat,nlon);    
      I_field = zeros(3,nt,nlev,nlat,nlon);    

      h_field(1,:,:,:,:) = -Re*(u.*sin(RLAT).*cos(RLON)-v.*sin(RLON)).*dm;
      h_field(2,:,:,:,:) = -Re*(u.*sin(RLAT).*sin(RLON)+v.*cos(RLON)).*dm;
      h_field(3,:,:,:,:) =  Re*(u.*cos(RLAT)).*dm;

      I_field(1,:,:,:,:) = -Re^2*(cos(RLAT).*sin(RLAT).*cos(RLON)).*dm;
      I_field(2,:,:,:,:) = -Re^2*(cos(RLAT).*sin(RLAT).*sin(RLON)).*dm;
      I_field(3,:,:,:,:) =  Re^2*(cos(RLAT).*cos(RLAT)).*dm;

      % integrate (sum) vertically so that we can compute covariances and correlations
      if savecomp == 1
        ww(:,k1:k2,:,:) = squeeze(sum(I_field,3));
        mm(:,k1:k2,:,:) = squeeze(sum(h_field,3));
      end

      % sum up all boxes and compute EAFs in rad/s.
      Xw(:,k1:k2) = prefac_h*ones(1,nt).*sum(sum(sum(h_field,3),4),5);
      Xm(:,k1:k2) = prefac_I*ones(1,nt).*sum(sum(sum(I_field,3),4),5);
      kold = k2;

      disp(['Finished   ',f])
    end % variables-not-screwed-up loop
  end % file exist if loop

end	% file loop


% this is the total AM 
X = Xw+Xm;


% save output
if savecomp == 1
  save 'temp_aam_massintegral_components.mat' Xw Xm MJD ww mm lat lon
else
  save 'temp_aam_massintegral.mat' Xw Xm MJD
end

