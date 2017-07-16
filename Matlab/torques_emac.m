%---aam_emac.m-------------------------:
% cycle through EMAC output files and compute atmospheric torques
% Lisa Neef	started 17 May 2010

%---choose start and stop times

y0 = 1960;
yf = 2000;

%---paths and file settings

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');

dir1 = '/dsk/nathan/lisa/Data/MESSy_input/EVAL/echam5_ini/add_spec/T42L90MA/';
dir2 = '/dsk/nathan/lisa/Data/ex/T42L90MA/data/';
ff1 = 'T42_19960101_surf.nc';
ff2 = 'T42L39____195801.01.nc';
pref = 'T42L39____';
suff  = '.01.nc';

%---start and stop times and time axis
mm = 1:1:12;

day0 = date2mjd(y0,1,1,0,0,0);
dayf = date2mjd(yf,12,31,0,0,0);
ndays = dayf - day0;
refday = date2mjd(1958,1,1,0,0,0);

%---load landmask data, for IB correction and to generate lat/lon grid (and convert to radians) 
dir_lm = '/home/ig/neef/Data/';
ff_lm = [dir_lm 'sample_landmask.nc'];
[data name] = r_netcdf(dir_lm,ff_lm);
n = size(data,2);
slm = data{4};
lon = data{1};
lat = data{2};
d2r = pi/180;
rlat = lat*d2r;
rlon = lon*d2r;
nlat = length(lat);
nlon = length(lon);
np = 27;
RLAT31 = zeros(nlon,nlat,np,31);
RLON31 = zeros(nlon,nlat,np,31);

for ilev=1:np
  for it = 1:31
    RLAT31(:,:,ilev,it) = ones(128,1)*rlat';
    RLON31(:,:,ilev,it) = rlon*ones(1,64);
  end
end

RLAT30 = RLAT31(:,:,:,1:30);
RLON30 = RLON31(:,:,:,1:30);


%---initialize arrays for 3 torques: TF, TM, and TG (foll. Lejenas et al. 1999)

TG = zeros(1,ndays)+NaN;
TF = zeros(1,ndays)+NaN;
TM = zeros(1,ndays)+NaN;

%---cycle through netcdf files by year and month and perform AAM integrals
kold = 0;

for y = y0:1:yf
  for m = 1:12 
    mstr = num2str(m);
    if size(mstr,2)==1, mstr=['0',mstr]; end
    datestr = [num2str(y),mstr];

    ff = [pref,num2str(y),mstr,suff];

    if exist([dir2,ff]) == 2 

      % load model output:
      [data name] = r_netcdf(dir2,ff);
      n = size(data,2);
      for i = 1:n
        eval([name{i},'=data{i};',]);
      end
      nt = length(time);
      k1 = kold+1;
      k2 = k1+nt-1;
      MJD(k1:k2) = time+refday;

      % (1) surface wind stress

      % (2) surface GW stress

      % (3) surface pressure

      % (4) orographic height

  
      % to integrate spherically, modulate u and v by sines and cosines.
      fw = zeros(3,nlon,nlat,np,nt);
      fm = zeros(3,nlon,nlat,nt_aps);
      RLAT = RLAT31(:,:,:,1:nt);
      RLON = RLON31(:,:,:,1:nt);
      RLATb = squeeze(RLAT(:,:,1,1:nt_aps));
      RLONb = squeeze(RLON(:,:,1,1:nt_aps));
      fw(1,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*cos(RLON)-v.*cos(RLAT).*sin(RLON);
      fw(2,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*sin(RLON)+v.*cos(RLAT).*cos(RLON);
      fw(3,:,:,:,:) = u.*cos(RLAT).*cos(RLAT);

      fm(1,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*cos(RLONb);
      fm(2,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*sin(RLONb);
      fm(3,:,:,:) = ps.*cos(RLATb).*cos(RLATb).*cos(RLATb);

      % motion term integral:
      fw_lon = trapz(rlon,fw,2);
      fw_lonlat = trapz(rlat,fw_lon,3);
      fw_lonlatp = trapz(lev,fw_lonlat,4);
    
      % mass term integral:
      fm_lon = trapz(rlon,fm,2);
      fm_lonlat = trapz(rlat,fm_lon,3);

      % compute EAFs in rad/s.
      Xw(:,k1:k2) = prefac_w*ones(1,nt).*squeeze(fw_lonlatp);
      Xm(:,k1:k2) = prefac_m*ones(1,nt_aps).*squeeze(fm_lonlat);
      kold = k2;


      disp(['Finished   ',num2str(y),'-',mstr])

    end % file exist if loop

  end	% month loop
end	% year loop


% take out bad flags
Xtemp = Xw+Xm;
Xwtemp = Xw;
Xmtemp = Xm;
mjdtemp = MJD;

b1 = find(abs(Xwtemp(1,:)) > 1e20);
b2 = find(abs(Xwtemp(2,:)) > 1e20);
b3 = find(abs(Xwtemp(3,:)) > 1e20);
bad = [b1, b2, b3];

for ii=1:3
  Xwtemp(ii,bad) = NaN;
  Xmtemp(ii,bad) = NaN;
end

Xw = Xwtemp;
Xm = Xmtemp;


% this is the total AM 
X = Xw+Xm;




save aam_emac_1960_2000_withIB_set3.mat
