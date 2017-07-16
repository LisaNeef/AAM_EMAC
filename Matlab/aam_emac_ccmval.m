function aam_emac_ccmval(dec)
%---aam_emac.m-------------------------:
% cycle through EMAC output files and compute atmospheric angular momentum
% Lisa Neef	started 16 Feb 2010
% see notes vol. 2 p. 22
% MODS:
%   19 mar 2010 - fixed my IB correction (notes vol 2 p 52)
%   12 aug 2010 - turn into a function, clean up code and make more like aam_emac_volintegral
%	here also fixed a sign mistake in the mass term!
%
% INPUTS:
%   dec: number given the decade for which we want to do the computation (faster than doing entire run)


%---paths and file settings
decade = num2str(dec);

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');

datadir = '/dsk/nathan/lisa/CCMval/Output/';
list_uv = dir([datadir,'*',decade(1:3),'*_uv.nc']);
list_ps = dir([datadir,'*',decade(1:3),'*_aps.nc']);
nf = size(list_uv,1);

%---define some constants

Re = 6.371e6;    	% radius of earth in meters
Q = 7.292115e-5; 	% rotation rate of the earth in rad/s
g = 9.81;        	% grav constant in m/s2
C = 7.1236e37;  	% principal moment of inertia of the mantle in kgm2
CminusA = 2.34e35;       % C minus eq. MOI, same units 
d2r = pi/180.;         %  degrees to radian

%---compute AAM prefactors (see Barnes et al., 1983, p. 48)
pm12 = -Re^4/(g*CminusA);
pw12 = -1.43*Re^3/(Q*g*CminusA);
pm3  = 0.7*Re^4/(g*C);
pw3  = Re^3/(Q*g*C);
prefac_w = [pw12; pw12; pw3];
prefac_m = [pm12; pm12; pm3];

%---load landmask data (for IB correction) 
ff_lm = '/home/ig/neef/Data/sample_landmask.nc';
slm = nc_varget(ff_lm,'slm');
lat = nc_varget(ff_lm,'lat');
lon = nc_varget(ff_lm,'lon');
rlat = lat*d2r;
rlon = lon*d2r;

%----these are the reference dates from which each file starts counting
refday_uv = date2mjd(1958,1,1,0,0,0);
refday_aps = date2mjd(1901,1,15,0,0,0);

%---set up grids and time axis
f0 = [datadir,list_uv(1).name];
ff = [datadir,list_uv(nf).name];

dum = nc_varget(ff,'var131');
[ntf,nlev,nlat,nlon] = size(dum);

% initial date of run
time0 = nc_varget(f0,'time');
day0 = time0(1)+refday_uv;

% final date of the times considered here.
timef = nc_varget(ff,'time');
dayf = timef(ntf)+refday_uv;
ndays = round(dayf - day0)+1;				% number of days in avail. files

% matricies of lats and lons that are the same grid as the u and v fields
% (make this array bigger than 31 days to leave from for any NaNs in the output)
RLAT31 = zeros(33,nlev,nlat,nlon);
RLON31 = RLAT31;
for ilev=1:nlev
  for it = 1:33
    RLAT31(it,ilev,:,:) = rlat*ones(1,128);
    RLON31(it,ilev,:,:) = ones(64,1)*rlon';
  end
end

%---initialize arrays for AM and LOD

Xw = zeros(3,ndays)+NaN;
Xm = zeros(3,ndays)+NaN;
MJD = zeros(1,ndays)+NaN;

%---cycle through netcdf files by year and month and perform AAM integrals
kold = 0;

for ifile = 1:nf
  f_uv = [datadir,list_uv(ifile).name];
  date_month = list_uv(ifile).name(15:20);
  f_ps = [datadir,'O3_',date_month,'_aps.nc'];

  if exist(f_uv) + exist(f_ps) == 4 

    % load wind data
    u = nc_varget(f_uv,'var131');
    v = nc_varget(f_uv,'var132');
    time = nc_varget(f_uv,'time');
    lev = nc_varget(f_uv,'lev');
    nt = length(time);

    % load surface pressure data
    % note that in these files the latitude dimension is reversed from the *uv.nc files:
    % deal with it by reversing that dimension of the ps array.
    % (see notes vol. 3 p. 48)
    ps_dum = nc_varget(f_ps,'APS');
    ps = flipdim(ps_dum,2);

    % note down the dates
    k1 = kold+1;
    k2 = k1+nt-1;
    MJD(k1:k2) = time+refday_uv;

    % generate approproiate lat/long arrays
    RLAT = RLAT31(1:nt,:,:,:);
    RLON = RLON31(1:nt,:,:,:);
    RLATb = squeeze(RLAT31(1:nt,1,:,:));
    RLONb = squeeze(RLON31(1:nt,1,:,:));

    %*******correct this later*****************
    %  % for surface pressure, apply Inverse barometer (IB) correction.
    %  LM = squeeze(slm(:,:,1));
    %  SM = zeros(nlon,nlat);            % craft a sea mask.
    %  SM(find(LM == 0)) = 1;
    %  dxyp = gridcellarea(nlat,nlon);
    %  dxyp2 = ones(nlon,1)*dxyp';
    %  sea_area = sum(sum(dxyp2.*SM));
    %  for ii=1:nt
    %    pstemp = ps(:,:,ii);
    %    ps_sea_ave = sum(sum(pstemp.*SM.*dxyp2))/sea_area;
    %    ps_ave = pstemp.*LM+ps_sea_ave*SM;
    %    ps(:,:,ii) = ps_ave;
    %  end


    % compute the contributions for each gridbox to the integrand.

    w_field = zeros(3,nt,nlev,nlat,nlon);    
    m_field = zeros(3,nt,nlat,nlon);    

    m_field(1,:,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*cos(RLONb);
    m_field(2,:,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*sin(RLONb);
    m_field(3,:,:,:,:) = ps.*cos(RLATb).*cos(RLATb).*cos(RLATb);

    w_field(1,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*cos(RLON)-v.*cos(RLAT).*sin(RLON);
    w_field(2,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*sin(RLON)+v.*cos(RLAT).*cos(RLON);
    w_field(3,:,:,:,:) = u.*cos(RLAT).*cos(RLAT);

    % integrate over each dimension to get the global value
    % note that the limits of integration for lat and lev are reversed from how the integral
    % is defined in Barnes et al (1983) - therefore add (-) before these trapz integrals.
    % (see notes vol. 3, p. 41)
    int_w_dp = -trapz(lev,w_field,3);
    int_w_dpdlat = -trapz(rlat,int_w_dp,4);
    int_w_dpdlatdlon = trapz(rlon,int_w_dpdlat,5);

    int_m_dlat = -trapz(rlat,m_field,3);
    int_m_dlatdlon = trapz(rlon,int_m_dlat,4);

    % multiply by prefactors to compute EAFs in rad/s.
    Xw(:,k1:k2) = prefac_w*ones(1,nt).*int_w_dpdlatdlon;
    Xm(:,k1:k2) = prefac_m*ones(1,nt).*int_m_dlatdlon;
    kold = k2;

    disp(['Finished   ',f_uv])

  end % file exist if loop
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




save 'temp_aam_volintegral_ccmval.mat' Xw Xm MJD
