%---aam_emac_components_eqA.m-------------------------:
% cycle through EMAC output files and compute 
% the individual AAM terms
% output sould have units rad/s (X1 and X2) and s (X3, i.e.  dLOD)
% Lisa Neef	started 27 Apr 2010
% version _eqA: compute the individual AM terms for gridboxes of the same area.
%  see notes vol. 2 p. 85

%---choose start and stop times
clear all;

y0 = 1995;
yf = 2000;
weight = 1;
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

%---define some constants

Re = 6.371e6;    	% radius of earth in meters
Q = 7.292115e-5; 	% rotation rate of the earth in rad/s
g = 9.81;        	% grav constant in m/s2
Cm = 7.1236e37;  	% principal moment of inertia of the mantle in kgm2
CminusA = 2.34e35;       % C minus eq. MOI, same units 
d2r = pi/180.;         %  degrees to radian

%---compute AAM prefactors (see Barnes et al., 1983, p. 48)
pw12 = -(1.43*Re^3)/(CminusA*Q*g); 
pw3  = (0.998*Re^3)/(Cm*Q*g); 
pm12 = -(1.00*Re^4)/(CminusA*g); 
pm3 = (0.753*Re^4)/(Cm*g); 
prefac_m = ([pm12; pm12; pm3]);
prefac_w = ([pw12; pw12; pw3]);

%---load landmask data (for IB correction) 

dir_lm = '/home/ig/neef/Data/';
ff_lm = [dir_lm 'sample_landmask.nc'];
[data name] = r_netcdf(dir_lm,ff_lm);
n = size(data,2);
slm = data{4};
lon = data{1};
lat = data{2};

% make a lon array from -180 to 180, for plotting
dum = find(lon > 180);
lon2 = lon;
lon2(dum) = lon(dum)-360;

% define the latitude bands on which we compute AM.

lb1 = zeros(23,2);
lb1(1,:) = [73.0, 90.0]; 
lb1(2,:) = [65.9, 73.0];
lb1(3,:) = [60.4, 65.9];
lb1(4,:) = [55.7, 60.4];
lb1(5,:) = [51.5, 55.7];
lb1(6,:) = [47.7, 51.5];
lb1(7,:) = [44.1, 47.7];
lb1(8,:) = [40.7, 44.1];
lb1(9,:) = [37.5, 40.7];
lb1(10,:) = [34.4, 37.5];
lb1(11,:) = [31.4, 37.5];
lb1(12,:) = [28.6, 31.4];
lb1(13,:) = [25.8, 28.6];
lb1(14,:) = [23.0, 25.8];
lb1(15,:) = [20.4, 23.0];
lb1(16,:) = [17.7, 20.4];
lb1(17,:) = [15.1, 17.7];
lb1(18,:) = [12.6, 15.1];
lb1(19,:) = [10.0, 12.6];
lb1(20,:) = [7.5, 10.0 ];
lb1(21,:) = [5.0, 7.5];
lb1(22,:) = [2.5, 5.0];
lb1(23,:) = [0.0, 2.5];

lb2 = sort(lb1,1);
lb3 = -[lb1(:,2), lb1(:,1)];
lband = [lb3;lb2];
nb = 46;

% generate lat/lon grid and convert to radians.

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

% compute latitude and longitude widths - note that dlat depends on the latbands
dlon = rlon(2)-rlon(1);
dlat = lband(:,2) - lband(:,1);

%---initialize arrays for AM and LOD

Xw = zeros(3,ndays)+NaN;
Xm = zeros(3,ndays)+NaN;
xxw = zeros(3,nlon,nb,ndays)+NaN;
xxm = zeros(3,nlon,nb,ndays)+NaN;
MJD = zeros(1,ndays)+NaN;

%---cycle through netcdf files by year and month and compute local AAM terms
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
      MJD(k1:k2) = data{4}+refday_uv;
      u = data{5};
      v = data{6};
      lev = data{3};

      % load surface pressure data
      % here variables 3 and 4 are time and surface pressure
      [data name] = r_netcdf(dir,ff_ps);
      nt_aps = length(data{3});
        MJD(k1:k2) = data{3}+refday_aps;	% the reference date is different for aps data
      ps = data{4};

      % for surface pressure, apply Inverse barometer (IB) correction.
      LM = squeeze(slm(:,:,1));
      SM = zeros(nlon,nlat);		% craft a sea mask.
      SM(find(LM == 0)) = 1;
      dxyp = gridcellarea(nlat,nlon);                                                                            
      dxyp2 = ones(nlon,1)*dxyp';
      sea_area = sum(sum(dxyp2.*SM));
      for ii=1:nt
        pstemp = ps(:,:,ii);
        ps_sea_ave = sum(sum(pstemp.*SM.*dxyp2))/sea_area;
        ps_ave = pstemp.*LM+ps_sea_ave*SM;
        ps(:,:,ii) = ps_ave;
      end
    
      % ------compute AEF terms 
      RLAT = RLAT31(:,:,:,1:nt);
      RLON = RLON31(:,:,:,1:nt);
      RLATb = squeeze(RLAT(:,:,1,1:nt_aps));
      RLONb = squeeze(RLON(:,:,1,1:nt_aps));
      fw = zeros(3,nlon,nlat,np,nt);
      fm = zeros(3,nlon,nlat,nt_aps);
      fw(1,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*cos(RLON)-v.*cos(RLAT).*sin(RLON);
      fw(2,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*sin(RLON)+v.*cos(RLAT).*cos(RLON);
      fw(3,:,:,:,:) = u.*cos(RLAT).*cos(RLAT);
      fm(1,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*cos(RLONb);
      fm(2,:,:,:) = ps.*sin(RLATb).*cos(RLATb).*cos(RLATb).*sin(RLONb);
      fm(3,:,:,:) = ps.*cos(RLATb).*cos(RLATb).*cos(RLATb);

      % loop over latitude bands: in each band, average over that band
      fwb = zeros(3,nlon,nb,np,nt);  
      fmb = zeros(3,nlon,nb,nt);  
      for ib = 1:nb
        band = find(lat > lband(ib,1) & lat < lband(ib,2));
        fwband = fw(:,:,band,:,:);
        fmband = fm(:,:,band,:);

        fwb(:,:,ib,:,:) = squeeze(mean(fwband,3));
        fmb(:,:,ib,:) = squeeze(mean(fmband,3));
        % integrate over latitude
      %  fwb(:,:,ib,:,:) = squeeze(trapz(lat(band),fwband,3));
      %  fmb(:,:,ib,:) = squeeze(trapz(lat(band),fmband,3));
      end


      % motion term: integrate over pressure levels
      fwb_p = squeeze(trapz(lev,fwb,4));
      
      % multiply EAF integral components for each new box by the prefactor and the dlon, dlat terms
      for ib = 1:nb
        xxw(1,:,ib,k1:k2) = dlon*dlat(ib)*pw12*fwb_p(1,:,ib,:);
        xxw(2,:,ib,k1:k2) = dlon*dlat(ib)*pw12*fwb_p(2,:,ib,:);
        xxw(3,:,ib,k1:k2) = dlon*dlat(ib)*pw3*fwb_p(3,:,ib,:);
        xxm(1,:,ib,k1:k2) = dlon*dlat(ib)*pm12*fmb(1,:,ib,:);
        xxm(2,:,ib,k1:k2) = dlon*dlat(ib)*pm12*fmb(2,:,ib,:);
        xxm(3,:,ib,k1:k2) = dlon*dlat(ib)*pm3*fmb(3,:,ib,:);
      end

      kold = k2;


      disp(['Finished   ',num2str(y),'-',mstr])

    end % file exist if loop

  end	% month loop
end	% year loop




% save the big state arrays
  save 'spatialAM_1995_2000_eqA.mat' xxw xxm y0 yf  MJD lon lat lband
