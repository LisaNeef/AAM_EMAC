%---aam_emac_volint_mm.m-------------------------:
% cycle through EMAC output files and compute atmospheric angular momentum
% then store MONTHLY MEANS of: AAM, u(300hPa), and ps for months 2,4,6,8,10,12
%
% see notes v3 p70
% Lisa Neef, 24 Aug 2010
%
% mods:
%    27 Sep 2010 - modify to accomodate CCMval data

clear all;

%---run-specific inputs 
%runid = 't7_T42L39';
%CCMval = 0;
CCMval = 1;                    % select if this is the data we want to use
ref_date_str = '19580101';      % date of start of simulation.

%---select months to survey
%months = [8,12];
months = 2:2:12;		
%months = 1:12;		
%months = 2;	

%---paths and file settings

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');

if CCMval == 1
  datadir = ['/dsk/nathan/lisa/CCMval/Output/'];
else
  datadir = ['/dsk/nathan/lisa/ex/',runid,'/uvps_pl/'];
end



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
slm_slice = squeeze(slm(1,:,:));
%LM = flipud(slm_slice);			% land mask
lat_sm = nc_varget(ff_lm,'lat');
lon_sm = nc_varget(ff_lm,'lon');


%---load the first file to get the grid (could be different from seamask - check for this!)
  % list of files corresponding to selected month
  if CCMval == 1
    flist0 = dir([datadir,'*_uv.nc']);
  else
    flist0 = dir([datadir,'*',month,'.01_uvps_ml.nc']);
  end
lat = nc_varget([datadir,flist0(1).name],'lat');
lon = nc_varget([datadir,flist0(1).name],'lon');
rlat = lat*d2r;
rlon = lon*d2r;
nlat = length(lat);
nlon = length(lon);

yr = str2num(ref_date_str(1:4));
mr = str2num(ref_date_str(5:6));
dr = str2num(ref_date_str(7:8));
ref_day = date2mjd(yr,mr,dr,0,0,0);		% initial date in avail. files

%---initialize arrays holding the mean fields for the selected months

nm = length(months);
u300_mean = zeros(nm, nlat, nlon);
ps_mean = zeros(nm, nlat, nlon);
Xm_mean = zeros(3,nm, nlat, nlon);
Xw_mean = zeros(3,nm, nlat, nlon);

%---cycle through netcdf files by month and perform AAM integrals
for im = 1:nm

  m = months(im);
  if length(num2str(m)) < 2, month = ['0',num2str(m)]; else month=num2str(m); end
  % list of files corresponding to selected month
  if CCMval == 1
    flist = dir([datadir,'scout02____dm_*',month,'_uv.nc']);
    flist_aps = dir([datadir,'O3_*',month,'_aps.nc']);
    nf = min([size(flist,1)-1,size(flist_aps,1)-1]);
  else
    flist = dir([datadir,'*',month,'.01_uvps_ml.nc']);
    flist_aps = flist;
    nf = size(flist,1)-1;
  end

  %--set up temporary arrays and time axis
  f0 = [datadir,flist(1).name];
  dum = nc_varget(f0,'var131');
  [nt,nlev,nlat,nlon] = size(dum);
  if m == 2, nt = 29; end		% allocate enough space for leap years


  %--create matricies of lats and lons (only on 1st month)
  if im == 1
    RLAT = zeros(31,nlev,nlat,nlon);
    RLON = RLAT;
    for ilev=1:nlev
      for it = 1:31
        RLAT(it,ilev,:,:) = rlat*ones(1,128);
        RLON(it,ilev,:,:) = ones(64,1)*rlon';
      end
    end
    RLATb = squeeze(RLAT(:,1,:,:));
    RLONb = squeeze(RLON(:,1,:,:));
  end


  %--allocate arrays for AAM terms, u(300hPa), and ps, covering all dates in all files for this month.

  pp  = zeros(nt*nf,nlat,nlon)+NaN;
  uu  = zeros(nt*nf,nlat,nlon)+NaN;
  ww = zeros(3,nt*nf,nlat,nlon)+NaN;
  mm = zeros(3,nt*nf,nlat,nlon)+NaN;

  %--cycle available files and compute / collect terms.
  k0 = 1;

  for ifile = 1:nf
    if CCMval == 1
      % (note that the list of uv files starts later - we ignore 1958/59
      f = [datadir,flist(ifile+2).name];
      f_aps = [f(1:31),'O3_',f(46:52),'aps.nc'];
    else
      f = [datadir,flist(ifile).name];
      f_aps = f;
    end
    u = nc_varget(f,'var131');

    % if the time-dimension isn't screwed up here, do the integral
    if size(size(u),2)==4 & exist(f_aps)==2
      v = nc_varget(f,'var132');
      %lsp = nc_varget(f,'lsp');
      time = nc_varget(f,'time');
      lev = nc_varget(f,'lev');
      nt = length(time);


     % remove bad points.
      no_data = find(time == 0);
      if length(no_data) > 0
        disp(['Found a few pad points in file  ',f])
        u(no_data) = NaN;
        v(no_data) = NaN;
        time(no_data) = NaN;
      end

      %  load surface pressure, and convert ln surface pressure to surface pressure
      %  note that in CCMval data, the latitude dimension is reversed in the O3 files.
      %    (so it needs to be flipped)
      if CCMval == 1
        ps_dum = nc_varget(f_aps,'APS');
        ps = flipdim(ps_dum,2);
      else
        lsp = nc_varget(f,'lsp');
        ps = exp(lsp);
      end

      % inverse barometer correction
      LM = slm_slice;
      SM = double(LM*0);
      SM(find(LM == 0)) = 1;                    % sea mask
      dxyp = gridcellarea(nlat,nlon);
      dxyp2 = dxyp*ones(1,nlon);
      sea_area = sum(sum(dxyp2.*SM));
      ps_ib = ps(1:nt,:,:);
      for ii=1:nt
        pstemp = double(squeeze(ps(ii,:,:)));
        ps_sea_ave = sum(sum(pstemp.*SM.*dxyp2))/sea_area;
        ps_ave = pstemp.*LM+ps_sea_ave*SM;
        ps_ib(ii,:,:) = ps_ave;
      end

      % resize RLAT to suit the right month
      RLAT2 = RLAT(1:nt,:,:,:);
      RLATb2 = RLATb(1:nt,:,:);
      RLON2 = RLON(1:nt,:,:,:);
      RLONb2 = RLONb(1:nt,:,:);


      % collect pressure and 300 hPa zonal winds (not p is in Pa)
      if im == 1, p300 = find(lev == 30000); end
      uu(k0:k0+nt-1,:,:) = squeeze(u(:,p300,:,:));
      pp(k0:k0+nt-1,:,:) = ps_ib;

      % compute the contributions for each gridbox to the integrand.

      w_field = zeros(3,nt,nlev,nlat,nlon);    
      m_field = zeros(3,nt,nlat,nlon);    

      m_field(1,:,:,:,:) = ps_ib.*sin(RLATb2).*cos(RLATb2).*cos(RLATb2).*cos(RLONb2);
      m_field(2,:,:,:,:) = ps_ib.*sin(RLATb2).*cos(RLATb2).*cos(RLATb2).*sin(RLONb2);
      m_field(3,:,:,:,:) = ps_ib.*cos(RLATb2).*cos(RLATb2).*cos(RLATb2);

      w_field(1,:,:,:,:) = u.*sin(RLAT2).*cos(RLAT2).*cos(RLON2)-v.*cos(RLAT2).*sin(RLON2);
      w_field(2,:,:,:,:) = u.*sin(RLAT2).*cos(RLAT2).*sin(RLON2)+v.*cos(RLAT2).*cos(RLON2);
      w_field(3,:,:,:,:) = u.*cos(RLAT2).*cos(RLAT2);

      % uncomment the below to get a snapshot of what we're looking at in this month
      % save 'temp.mat' w_field u m_field ps lat lon ps_ib ps_dum
      % disp('saved a file')

      % integrate the wind term vertically
      % note that the limits of integration for lev are reversed from how the integral
      % is defined in Barnes et al (1983) - therefore add (-) to integral.
      % (see notes vol. 3, p. 41)
      int_w_dp = -trapz(lev,w_field,3);

      % store the contributions to the integral
      ww(:,k0:k0+nt-1,:,:) = squeeze(int_w_dp);
      mm(:,k0:k0+nt-1,:,:) = m_field;

      k0 = k0+nt;

      disp(['Finished   ',f])
    end % variables-not-screwed-up loop

  end	% file loop

  % truncate the arrays to remove the empty spaces at the end
  dum = pp(:,1,1);
  cutoff = min(find(isfinite(dum) == 0))-1;
  if isfinite(cutoff) == 1
    u300 = uu(1:cutoff,:,:);
    ps   = pp(1:cutoff,:,:);
    Xm   = mm(:,1:cutoff,:,:);
    Xw   = ww(:,1:cutoff,:,:);
  else
    u300 = uu;
    ps   = pp;
    Xm   = mm;
    Xw   = ww;
  end


  %---compute means for this month
  u300_mean(im,:,:) = squeeze(mean(u300,1));
  ps_mean(im,:,:) = squeeze(mean(ps,1));
  Xm_mean(:,im,:,:) = squeeze(mean(Xm,2));
  Xw_mean(:,im,:,:) = squeeze(mean(Xw,2));

end % month loop


save 'temp_aamcomponents_mm.mat' months u300_mean ps_mean Xm_mean Xw_mean lat lon RLAT RLON
