%---aam_emac_volint_regional_mm.m-------------------------:
% cycle through EMAC output files and compute atmospheric angular momentum
% then store MONTHLY MEANS of: AAM, u(300hPa), and ps for months 2,4,6,8,10,12
%
% see notes v3 p70
% Lisa Neef, 24 Aug 2010
%
% mods:
%    27 Sep 2010 - modify to accomodate CCMval data
%    21 Oct 2010 - instead of defining constats here, load from common function.
%		 - customize so that other other resolutions also work
%     9 Nov 2010 - add possibility of procesing ERAinterim data, clean up code.
%    17 Nov 2010 - add option of processing only last 20 years, to make it comparable to ERA interim.
%    17 Jan 2011 - option to not save the average, but individual years (see notes vol. 4, p. 69)

clear all;

%---run-specific inputs 
runid = 'ERAint';               HRES = 'T42';		year0 = 1989;
%runid = 'ERA40';               HRES = 'T42';		year0 = 1989;
%runid = 'CCMval';               HRES = 'T42';		year0 = 1958;
%runid = 'ref2_T31L39';          HRES = 'T31';		year0 = 1958;
%runid = 't7_T42L39';            HRES = 'T42';		year0 = 1958;
%runid = 'ref_T63L39';            HRES = 'T63';		year0 = 1958;

ref_date_str = '19580101';      % date of start of simulation.

if runid(1:3) == 'CCM', CCMval = 1; else CCMval = 0; end	
if runid(1:3) == 'ERA', ERAint = 1; else ERAint = 0; end	

%---define months to survey
months = 1:12;		

%---paths and file settings

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');

datadir = ['/dsk/nathan/lisa/ex/',runid,'/uvps_pl/'];
if CCMval; datadir = ['/dsk/nathan/lisa/ex/CCMval/Output/']; end
if ERAint; datadir = ['/dsk/nathan/lisa/ERAinterim/nc/']; end

%  load constants and compute AAM prefactors
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)

% prefactors (Gross 09)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));
pm3  = 0.997*0.750*Re^4/(g*Cm);
pw3  = 0.997*Re^3/(Q*g*Cm);
pm = [pm12; pm12; pm3];
pw = [pw12; pw12; pw3];


%---load landmask data (for IB correction)
ff_lm = ['/home/ig/neef/Data/sample_landmask_',num2str(HRES),'.nc'];
slm = nc_varget(ff_lm,'slm');
slm_slice = squeeze(slm(1,:,:));
lat_sm = nc_varget(ff_lm,'lat');
lon_sm = nc_varget(ff_lm,'lon');


%---load the first file to get the grid (could be different from seamask - check for this!)
flist0 = dir([datadir,'*.01_uvps_pl.nc']);
if size(flist0,1) == 0, ml_flag = 1; else ml_flag = 0; end
if ml_flag, flist0 = dir([datadir,'*.01_uvps_ml.nc']); end  
if CCMval,  flist0 = dir([datadir,'*_uv.nc']); end
if ERAint, flist0 = dir([datadir,'*131*.nc']); end

lat = nc_varget([datadir,flist0(1).name],'lat');
lon = nc_varget([datadir,flist0(1).name],'lon');
rlat = lat*d2r;
rlon = lon*d2r;
nlat = length(lat);
nlon = length(lon);

%---find out how many years are available for this run 
  dum = char(flist0.name);
  nfiles = size(dum,1);
  ny = floor(nfiles/12);

%---initialize arrays holding the mean fields for the selected months

nm = length(months);
u300_mean = zeros(nm, ny, nlat, nlon)+NaN;
ps_mean = zeros(nm, ny, nlat, nlon)+NaN;
Xm_mean = zeros(3,nm, ny, nlat, nlon)+NaN;
Xw_mean = zeros(3,nm, ny, nlat, nlon)+NaN;

%---cycle through netcdf files by month and perform AAM integrals
for im = 1:nm

  m = months(im);
  if length(num2str(m)) < 2, month = ['0',num2str(m)]; else month=num2str(m); end

  %--set up temporary arrays and time axis appropriate for this month.
  flist_u = dir([datadir,'*',month,'.01_uvps_pl.nc']);
  if ml_flag, flist_u = dir([datadir,'*',month,'.01_uvps_ml.nc']); end
  if CCMval, flist_u = dir([datadir,'scout02____dm_*',month,'_uv.nc']); end
  if ERAint, flist_u = dir([datadir,'eia6_*131.pl_T42_',month,'.nc']); end

  f0 = [datadir,flist_u(1).name];
  dum = nc_varget(f0,'var131');
  [nt,nlev,nlat,nlon] = size(dum);
  if m == 2, nt = 29; end		% allocate enough space for leap years

  %--create matrices of lats and lons (only on 1st month)
  if im == 1
    RLAT = zeros(31,nlev,nlat,nlon);
    RLON = RLAT;
    for ilev=1:nlev
      for it = 1:31
        RLAT(it,ilev,:,:) = rlat*ones(1,nlon);
        RLON(it,ilev,:,:) = ones(nlat,1)*rlon';
      end
    end
    RLATb = squeeze(RLAT(:,1,:,:));
    RLONb = squeeze(RLON(:,1,:,:));
  end


  %--allocate arrays for AAM terms, u(300hPa), and ps, covering all dates in all files for this month.

  pp  = zeros(nt,ny,nlat,nlon)+NaN;
  uu  = zeros(nt,ny,nlat,nlon)+NaN;
  ww = zeros(3,nt,ny,nlat,nlon)+NaN;
  mm = zeros(3,nt,ny,nlat,nlon)+NaN;

  %--cycle through all years and collect monthly mean field values

  for iyear = 1:ny
    year = num2str(year0+iyear-1);
    if length(runid) >= 9, runid_short = runid(1:9); else runid_short = runid; end
    fu = [datadir,runid_short,'_',year,month,'.01_uvps_pl.nc'];
    if ml_flag, fu = [datadir,runid_short,'_',year,month,'.01_uvps_ml.nc']; end
    fv = fu;
    fp = fu;
    if CCMval 
      % (note that the list of uv files starts later - we ignore 1958/59
      fu = [datadir,'scout02____dm_',year,month,'_uv.nc'];
      fv = fu;
      fp = [fu(1:34),'O3_',fu(49:55),'aps.nc'];
    end
    if ERAint
      fu = [datadir,'eia6_',year,'.131.pl_T42_',month,'.nc'];
      fv = [datadir,'eia6_',year,'.132.pl_T42_',month,'.nc'];
      fp = [datadir,'eia6_',year,'.152_T42_',month,'.nc'];
    end

    if exist(fu) == 0, disp(['Cant find file  ' fu]), end
    if exist(fv) == 0, disp(['Cant find file  ' fv]), end
    if exist(fp) == 0, disp(['Cant find file  ' fp]), end

    % if the files all exist and time-dimension isn't screwed up here, do the integral
    exist_check = exist(fu)+exist(fv)+exist(fp);
    if exist(fu) == 2, u = nc_varget(fu,'var131'); else u = NaN; end
    if size(size(u),2)==4 & exist_check == 6
      v = nc_varget(fv,'var132');
      time = nc_varget(fu,'time');
      lev = nc_varget(fu,'lev');
      nt = length(time);
      % ignore 29th of feb if it's era
      if ERAint & m == 2
         nt = 28;
         v = v(1:nt,:,:,:);
         u = u(1:nt,:,:,:);
      end

     % remove bad points.
      no_data = find(time == 0);
      if length(no_data) > 0
        disp(['Found a few pad points in file  ',fu])
        u(no_data) = NaN;
        v(no_data) = NaN;
        time(no_data) = NaN;
      end

      %  load surface pressure, and convert ln surface pressure to surface pressure
      %  note that in CCMval data, the latitude dimension is reversed in the O3 files.
      %    (so it needs to be flipped)
      if CCMval 
        ps_dum = nc_varget(fp,'APS');
        ps = flipdim(ps_dum,2);
      end
      if ERAint
        % var152 is ln(ps)
        lsp = nc_varget(fp,'var152');
        ps = exp(lsp);
      end
      if CCMval == 0 & ERAint == 0
        lsp = nc_varget(fp,'lsp');
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
      uu(1:nt,iyear,:,:) = squeeze(u(:,p300,:,:));
      pp(1:nt,iyear,:,:) = ps_ib;

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
      % resize RLAT to suit the right month
      RLAT2 = RLAT(1:nt,:,:,:);
      RLATb2 = RLATb(1:nt,:,:);
      RLON2 = RLON(1:nt,:,:,:);
      RLONb2 = RLONb(1:nt,:,:);

      % integrate the wind term vertically
      int_w_dp = -trapz(lev,w_field,3);
      % note that pressure levels are sorted from top to bottom for ERAinterim data
      if ERAint, int_w_dp = trapz(lev,w_field,3); end

      % store the contributions to the integral
      ww(:,1:nt,iyear,:,:) = squeeze(int_w_dp);
      mm(:,1:nt,iyear,:,:) = m_field;

      disp(['Finished   ',fu])
    end % variables-not-screwed-up loop

    % within this year, compute monthly-means of all fields.
    dum = squeeze(uu(:,iyear,1,1));
    g = find(isfinite(dum));
    if length(g) > 0
      u300_mean(im,iyear,:,:) = squeeze(mean(uu(g,iyear,:,:),1));
      ps_mean(im,iyear,:,:) = squeeze(mean(pp(g,iyear,:,:),1));
      Xw_mean(:,im,iyear,:,:) = squeeze(mean(ww(:,g,iyear,:,:),2));
      Xm_mean(:,im,iyear,:,:) = squeeze(mean(mm(:,g,iyear,:,:),2));
    end

  end	% year loop

end % month loop

save 'dump_am_emac_volint_regional_mm.mat'

save 'regional_comp_mm_RUNID_allyears.mat' months u300_mean ps_mean Xm_mean Xw_mean lat lon RLAT RLON
