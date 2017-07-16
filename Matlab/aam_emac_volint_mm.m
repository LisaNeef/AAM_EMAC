%---aam_emac_volint_mm.m-------------------------:
% cycle through EMAC output files and compute atmospheric angular momentum
% then store MONTHLY MEANS of: AAM, u(300hPa), and ps 
% for aam, integrate over the whole globe and save all years of the run.
%
% see notes v3 p101
% Lisa Neef, 22 Sep 2010
%
% mods:
%  18 Oct 2010: adapted slightly to handle other data resolutions.


clear all;

%---run-specific inputs 
%runid = 'ref_T63L39';		HRES = 'T63';
runid = 'ref2_T31L39';		HRES = 'T31';
%runid = 't7_T42L39';		HRES = 'T42';
%runid = 'CCMval';		HRES = 'T42';

if runid(1:3) == 'CCM', CCMval = 1; else CCMval = 0; end

%---paths and file settings/slm
d2r = pi/180.;         %  degrees to radian
months = 1:12;		

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');

datadir = ['/dsk/nathan/lisa/ex/',runid,'/uvps_pl/'];
if CCMval,datadir = ['/dsk/nathan/lisa/ex/CCMval/Output/']; end

%---load landmask data (for IB correction)
ff_lm = ['/home/ig/neef/Data/sample_landmask_',HRES,'.nc'];
slm = nc_varget(ff_lm,'slm');
lat = nc_varget(ff_lm,'lat');
lon = nc_varget(ff_lm,'lon');
rlat = lat*d2r;
rlon = lon*d2r;
nlat = length(lat);
nlon = length(lon);

%---compute how many years this run has by computing the nr of januaries
if CCMval == 1
  flist_uv = dir([datadir,'*01_uv.nc']);
  flist_aps = dir([datadir,'*01_aps.nc']);
  nyears = min([size(flist_uv,1)-1,size(flist_aps,1)-1]);
else
  flist0 = dir([datadir,'*01.01_uvps_pl.nc']);
  nyears = size(flist0,1)-1;
end

%---initialize arrays holding the mean AAM terms for the selected months and years

nm = length(months);
Xm_mm = zeros(3,nm, nyears);
Xw_mm = zeros(3,nm, nyears);


%---cycle through netcdf files by months and years and perform AAM integrals
for im = 1:nm

  m = months(im);
  if length(num2str(m)) < 2, month = ['0',num2str(m)]; else month=num2str(m); end
  % list of files corresponding to selected month
  if CCMval == 1 
    flist = dir([datadir,'scout02____dm_*',month,'_uv.nc']);
    flist_aps = dir([datadir,'O3_*',month,'_aps.nc']);
  else
    flist = dir([datadir,'*',month,'.01_uvps_pl.nc']);
    flist_aps = flist;
  end

  %--set up temporary arrays and time axis
  f0 = [datadir,flist(2).name];
  dum = nc_varget(f0,'var131');
  [nt,nlev,nlat,nlon] = size(dum);

  %--create matricies of lats and lons (only on 1st month)
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

  % resize RLAT to suit the right month
  if nt == 29, nt = 28;	end		% for convenience, ignore Feb 29ths.
  RLAT2 = RLAT(1:nt,:,:,:);
  RLATb2 = RLATb(1:nt,:,:);
  RLON2 = RLON(1:nt,:,:,:);
  RLONb2 = RLONb(1:nt,:,:);

  %--allocate arrays for AAM terms, all years for this month 
  ww = zeros(3,nt,nyears)+NaN;
  mm = zeros(3,nt,nyears)+NaN;

  %--cycle available files and compute / collect terms.
  for iyear = 1:nyears
    if CCMval == 1
      % (note that the list of uv files starts later - we ignore 1958/59
      f = [datadir,flist(iyear+2).name];
      f_aps = [f(1:31),'O3_',f(46:52),'aps.nc'];
    else  
      f = [datadir,flist(iyear).name];
      f_aps = f;
    end
    u = nc_varget(f,'var131');

    % if the time-dimension isn't screwed up here, and all files exist, do the integral
    if size(size(u),2)==4 & exist(f_aps)==2
      v = nc_varget(f,'var132');
      time = nc_varget(f,'time');
      lev = nc_varget(f,'lev');
      nt = length(time);
      if nt == 29, nt = 28;	end		% for convenience, ignore Feb 29ths.
      u = u(1:nt,:,:,:);
      v = v(1:nt,:,:,:);

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
      LM = double(squeeze(slm(1,:,:)));			% land mask
      SM = double(LM*0);            
      SM(find(LM == 0)) = 1;			% sea mask
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

      % compute the contributions for each gridbox to the integrand.

      w_field = zeros(3,nt,nlev,nlat,nlon);    
      m_field = zeros(3,nt,nlat,nlon);    

      m_field(1,:,:,:,:) = ps_ib.*sin(RLATb2).*cos(RLATb2).*cos(RLATb2).*cos(RLONb2);
      m_field(2,:,:,:,:) = ps_ib.*sin(RLATb2).*cos(RLATb2).*cos(RLATb2).*sin(RLONb2);
      m_field(3,:,:,:,:) = ps_ib.*cos(RLATb2).*cos(RLATb2).*cos(RLATb2);

      w_field(1,:,:,:,:) = u.*sin(RLAT2).*cos(RLAT2).*cos(RLON2)-v.*cos(RLAT2).*sin(RLON2);
      w_field(2,:,:,:,:) = u.*sin(RLAT2).*cos(RLAT2).*sin(RLON2)+v.*cos(RLAT2).*cos(RLON2);
      w_field(3,:,:,:,:) = u.*cos(RLAT2).*cos(RLAT2);

      % integrate the wind term vertically
      int_w_dp = -trapz(lev,w_field,3);

      % integrate both the wind and mass terms horizontally.
      int_m_dlat = -trapz(rlat,m_field,3);
      int_m_dlat_dlon = trapz(rlon,int_m_dlat,4);
      int_w_dp_dlat = -trapz(rlat,int_w_dp,4);
      int_w_dp_dlat_dlon = trapz(rlon,int_w_dp_dlat,5);


      % store the contributions to the integral, for that month, year by year
      ww(:,:,iyear) = squeeze(int_w_dp_dlat_dlon);
      mm(:,:,iyear) = squeeze(int_m_dlat_dlon);
      disp(['Finished month',month,',  Year ', num2str(iyear)])
    end % variables-not-screwed-up loop

  end	% year loop

  %---compute means for this month
  Xm_mm(:,im,:) = squeeze(nanmean(mm,2));
  Xw_mm(:,im,:) = squeeze(nanmean(ww,2));

end % month loop


save 'aam_vol_mm_RUNID.mat' months Xm_mm Xw_mm lat lon RLAT RLON
