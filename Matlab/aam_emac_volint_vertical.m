%function aam_emac_volint_vertical(runid,dec)
%
% Cycle through EMAC output files and compute atmospheric angular momentum.
% The integral is computed over volume elements; 
% this requires output on pressure levels (i.e. from afterburner).
%
%
% This script produces timeseries of chi_1, chi_2, and chi_3 over a given decade 
% (as it typically takes really long to do an entire run).
% The unit is radians for chi_1 and chi_2, and rad/2 for chi_3.
%
% Lisa Neef
% started 20 Dec 2010
% In this version  we save vertical slices, integrated only over lat and lon.
% (See notes vol. 4, p. 65b.)
%
%  SYNTAX:  aam_emac_volintegral(runid,decade)
%  ALTERNATIVELY: put the inputs in by hand below, comment out the "function" line above, and run like a regular program.
%
%  INPUT:
%	runid = name of the run;
%	decade = year at which to start decadal loop (eg 1970);
%
%  OUTPUT:
%	Xm, Xw: the global integrated AAM exciation functions (not yet nondimensionalized)
%	ww: the vertical contributions, after integrating over lat and lon


clear all;

%--------- inputs-----[comment out to run this as a function]--------------------------

% select the name of a run and its horizontal resolution
runid = 'CCMval';		HRES = 'T42';
%runid = 'ref2_T31L39';		HRES = 'T31';
%runid = 't7_T42L39';		HRES = 'T42';		
%runid = 'ref_T63L39';		HRES = 'T63';		
%runid = 'noQBO_T31L39';		HRES = 'T31';

ref_date_str = '19580101';	% first date of simulation or dataset
dec = 1980;			% decade to process here

if runid(1:3) == 'CCM', CCMval = 1; else CCMval = 0; end	% flag for CCMVal runs

ib = 1;			% set to 1 to include IB approximation
savecomp = 1;		% select to save the constituents of the AAM integral (for covariances)

%--------------------------------------------------------------------------------------

%---paths ----------------------------------------------------------------------------------

addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');

%--------------------------------------------------------------------------------------
decade = num2str(dec);
if CCMval 
  datadir = ['/dsk/nathan/lisa/ex/CCMval/Output/'];
  flist = dir([datadir,'*',decade(1:3),'*_uv.nc']);
else
  datadir = ['/dsk/nathan/lisa/ex/',runid,'/uvps_pl/'];
  flist = dir([datadir,'*',decade(1:3),'*_uvps_pl.nc']);
  % fix for old filename that say ml (even tho it's pressure levels.)
  if length(flist) == 0, flist = dir([datadir,'*',decade(1:3),'*_uvps_ml.nc']); end
end
nf = size(flist,1);

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
lat = nc_varget(ff_lm,'lat');
lon = nc_varget(ff_lm,'lon');
rlat = lat*d2r;
rlon = lon*d2r;

% compute the reference day for this run
yr = str2num(ref_date_str(1:4));
mr = str2num(ref_date_str(5:6));
dr = str2num(ref_date_str(7:8));
ref_day = date2mjd(yr,mr,dr,0,0,0);		% initial date in avail. files

%---set up grids and time axis

f0 = [datadir,flist(1).name];
ff = [datadir,flist(nf).name];

dm = nc_varget(ff,'var131');
[ntf,nlev,nlat,nlon] = size(dm);

% specify initial date of files considered here.
time0 = nc_varget(f0,'time');
day0 = time0(1)+ref_day;

% final date of the times considered here.
timef = nc_varget(ff,'time');
dayf = timef(ntf)+ref_day;

ndays = round(dayf - day0);			% number of days in avail. files
opd = round(ntf/30);				% number of obs (or samples)  per day

% matricies of lats and lons that are the same grid as the u and v fields
% (make this array bigger than 31 days to leave from for any NaNs in the output)
RLAT31 = zeros(33*opd,nlev,nlat,nlon);
RLON31 = RLAT31;
for ilev=1:nlev
  for it = 1:33*opd
    RLAT31(it,ilev,:,:) = rlat*ones(1,nlon);
    RLON31(it,ilev,:,:) = ones(nlat,1)*rlon';
  end
end

%---initialize arrays for AM and LOD

ww = zeros(3,ndays*opd,nlev)+NaN;
Xw = zeros(3,ndays*opd)+NaN;
Xm  = zeros(3,ndays*opd)+NaN;
MJD = zeros(1,ndays*opd)+NaN;


%---cycle through netcdf files by year and month and perform AAM integrals
kold = 0;

for ifile = 1:nf
  f = [datadir,flist(ifile).name];
  f_aps = f;
  if CCMval 
    dum = explode(f,'_');
    dum_aps = dir([char(datadir),'O3_',char(dum(3)),'*aps.nc']); 
    if length(dum_aps) ~= 0
      f_aps = [datadir,dum_aps(1).name];
    else
      disp(['cant find pressure file for CCMVal',dum(3)])
      f_aps = '';
    end
  end

  % if the u-file exist, check it out.
  if exist(f) == 2  
    u = nc_varget(f,'var131');

    % if the time-dimension isn't screwed up here, and all files exist, do the integral
    if size(size(u),2)==4 & exist(f_aps)==2
      v = nc_varget(f,'var132');
      time = nc_varget(f,'time');
      lev = nc_varget(f,'lev');
      nt = length(time);

      % load surface pressyre; convert ln surface pressure to surface pressure
      if CCMval 
        ps_dum = nc_varget(f_aps,'APS');
        ps = flipdim(ps_dum,2);
      else
        lsp = nc_varget(f,'lsp');
        ps = exp(lsp);
      end

      % remove bad points.
      no_data = find(time == 0);
        if length(no_data) > 0, disp('Found a few pad points'), end
        u(no_data) = NaN;
        v(no_data) = NaN;
        mjd(no_data) = NaN;
        dm(no_data) = NaN;
        time(no_data) = NaN;
     
      % note down the dates
      k1 = kold+1;
      k2 = k1+nt-1;
      MJD(k1:k2) = time+ref_day;

      % generate approproiate lat/long arrays
      RLAT = RLAT31(1:nt,:,:,:);
      RLON = RLON31(1:nt,:,:,:);
      RLATb = squeeze(RLAT31(1:nt,1,:,:));
      RLONb = squeeze(RLON31(1:nt,1,:,:));

      % inverse barometer correction
      LM = double(squeeze(slm(1,:,:)));                 % land mask
      SM = double(LM*0);
      SM(find(LM == 0)) = 1;                    % sea mask
      dxyp = gridcellarea(nlat,nlon);
      dxyp2 = dxyp*ones(1,nlon);
      sea_area = sum(sum(dxyp2.*SM));
      %ps_ib = ps(1:nt,:,:);
      ps_ib = ps;
      if ib
        for ii=1:nt
          pstemp = double(squeeze(ps(ii,:,:)));
          ps_sea_ave = sum(sum(pstemp.*SM.*dxyp2))/sea_area;
          ps_ave = pstemp.*LM+ps_sea_ave*SM;
          ps_ib(ii,:,:) = ps_ave;
        end
      end

      % compute the contributions for each gridbox to the integrand.

      w_field = zeros(3,nt,nlev,nlat,nlon);    
      m_field = zeros(3,nt,nlat,nlon);    

      m_field(1,:,:,:) = ps_ib.*sin(RLATb).*cos(RLATb).*cos(RLATb).*cos(RLONb);
      m_field(2,:,:,:) = ps_ib.*sin(RLATb).*cos(RLATb).*cos(RLATb).*sin(RLONb);
      m_field(3,:,:,:) = ps_ib.*cos(RLATb).*cos(RLATb).*cos(RLATb);

      w_field(1,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*cos(RLON)-v.*cos(RLAT).*sin(RLON);
      w_field(2,:,:,:,:) = u.*sin(RLAT).*cos(RLAT).*sin(RLON)+v.*cos(RLAT).*cos(RLON);
      w_field(3,:,:,:,:) = u.*cos(RLAT).*cos(RLAT);


      % integrate both the wind and mass terms horizontally.
      int_m_dlat = -trapz(rlat,m_field,3);
      int_m_dlat_dlon = trapz(rlon,int_m_dlat,4);
      int_w_dlat = -trapz(rlat,w_field,4);
      int_w_dlat_dlon = trapz(rlon,int_w_dlat,5);
      int_w_dlat_dlon_dp = -trapz(lev,int_w_dlat_dlon,3);

      ww(:,k1:k2,:) = squeeze(int_w_dlat_dlon);

      % multiply by prefactors to compute EAFs in rad/s.
      Xw(:,k1:k2) = pw*ones(1,nt).*int_w_dlat_dlon_dp;
      Xm(:,k1:k2) = pm*ones(1,nt).*int_m_dlat_dlon;
      kold = k2;

      disp(['Finished   ',f])
    end % variables-not-screwed-up loop
  end % file exist if loop

end	% file loop


% this is the total AM 
X = Xw+Xm;

% if option is selected, save the 2-d contributing fields
% (i.e. with vertical integration in wind term)


% save output: 
% this makes a mat file with a generic name, which should be written over with the appropriate runid and decade.
savename = ['aam_vol_vert_slice',runid,'_',decade];
if ib == 0, savename = [savename,'_noib']; end
if savecomp, savename = [savename,'_long']; end
filename = [savename,'.mat'];

save(filename, 'Xw', 'Xm', 'MJD', 'ww', 'lev')
