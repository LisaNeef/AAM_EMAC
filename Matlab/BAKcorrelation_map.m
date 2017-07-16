function correlation_map( runid, tscale, comp, var_name, fil_order   )
%
% Load pre-computed AAM spatial data, filter to a desired timescale,
% and compute the covariance and correlation with the total.
%  9 Dec 1010
%  MODS:
%   26 Jan 2011: separate into mass and wind components
%    3 Feb 2011: added filter order as an input.
%    7 Feb 2011: instead of max corr/cov, look for max abs value
%	(so that we also get the negative correlations)
%   20 May 2011: option to compute correlation with respect to variables.
%
%  INPUTS
%    runid: name of the dataset
%    tscale: which timescale
%       tscale = 1: 30-90 day (subseasonal) variations.
%       tscale = 2: 350-370d (annual cycle)
%       tscale = 3: 2-7 year variation.
%       tscale = 4: 7-20 year variations.
%    comp: which vector component of AAM (1,2,3) - set to zero to do a variable.
%    var_name: u, v, ps, or 'aam' to just do the local-to-global corr for that aam component..
%    fil_order : order of butterworth filter.  so far have been using 2.
%
%  OUTPUTS
%    RHO: field of correlations
%    SIG: field of covariances
%    LAG: field of lags where the above maximize.
%------------------------------------------------------------------

%---temp inputs-----------

%clear all;
%runid = 'CCMval';
%tscale = 2;
%comp = 3;
%var_name = 'u';
%fil_order = 1;

%---special run settings
filtype = 1;

if strcmp(runid,'CCMVal') | strcmp(runid,'CCMval'); HRES = 'T42'; end

% set flag for simple variable correlation
if strcmp(var_name,'aam') | strcmp(var_name,'AAM'); simple_var = 0; else simple_var = 1; end

%---filtering settings

if tscale == 1
  T = [30,90]; Ty = (1/365)*T;
  date_start = [1989,1,1];
  date_stop = [1999,12,1];
  maxlag = 2*30;        % 2 MONTHS;
end
if tscale == 2
  T = [350,370];
  date_start = [1960,1,1];
  date_stop = [1999,12,1];
  maxlag = round(180);	% ~6 months
end
if tscale == 3
  Ty = [2,7]; T = Ty*365;
  date_start = [1960,1,1];
  date_stop = [1999,12,1];
  maxlag = round(1.5*365);	% 1.5 years
end
if tscale == 4
  Ty = [7,20]; T = Ty*365;
  date_start = [1960,1,1];
  date_stop = [1999,12,1];
  maxlag = 2*365;       % 2 years
end

if tscale == 1, tscalename = 'Subseasonal'; end
if tscale == 2, tscalename = 'Annual'; end
if tscale == 3, tscalename = 'Interannual'; end
if tscale == 4, tscalename = 'Long-Term'; end

%--- load mat files containing integrand components

% make a time array for the dates considered.
mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
MJD_big = mjd0:1:mjdf;


% prepare output arrays 

if HRES(1:3) == 'T42'
  nlat = 64;
  nlon = 128;
end

if simple_var
  RHO = zeros(nlat,nlon)+NaN;
  SIG = zeros(nlat,nlon)+NaN;
  LAG = zeros(nlat,nlon)+NaN;
  LAG2 = zeros(nlat,nlon)+NaN;
else
  RHO = zeros(3,nlat,nlon)+NaN;
  SIG = zeros(3,nlat,nlon)+NaN;
  LAG = zeros(3,nlat,nlon)+NaN;
  LAG2 = zeros(3,nlat,nlon)+NaN;
end
% load conversions
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)

% prefactors (Gross 09)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));
pm3  = 0.997*0.750*Re^4/(g*Cm);
pw3  = 0.997*Re^3/(Q*g*Cm);
pm = [pm12; pm12; pm3];
pw = [pw12; pw12; pw3];


%--- cycle through giant array and collect regional contributions
% here also nondimensionalize using prefactors (notes vol. 4, p. 52c)

dec0 =  10*floor(date_start(1)/10);
decf =  10*floor(date_stop(1)/10);
dec = dec0:10:decf;

sdir = ['/dsk/nathan/lisa/EMAC_ex/',char(runid),'/mat/'];
d2r = pi/180.;                   %  degrees to radian

if strcmp(var_name,'aam') | strcmp(var_name,'AAM')
  X = zeros(3,nd,nlat,nlon)+NaN;
  Xint = zeros(3,nd)+NaN;
else
  X = zeros(nd,nlat,nlon)+NaN;
  Xint = zeros(1,nd)+NaN;
end

for idec = 1:length(dec)
  % load either an AAM field or variable fields
  if strcmp(var_name,'aam') | strcmp(var_name,'AAM')
    suffix = 's_long.mat';
    dum = [sdir,'*vol*',num2str(dec(idec)),suffix];
    disp(dum)
    ff = dir(dum);
    fname = [sdir,ff.name];
  else
    dum = [sdir,'*uvpfields*.mat'];  
    ff = dir(dum);
    fname = [sdir,ff.name];
  end
  disp(['Loading file  ',fname])
  if exist(fname) == 2
    load(fname,'lat','lon','ww','mm','MJD')
    if strcmp(var_name,'aam') | strcmp(var_name,'AAM')
      xx = zeros(3,length(MJD),nlat,nlon)+NaN;
      % compute wind, mass, and total AAM contributions per grid cell in rad (or rad/s)
      xx(1,:,:,:) = squeeze(pw(comp)*ww(comp,:,:,:));
      xx(2,:,:,:) = squeeze(pm(comp)*mm(comp,:,:,:));
      xx(3,:,:,:) = squeeze(pm(comp)*mm(comp,:,:,:)+pw(comp)*ww(comp,:,:,:));
      rlat = lat*d2r;
      rlon = lon*d2r;
      for iday = 1:nd
        targ = find(round(MJD) == MJD_big(iday));
        if isempty(targ) == 0
          xt = xx(:,targ(1),:,:);
          X(:,iday,:,:) = xt;
          int_dlat = -trapz(rlat,xt,3);
          int_dlat_dlon = trapz(rlon,int_dlat,4);
          Xint(:,iday) = int_dlat_dlon;
        end
      end
    else
      % If computing variable correlations, simply sort u, v, and ps into x.
      [c,targ,targ_big] = intersect(round(MJD),MJD_big);
      if strcmp(var_name,'u'), X(targ_big,:,:) = permute(U300(:,:,targ),[3,2,1]); end
      if strcmp(var_name,'v'), X(targ_big,:,:) = permute(V300(:,:,targ),[3,2,1]); end
      if strcmp(var_name,'ps'), X(targ_big,:,:) = permute(PS(:,:,targ),[3,2,1]); end


      % in this case we also need to load the total AAM timeseries.
      % note that these are already nondimensionalized.
      dum = [sdir,'*vol*',num2str(dec(idec)),'s.mat'];
      ff = dir(dum);
      fname = [sdir,ff.name];
      load(fname); 
      [c,targ,targ_big] = intersect(round(MJD),MJD_big);
      if comp < 3
        Xint(targ_big) = Xm(comp,targ)+Xw(comp,targ);
      else
        Xint(targ_big) = Xm(comp,targ)+Xw(comp,targ);
      end  
    end
  else
    disp(['Cant locate the necessary AAM input file  ' fname])
  end
end


%--- scale, detrend, and filter the integrated timeseries to the timescale at hand

F = Xint*0+NaN;
Tbot = T(1);
Ttop = T(2);
filter_interval = 1;

for ii = 1:size(Xint,1)
  g = find(isfinite(Xint(ii,:)));
  if comp < 3
    Xint_scaled = rad2mas*(Xint(ii,:) - nanmean(Xint(ii,:)));
  else
    Xint_scaled = LOD0_ms*(Xint(ii,:) - nanmean(Xint(ii,:)));
  end
  if isempty(g) == 0
    xfil = cfilter(Xint_scaled(g),filter_interval,Tbot,Ttop,fil_order,'days');
    F(ii,g) = xfil(:,filtype);
  end
end

% (note that now the integral is already in mas or ms)


%--- loop through lats and lons and compute covariance with total in each gridbox
for ii = 1:size(Xint,1)
  for ilat = 1:nlat
    for ilon = 1:nlon
      if simple_var, g = find(isfinite(X(:,1,1))); else g = find(isfinite(Xint(ii,:))); end
      if simple_var, x1 = squeeze(X(g,ilat,ilon))'; else x1 = squeeze(X(ii,g,ilat,ilon))'; end
      % Q: does it make a diff if we detrend the individual tseries
      x2 = detrend(x1,'constant');
      % if the variable we're corellating is aam, scale to correct units.
      if simple_var == 0
        if comp < 3
          x = rad2mas*x2; 
        else
          x = LOD0_ms*x2; 
        end
      else
        x = x1;
      end 
      [c,lag] = xcorr(x,F(ii,g),maxlag,'coeff');    % correlation
      [c2,lag2] = xcov(x,F(ii,g),maxlag);           % covariance
      % how to compute the correlation and covariance -- zero lag, or at
      % maximum?
      % these are the maximal values
      %[a,b] = max(abs(c));
      %[a2,b2] = max(abs(c2));
      % these are the central values.
      b = find(lag == 0);
      b2 = find(lag2 == 0);
      
      if simple_var
        RHO(ilat,ilon) = c(b);
        SIG(ilat,ilon) = c2(b2);
        LAG(ilat,ilon) = lag(b);
        LAG2(ilat,ilon) = lag2(b2);
      else
        RHO(ii,ilat,ilon) = c(b);
        SIG(ii,ilat,ilon) = c2(b2);
        LAG(ii,ilat,ilon) = lag(b);
        LAG2(ii,ilat,ilon) = lag2(b2);
      end
    end
  end
end

%--- shift the lons over to get it right for the map 
% (lons should go -180 to 180)
if min(lon) <= 0
  dum = find(lon > 180);
  lon2=lon;
  lon2(dum) = -(360 - lon(dum));
  [a,b] = sort(lon2);
  lon = a;
end

% don't forget to also shift the data!
if simple_var
  RHO = RHO(:,b);
  SIG = SIG(:,b);
  LAG = LAG(:,b);
  LAG2 = LAG2(:,b);
else
  RHO = RHO(:,:,b);
  SIG = SIG(:,:,b);
  LAG = LAG(:,:,b);
  LAG2 = LAG2(:,:,b);
end

%--- save the computed dataset

% save output:
% this makes a mat file with a generic name, which should be written over with the appropriate runid and decade.


savename1 = ['corr_',var_name,'_X',num2str(comp),'_',runid,'_',tscalename,'_filorder',num2str(fil_order)];
savename = [savename1,num2str(date_start(1)),'-',num2str(date_stop(1)),'.mat'];
save(savename, 'RHO', 'LAG', 'SIG', 'LAG2','lat', 'lon')

