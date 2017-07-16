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
%   22 May 2011: add option to compute correlation to observations
%   31 Aug 2011: change the filltering of the annual cycle to a more
%   generous range (230-450 d), to match the Nastula et al., 2009 paper.
%
%  INPUTS
%    runid: name of the dataset
%    tscale: which timescale
%       tscale = 1: 30-90 day (subseasonal) variations.
%       tscale = 2: 350-370d (annual cycle)
%       tscale = 3: 2-7 year variation.
%       tscale = 4: 7-20 year variations.
%    comp: which vector component of AAM (1,2,3) - set to zero to do a variable.
%    var_name: u, v, ps, or 'EF' to just do the local-to-global corr for that EF component..
%    total_name: name of the global variable to compute correlation for.
%    Choose 'EF' for excitation function, or 'ERP' for corresponding obs.
%    parameter.
%    fil_order : order of butterworth filter.  so far have been using 2.
%
%  OUTPUTS
%    RHO: field of correlations
%    SIG: field of covariances
%    LAG: field of lags where the above maximize.
function correlation_map( runid, tscale, comp, var_name, total_name, fil_order   )
%------------------------------------------------------------------

%---temp inputs-----------

%clear all;
%runid = 'CCMval';
%tscale = 1;
%comp = 3;
%var_name = 'EF';
%total_name = 'ERP';
%fil_order = 1;

%---special run settings
filtype = 1;

if strcmp(runid,'CCMVal') | strcmp(runid,'CCMval'); HRES = 'T42'; end

% set flag for simple variable correlation
if strcmp(var_name,'EF'), simple_var = 0; else simple_var = 1; end

% set flag for correlation with aam (rather than observations)
if strcmp(total_name,'ERP')
    corr_with_aam = 0;
    corr_with_ERP = 1;
else
    corr_with_aam = 1;
    corr_with_ERP = 0;
end

%---filtering settings

switch tscale350-
    case 0
        T = [0,0];      % (no filtering)
        date_start = [1960,1,1];
        date_stop = [1999,12,1];
        maxlag = 2*30;        % 2 MONTHS;
    case 1
        T = [30,90]; Ty = (1/365)*T;
        date_start = [1989,1,1];
        date_stop = [1999,12,1];
        maxlag = 2*30;        % 2 MONTHS;
        tscalename = 'Subseasonal'; 
    case 2
      T = [230-450];
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
     maxlag = round(180);	% ~6 months
     tscalename = 'Annual';
    case 3
      Ty = [2,7]; T = Ty*365;
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
      maxlag = round(1.5*365);	% 1.5 years
       tscalename = 'Interannual'
    case 4
      Ty = [7,20]; T = Ty*365;
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
      maxlag = 2*365;       % 2 years
      'Long-Term';
end


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

% initialize correlation / covariance matrix.
% for correlations between local EFs and a global number, RHO has 3 levels:
% mass, wind, and total.  For correlations with simple variables or
% observations, RHO is is simply latxlon.  Same for SIG and LAG.
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

if simple_var
    X = zeros(nd,nlat,nlon)+NaN;
    Xint = zeros(1,nd)+NaN;
else
  X = zeros(3,nd,nlat,nlon)+NaN;
  Xint = zeros(3,nd)+NaN;
end

for idec = 1:length(dec)
  % load either an AAM field or variable fields
  if simple_var
    dum = [sdir,'*uvpfields*.mat'];  
    ff = dir(dum);
    fname = [sdir,ff.name];
  else
    suffix = 's_long.mat';
    dum = [sdir,'*vol*',num2str(dec(idec)),suffix];
    disp(dum)
    ff = dir(dum);
    fname = [sdir,ff.name];
  end
  disp(['Loading file  ',fname])
  if exist(fname) == 2
    load(fname,'lat','lon','ww','mm','MJD')
    if simple_var == 0
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
          
          % if computing correlation with total, integrate here, and also
          % convert to equivalent LOD and PM.
          if corr_with_aam
            int_dlat = -trapz(rlat,xt,3);
            int_dlat_dlon = trapz(rlon,int_dlat,4);
            if comp < 3
                Xint(:,iday) = rad2mas*int_dlat_dlon;
            else
                Xint(:,iday) = LOD0_ms*int_dlat_dlon;
            end
          end
        end
      end
    else
        
      % If computing variable correlations, simply sort u, v, and ps into x.
      [c,targ,targ_big] = intersect(round(MJD),MJD_big);
      if strcmp(var_name,'u'), X(targ_big,:,:) = permute(U300(:,:,targ),[3,2,1]); end
      if strcmp(var_name,'v'), X(targ_big,:,:) = permute(V300(:,:,targ),[3,2,1]); end
      if strcmp(var_name,'ps'), X(targ_big,:,:) = permute(PS(:,:,targ),[3,2,1]); end


      % in this case we also need to load the total AAM timeseries, if
      % that's what we are computing the correlation with.
      % note that these are already nondimensionalized, so just multiply
      % by conversion to either equivalent PM or LOD.
      if corr_with_ERP
        dum = [sdir,'*vol*',num2str(dec(idec)),'s.mat'];
        ff = dir(dum);
        fname = [sdir,ff.name];
        load(fname); 
        [c,targ,targ_big] = intersect(round(MJD),MJD_big);
        if comp < 3
         Xint(targ_big) = rad2mas*(Xm(comp,targ)+Xw(comp,targ));
        else
         Xint(targ_big) = LOD_0*(Xm(comp,targ)+Xw(comp,targ));
        end  
      end
    end
  else
    disp(['Cant locate the necessary AAM input file  ' fname])
  end
end



% if computing correation or covariance with the observations, read these
% in here, then select only the overlap with the time axis loaded above.
% also note that the aam integrand components are not in equivalent PM or
% LOD units yet, but the conversion is performed below.
if corr_with_ERP
    Xint = zeros(1,length(MJD_big));     % empty array for obs
    [X1_obs,X2_obs,X3_obs,mjd_obs,eX1_obs,eX2_obs,eX3_obs] = read_eops;
    [c,targ,targ_big] = intersect(round(mjd_obs),MJD_big);
    switch comp
        case 1
            Xint(targ_big) = X1_obs(targ);
        case 2
            Xint(targ_big) = X2_obs(targ);
        case 3 
            Xint(targ_big) = X3_obs(targ);
    end
end
      if strcmp(var_name,'u'), X(targ_big,:,:) = permute(U300(:,:,targ),[3,2,1]); end


%--- detrend and filter the global timeseries (Xint) to the timescale at hand

F = Xint*0+NaN;
Tbot = T(1);
Ttop = T(2);
filter_interval = 1;

for ii = 1:size(Xint,1)
  g = find(isfinite(Xint(ii,:)));
  if comp < 3
    Xint_DT = Xint(ii,:) - nanmean(Xint(ii,:));
  else
    Xint_DT = Xint(ii,:) - nanmean(Xint(ii,:));
  end
  if isempty(g) == 0
    xfil = cfilter(Xint_DT(g),filter_interval,Tbot,Ttop,fil_order,'days');
    F(ii,g) = xfil(:,filtype);
  end
end

% (note that now the integral is already in mas or ms)


%--- loop through lats and lons and compute covariance with global value 
% for each gridbox.
if simple_var, nlevels =1; else nlevels = 3; end
for ii = 1:nlevels
  for ilat = 1:nlat
    for ilon = 1:nlon
      if simple_var
          g = find(isfinite(X(:,1,1))==1); 
          x1 = squeeze(X(g,ilat,ilon))';
      else
         g = find(isfinite(X(3,:,1,1))==1);  
         x1 = squeeze(X(ii,g,ilat,ilon))';
      end
      % Q: does it make a diff if we detrend the individual tseries
      x2 = detrend(x1,'constant');
      % if the local variable we're corellating is aam, scale to correct units.
      if simple_var == 0
        if comp < 3
          x = rad2mas*x2; 
        else
          x = LOD0_ms*x2; 
        end
      else
        x = x1;
      end 
      
      if corr_with_ERP
        [c,lag] = xcorr(x,F(g),maxlag,'coeff');    % correlation
        [c2,lag2] = xcov(x,F(g),maxlag);           % covariance
      else
        [c,lag] = xcorr(x,F(ii,g),maxlag,'coeff');    % correlation
        [c2,lag2] = xcov(x,F(ii,g),maxlag);           % covariance
      end
      
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


if corr_with_ERP
    savename1 = ['corr_',var_name,'_obsX',num2str(comp),'_',runid,'_',tscalename,'_filorder',num2str(fil_order)];
else
   savename1 = ['corr_',var_name,'_X',num2str(comp),'_',runid,'_',tscalename,'_filorder',num2str(fil_order)]; 
end
savename = [savename1,num2str(date_start(1)),'-',num2str(date_stop(1)),'.mat'];
save(savename, 'RHO', 'LAG', 'SIG', 'LAG2','lat', 'lon')

