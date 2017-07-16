%
% Load pre-computed AAM spatial data, filter to a desired timescale,
% and compute the covariance and correlation with the total.
% In this version, use a bootstrap method to compute the correlations,
% since they aren't really Gaussian.
% 6 July 2011
%
%  MODS
%   31 Aug 2011: Change the Annual Cycle filter to the more generous
%   230-450 days, to make it consistent with Nastula et al., 2009.
%   31 Aug 2011: also, do a block-bootstrap based on decorrelation time.
%   See notes in revision journal for this day.
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
%    nboot: number of bootstrap samples.
%
%  OUTPUTS - saved to a file.
%    RHO: field of correlations
%    R_LO: lower bound of the correlation 98% confidence interval
%    R_HI: upper bound of the correlation 98% confidence interval
function correlation_map_bootstrap( runid, tscale, comp, var_name, total_name, fil_order, nboot, alpha   )
%------------------------------------------------------------------

%---temp inputs-----------

clear all;
runid = 'CCMval';
tscale = 1;
comp = 3;
var_name = 'EF';
total_name = 'EF';
fil_order = 1;
nboot = 500;
alpha = 0.05;



%---special run settings
filtype = 1;

if strcmp(runid,'CCMVal') || strcmp(runid,'CCMval'); HRES = 'T42'; end

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

switch tscale
    case 0
        T = [0,0];      % (no filtering)
        date_start = [1960,1,1];
        date_stop = [1999,12,1];
    case 1
        T = [30,90]; 
        date_start = [1989,1,1];
        date_stop = [1999,12,1];
        tscalename = 'Subseasonal'; 
    case 2
      T = [230,450];
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
     tscalename = 'Annual';
    case 3
      Ty = [2,7]; 
      T = Ty*365;
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
       tscalename = 'Interannual';
    case 4
      Ty = [7,20]; 
      T = Ty*365;
      date_start = [1960,1,1];
      date_stop = [1999,12,1];
      'Long-Term';
end

% choose assumed decorrelation time (in days) based on component and timescale
% note that when we apply this below, the data HAVE to be daily.
switch tscale
    case 0
        tau_eq = 1;
        tau_ax = 1;
    case 1
        tau_eq = 11;
        tau_ax = 13;
    case 2
        tau_eq = 80;
        tau_ax = 90;
    case 3
        tau_eq = 370;
        tau_ax = 160;
end
if comp == 1||comp==2, tau = tau_eq; end
if comp == 3, tau = tau_ax; end

% make a time array for the dates considered.
mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
MJD_big = mjd0:1:mjdf;


% prepare output arrays 

if strcmp(HRES(1:3),'T42')
  nlat = 64;
  nlon = 128;
end

% initialize correlation / covariance matrix.
% for correlations between local EFs and a global number, RHO has 3 levels:
% mass, wind, and total.  For correlations with simple variables or
% observations, RHO is is simply latxlon.  Same for the others.
if simple_var
  RHO = zeros(nlat,nlon)+NaN;
  R_LO = zeros(nlat,nlon)+NaN;
  R_HI = zeros(nlat,nlon)+NaN;
else
  RHO = zeros(3,nlat,nlon)+NaN;
  R_LO = zeros(3,nlat,nlon)+NaN;
  R_HI = zeros(3,nlat,nlon)+NaN;
end

% load conversions
load aam_constants_gross
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
  if exist(fname,'file') == 2
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
      [~,targ,targ_big] = intersect(round(MJD),MJD_big);
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
        [~,targ,targ_big] = intersect(round(MJD),MJD_big);
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
    [X1_obs,X2_obs,X3_obs,mjd_obs] = read_eops;
    [~,targ,targ_big] = intersect(round(mjd_obs),MJD_big);
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
    disp('Computing correlations for component')
    disp(ii)
  for ilat = 1:nlat
      disp(ilat)
    for ilon = 1:nlon
      if simple_var
          g = find(isfinite(X(:,1,1))==1); 
          x1 = squeeze(X(g,ilat,ilon))';
      else
         g = find(isfinite(X(3,:,1,1))==1);  
         x1 = squeeze(X(ii,g,ilat,ilon))';
      end
      % Q: does it make a diff if we detrend the individual tseries
      %x2 = detrend(x1,'constant');
      x1_fil = cfilter(x1',filter_interval,Tbot,Ttop,fil_order,'days');
      x2 = x1_fil(:,filtype);
      
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
          % correlation with observations
          y = F(g)';
      else
          % correlation with global AEF - respective wind or mass term
          % only.
          y = F(ii,g)';
      end
      
      
      
      % bootstrap confidence interval:
      % [B] = bootstrp(nboot,@corr,x,y);
      % -----Option 1: simple subsampling bootstrap.
      %[B] = bstrap(nboot,4*tau/length(x),'corr',[x,y]); 
      %B3 = squeeze(B(1,2,:));
      %Bsort = sort(B3);  %Rank the coefficients from low to high
      
      % -----Option 2: stationary bootstrap (blocks of varying sizes)
      [bsdata,~] = stationary_bootstrap((1:length(x))',nboot,4*tau);
      B = zeros(nboot,1);
      for iboot = 1:nboot
          B(iboot) = corr(x(bsdata(:,iboot)),y(bsdata(:,iboot)));
      end
      Bsort = sort(B);
      
      
      corr_ci = [Bsort(round((alpha/2)*nboot),:); Bsort(round((1-alpha/2)*nboot),:)] ; %Report confidence intevals
      corr_out = mean(Bsort);
      %corr_std = std(B(:,2));
            
      if simple_var
        RHO(ilat,ilon) = corr_out;
        R_LO(ilat,ilon) = corr_ci(1);
        R_HI(ilat,ilon) = corr_ci(2);
      else
        RHO(ii,ilat,ilon) = corr_out;
        R_LO(ii,ilat,ilon) = corr_ci(1);
        R_HI(ii,ilat,ilon) = corr_ci(2);        
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
  R_LO = R_LO(:,b);
  R_HI = R_HI(:,b);
else
  RHO = RHO(:,:,b);
  R_LO = R_LO(:,:,b);
  R_HI = R_HI(:,:,b);
end

%--- save the computed dataset

% save output:
% this makes a mat file with a generic name, which should be written over with the appropriate runid and decade.


if corr_with_ERP
    savename1 = ['corr_',var_name,'_obsX',num2str(comp),'_',runid,'_',tscalename,'_filorder',num2str(fil_order)];
else
   savename1 = ['corr_',var_name,'_X',num2str(comp),'_',runid,'_',tscalename,'_filorder',num2str(fil_order)]; 
end
savename = [savename1,num2str(date_start(1)),'-',num2str(date_stop(1)),'_bootstrap',num2str(nboot),'.mat'];
save(savename, 'RHO', 'R_LO', 'R_HI', 'lat', 'lon','nboot','alpha')

