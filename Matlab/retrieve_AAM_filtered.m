function [X,XF,MJDout] = retrieve_AAM_filtered(runid,T,filtype,fil_order,date_start,date_stop, emac_flag,no_IB,ERP_var_type)
%
% select an excitation functiond data set over a given timespan,
% filtered over desired period.
% started 1 Dec 2010
% MODS:
%  6 Dec 2010: add no_IB flag (so far only for emac datasets)
%  1 Feb 2011: add details for OAM from ECCO
% 14 Mar 2011: modify paths to EMAC run data.
% 23 Mar 2011: fixed a small bug in loading the ECCO OAM dataset
% 31 Mar 2011: add the option of lo-pass filtering only
% 13 May 2011: add the option of outputting uncertainties (for the obs sets)
%
% INPUTS:
%  runid : choose EMAC run name, CCMVal, ERAin, ERA40, OAMin, OAM40, HAMin, HAM40, 
%            CAMt1, CAMt4, CAMs1, CAMs4, GEO
%  T: array holding bottom and top cutoff periods.
%  filtype: 1 for Butterworth, 2 for Chebyeshev-1
%  fil_order: order of the filter.
%  date_start = [y,m,d] where we begin
%  date_stop = [y,m,d] where we cut off
%  emac_flag: set to 1 if the dataset is an EMAC run (incl. ccmval)
%  no_IB: set to 1 to load the dataset were we have no IB approximation
%  ERP_var_type: set to 'p' to read the ERPs themselves, or 'u' to read in their uncertaintes.
%	(for other datasets, the input for this variable doesn't matter.)
%
% OUTPUTS:
%  X: 3x3xndays array.  
%	col 1 is the wind term, mass term, and total
%	col 2 is for chi1, chi2 (both in mas) and dlod (in ms)
%  MJDout: time array out in modified julian day.
%
% TO-DO:
%    - add no_IB flag to non-EMAC datasets


%-------------------------------------------------

% make a time array for the dates considered.
mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
MJDout = mjd0:1:mjdf;

% make output arrays of filtered and unfiltered AAM
X = zeros(3,3,nd)+NaN;
XF = zeros(3,3,nd)+NaN;

% load the data
run_name = char(runid);
data_type = 'no data type found';
if emac_flag, data_type = 'EMAC'; end
if strcmp(run_name(1:3), 'ERA'), data_type = 'RA'; end
if strcmp(run_name(1:3), 'OAM'), data_type = 'RA'; end
if strcmp(run_name(1:3), 'HAM'), data_type = 'RA'; end
if strcmp(run_name, 'GEO'), data_type = 'GEO'; end
if strcmp(run_name, 'NET'), data_type = 'GEO'; end
if strcmp(run_name(1:3), 'CAM'), data_type = 'CAM'; end
if strcmp(run_name,'UNC'), data_type = 'GEO'; end
disp(data_type)

% load conversions
aam_constants_gross

switch data_type

  % (1) -----------------EMAC runs-----------------------------
  case 'EMAC'

    % load all the relevant decades and slap together into a timeseries
    dec0 =  10*floor(date_start(1)/10);
    decf =  10*floor(date_stop(1)/10);
    dec = dec0:10:decf;

    xm = zeros(3,1)+NaN;
    xw = zeros(3,1)+NaN;
    mjd = NaN;
    sdir = ['/dsk/nathan/lisa/EMAC_ex/',char(runid),'/mat/'];

    for idec = 1:length(dec)
      if no_IB, suffix = 's_noIB.mat'; else suffix = 's.mat'; end
      ff = dir([sdir,'*vol*',num2str(dec(idec)),suffix]);
      fname = [sdir,ff.name];
      if exist(fname) == 2
        load(fname)
        xm = [xm,Xm];
        xw = [xw,Xw];
        mjd = [mjd,MJD];
      end
    end

    % error message if no data is there.
    if size(xm,2) == 1, disp(['Unable to find any files for run ' runid]), return, end
 
    % take out NaNs and detrend
    gw = find(isfinite(xw(1,:)) == 1);
    gm = find(isfinite(xm(1,:)) == 1);
    xw_DT = xw*0; 	xm_DT = xm*0;
    for ii = 1:3
      xw_DT(ii,gw) = detrend(xw(ii,gw),'constant');
      xm_DT(ii,gm) = detrend(xm(ii,gm),'constant');
    end

    % fill in arrays for the dates specified
    for iday = 1:nd
      target = find(floor(mjd) == MJDout(iday));
      if length(target) > 0
        tt = target(1);
        X(1,1:2,iday) = rad2mas*xw_DT(1:2,tt);
        X(1,3,iday) = LOD0_ms*xw_DT(3,tt);
        X(2,1:2,iday) = rad2mas*xm_DT(1:2,tt);
        X(2,3,iday) = LOD0_ms*xm_DT(3,tt);
      end
    end
    X(3,:,:) = X(1,:,:)+X(2,:,:);

  % (2) -----------------Reanalysis and other "outside" datasets-----------------------------
  case 'RA'
    
    % choose ERAinterim or ERA40.
    if run_name(4:5) == 'in'; RA = 'ERAinterim'; end
    if run_name(4:5) == '40'; RA = 'ERA40'; end
    if run_name(4:5) == 'ec'; RA = 'ECCO'; end

    if run_name(1:3) == 'ERA', AM = 'aam'; end
    if run_name(1:3) == 'OAM', AM = 'oam'; end
    if run_name(1:3) == 'HAM', AM = 'ham'; end

    dailyvalues = 1;
    % special case: to read in ECCO for long timescales, use the 10-day average data.
    % kind of a messy trick, this dataset is called by setting dailyvalues to 0.
    % more messiness: when we do the entire spectrum (T(2) = 20*365) then go back to daily vals.
    if RA(1:4) == 'ECCO' & T(1) > 700 & T(2) < 20*300, dailyvalues = 0; else, dailyvalues = 1; end
    [xw,xm,mjd_RA]= read_EFs(AM,RA,dailyvalues);

    % detrend!
    gw = find(isfinite(xw(1,:)) == 1);
    gm = find(isfinite(xm(1,:)) == 1);
    xw_DT = xw*0; 	xm_DT = xm*0;
    for ii = 1:3
      xw_DT(ii,gw) = detrend(xw(ii,gw),'constant');
      xm_DT(ii,gm) = detrend(xm(ii,gm),'constant');
    end

    % fill in arrays for the dates specified
    for iday = 1:nd
      target = find(floor(mjd_RA) == MJDout(iday));
      if length(target) > 0
        tt = target(1);
        X(1,1:2,iday) = rad2mas*xw_DT(1:2,tt);
        X(1,3,iday) = LOD0_ms*xw_DT(3,tt);
        X(2,1:2,iday) = rad2mas*xm_DT(1:2,tt);
        X(2,3,iday) = LOD0_ms*xm_DT(3,tt);
      end
    end
    X(3,:,:) = X(1,:,:)+X(2,:,:);

  % (3) -----------------ERP observations-----------------------------
  case 'GEO'
    
    % read in the data  (already has the right units)
    [x1,x2,dlod,mjd_geo,ex1,ex2,edlod] = read_eops;

    % if the real ERPs, take out the long-term mean value first.
    if strcmp(ERP_var_type,'p')
      x1_DT = detrend(x1,'constant');
      x2_DT = detrend(x2,'constant');
      dlod_DT = detrend(dlod,'constant');
    end

    % fill in arrays for the dates specified
    % (since we only have the 'total' for GEO, leave rows 1 and 2 as NaNs.)
    for iday = 1:nd
      target = find(floor(mjd_geo) == MJDout(iday));
      if length(target) > 0
        tt = target(1);
        % either read in the detrended parameters....
        if strcmp(ERP_var_type,'p')
          X(3,1,iday) = x1_DT(tt);
          X(3,2,iday) = x2_DT(tt);
          X(3,3,iday) = dlod_DT(tt);
        end
        % ... or read in their uncertainties instead.
        if strcmp(ERP_var_type,'u')
          X(3,1,iday) = ex1(tt);
          X(3,2,iday) = ex2(tt);
          X(3,3,iday) = edlod(tt);
        end
      end
    end

  % (4) -----------------CAM from Jan Hagedroon-----------------------------
  case 'CAM'
    
    % choose which CAM timeseries
    topog = run_name(4); 
    cond = run_name(5);
    
    % read in.
    [x1,x2,dlod,mjd_cam,date_cam] = read_CAM_JH(topog,cond);

    % detrend
    x1_DT = detrend(x1,'constant');
    x2_DT = detrend(x2,'constant');
    dlod_DT = detrend(dlod,'constant');

    % fill in arrays for the dates specified
    % (since we only have the 'total' for CAM, leave rows 1 and 2 as NaNs.)
    for iday = 1:nd
      target = find(floor(mjd_cam) == MJDout(iday));
      if length(target) > 0
        tt = target(1);
        X(3,1,iday) = x1_DT(tt);
        X(3,2,iday) = x2_DT(tt);
        X(3,3,iday) = dlod_DT(tt);
      end
    end

end


%--now filter the timeseries to the desired timescale.

XF = X*0+NaN;
Tbot = T(1);
if T(1) == 0, Tbot = []; end
Ttop = T(2);
if T(2) == 0, Ttop = []; end
filter_interval = 1/(MJDout(2)-MJDout(1));

for iterm = 1:3
  for ivec = 1:3
    x = squeeze(X(iterm,ivec,:));
    g = find(isfinite(x));
    if length(g) > 0
      xf = cfilter(x(g),filter_interval,Tbot,Ttop,fil_order,'days');
      XF(iterm,ivec,g) = xf(:,filtype);
    end
  end
end
