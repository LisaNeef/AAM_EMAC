function [max_corr,max_corr_lag,Nsamples, tau,runs_out,C,L,R2, corr_sig] = correlations_aam_filtered_v2(comp,tscale,exp_type)
% Compute the correlations between different AAM datasets 
% for a given component.
% v2: retrieve the timescale-filtered runs with an external function.
% 6 Dec 2010
% Mods:
%  8 Dec 2010: apply Lillietest to find out whether samples are Gaussian-distributed (they generally aren't)
% 14 Dec 2010: added exp_type 2
% 16 Dec 2010: eradicated some problems with NaNs.
% 20 Dec 2010: estimate lowest significant correlation for each case (corr_sig) and add to output.  	
% 10 Mar 2011: compute ALL correlations across the max period possible (up to 2008 for ERAinterim)
% 11 Mar 2011: instead of outputting significance limit, outout decorrelation time per run, and # of samples
%		per run.
%		Also, simplify by removing the central correlation as an output.
% 28 Jun 2011: change Nsamples to be the *effective* number of independent
% samples after accounting for lag and autocorrelation.  modify the
% significance limit to reflect both autocorrelation AND lag.
%
% TO DO:
%  - check how decorrelation times compare with lags in each case -- how does Neff vary from case to case?
%	(see notes vol. 4, p. 64c)
%
% OUTPUTS:
%   max_corr: max correlation found within the prescribed lags.
%   max_corr_lag: lag corresponding to the max correlation
%   Nsamples: effective number of samples after adjusting for
%   autocorrelation and lag.
%   tau: decorrelation time computed for each dataset & filtering
%   timescale.
%   runs_out: names of the runs for which correlations were computed.
%   C: correlations for each run as a function of lags
%   L: vector of lags.
%   R2: pointers to the runs for which correlations were computed.
%   corr_sig: the minimum correlation required in order to be statistically
%   significant. Has a -1 if no number is possible because the optimal lag
%   was greater than the autocorrelation.
%
% INPUTS:
%   comp: 1, 2, or 3
%   tscale: which timescale to look at.
%	tscale = 1: 30-90 day (subseasonal) variations.
%	tscale = 2: 24-30 month variations (QBO).
%	tscale = 3: 2-7 year variation.
%	tscale = 4: 7-20 year variations.
%  exp_type: which correlation do we want to compute?
%	1: correlations between EMAC runs and ERA datasets and GEO.
%       2: correlations between CAM datasets and GEO (notes vol. 4, p. 50)
%	3: correlations between all datasets and GEO minus CAM (for long-term)
%--------------------------------------------------------------------------



% select the timescales to examine here,
% plus the runs we want to focus on for each timescale.
if tscale == 1	% subseasonal
  T = [30,90]; Ty = (1/365)*T;
  maxlag = 2*30;	% 2 MONTHS;
end
if tscale == 2   % QB
  Tm = [23,34];
  Ty = Tm*(1/12); T = Ty*365;
  maxlag = 2*365;   % 2 years.

end
if tscale == 3 % interannual
  Ty = [2,7]; T = Ty*365;
  maxlag = 2*365;   % 2 years.
end
if tscale == 4
  Ty = [7,20]; T = Ty*365;
  maxlag = 3*365;	% 3 years
end

date_start = [1960,1,1];
date_stop = [2008,12,31];

%---load EMAC runs to compare to GEO
% and select only the ones that are relevant.
[runs,names,ef,no_ib] = aam_paper_runs;
nf = length(runs);
no_filter = 0;
switch exp_type
  case 1
    R = [1:8,10];	% runs to retrieve
    R2 = [1:8];		% runs to compute correlation for.
  case 2
    R2 = [9,15:17];
    R = [R2, 10];
    no_filter = 1;
  case 3 
    R = [1:11];
    R2 = [1:9];
end


% set pointers to where OAM and HAM and CAM are.
% for consistency, use ERA-40 here.
  oam = 7;
  ham = 8;
  cam = 9;
  net = 11;
  geo = 10;

% cycle through runs, load data, and filter out interannual
fil_order = 2;		% greater than 2 seems to not work for ERA data

% choose type of filter to plot
filtype = 1;		% butterworth filter
%filtype = 2;		% chebyeshev-1 filter
%filtype = 3;		% chebyeshev-2 filter


% initialize empty array

mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
XX = zeros(nf,nd);
XXF = zeros(nf,nd);

% retrieve the filtered runs

for irun = R
  [X,XF,MJD] = retrieve_AAM_filtered(runs(irun),T,filtype,fil_order,date_start,date_stop, ef(irun),no_ib(irun),'p');
  XX(irun,:) = squeeze(X(3,comp,:));
  XXF(irun,:) = squeeze(XF(3,comp,:));
end

%----single out the ERP variations ("GEO")

XFG = squeeze(XXF(geo,:));
XG = squeeze(XX(geo,:));

% for plot type 3, remove CAM
if exp_type == 3
    XFG =squeeze(XXF(geo,:))  - XXF(cam,:); 
    XG = squeeze(XX(geo,:)) - XX(cam,:); 
end

%--decide whether to use filtered or unfiltered timeseries
if no_filter
  XFG = XG;
  XXF = XX;
end

%---compute correlations with the geodetic residual here---------------

% initialize output arrays
% these only hold statistics for the runs that aren't GEO
nR2 = length(R2);
max_corr = zeros(nR2,1);
max_corr_lag = zeros(nR2,1);
tau = zeros(nR2,1);
Nsamples = zeros(nR2,1);

nc = maxlag*2+1; 

C = zeros(nR2,nc);
L = zeros(nR2,nc);	


for ii = 1:nR2
  irun = R2(ii);
  Xrun = squeeze(XXF(irun,:));		% isolate the run at hand
  % ..but for CAM runs, don't use the filtered series.
  if ((R2(ii) == 9) | (R2(ii) == 15) | (R2(ii) == 16) | (R2(ii) == 17)), 
    Xrun = squeeze(XX(irun,:));
  end
  finite_test = isfinite(Xrun) + isfinite(XFG);
  g = find(finite_test == 2);
  x = Xrun(g);
  y = XFG(g);
  [c,lag] = xcorr(x,y,maxlag,'coeff');
  C(ii,:) = c;
  L(ii,:) = lag;  
  
  % find local maxima in the correlation.
  [pks,locs] = findpeaks(abs(c));
  [pks2,locs2] = sort(pks);
  
  % choose the largest max, unless it is at the boundary (=maxlag), in
  % which case take the next largest one.
  %npeaks = length(locs2);
  %pmax = max(pks2);
  %lag_at_max = lag(locs(max(locs2)));
  %[pmax,locmax] = max(pks);
  %if pmax == maxlag
  %    pmax = pks(npeaks-1);
  %    lag_at_max = lag(locs2(npeaks-1));
  %end
  %max_corr(ii) = pmax;
  %max_corr_lag(ii) = lag_at_max;
 
  % alternative: choose the local max with the smallest lag.  that will be
  % more statistically significant, though the correlation may be small.
  [smallest_lag,dum] = min(abs(lag(locs)));
  rho_smallest_lag = pks(dum);
  max_corr(ii) = rho_smallest_lag;
  max_corr_lag(ii) = lag(locs(dum));

  %  compute effective nr of DOF and estimate min corr. for significance
  [c_auto,lag_auto] = xcorr(x,x,maxlag,'coeff');
  [a2,b2] = min(abs(c_auto));
  decorrelation_time = abs(lag_auto(b2));
  tau(ii) = decorrelation_time;
  
  % note that the decorrelation timesteps are months for tscale = 4, not days
  if tscale < 4, nd = decorrelation_time; else nd = decorrelation_time/30; end
  N = length(g);
  disp(names(R(ii)));
  disp(['decorrelation time estimate ' num2str(nd)])
  disp(['sample size: '  num2str(N)])
  %Neff = N/abs(nd) - 2*max_corr_lag(ii);
  
  % effective number of independent samples: remove the lags and divide out
  % the decorrelation time of that dataset.
  Neff = (N - 2*max_corr_lag(ii))/abs(nd);
  Nsamples(ii) = Neff;
  disp(['effective sample size: '  num2str(Neff)])
  corr_sig(ii) = erfinv(0.98)*sqrt(2/Neff);
  
  % if the lag that corresponds to the max correlation exceeds the
  % autocorrelation timescale, we don't have a significant correlation.  In
  % that case, set corr_sig to -1, to indicate no significant correlation
  % possible.
  if Neff < 0, corr_sig(ii) = -1; end
end

% for CAM series, the lag is in months (because data are monthly), so 
% convert back to days
days2months = 30;
cam_runs =  find ((R2 == 9) | (R2 == 15) | (R2 == 16) | (R2 == 17));
L(cam_runs,:) = days2months*L(cam_runs,:);
max_corr_lag(cam_runs) = days2months*max_corr_lag(cam_runs);


runs_out = names(R2);



