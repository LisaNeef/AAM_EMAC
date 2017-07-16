function [corr_out,corr_ci,corr_std,Nsamples, tau,runs_out,C,L,R2] = correlations_aam_filtered_bootstrap(comp,tscale,alpha)
% Compute the correlations between different AAM datasets 
% for a given component.
% 
% v3: instead of computing a correlation significance limit, use a
% bootstrap approach to compute the confidence interval of each
% correlation.
% 1 July 2011
%
% MODS:
%   1 Sep 2011: add option of annual variations, now that we've channged
%   the filter for annual variations, so that it's consistent with the window used in
%   Nastula et al., 2009
%   1 Sep 2011: do a block-bootstrap based on decorrelation time for
%   each timescale (See notes for Aug 31st).
%
% OUTPUTS:
%   corr_out: max correlation found within the prescribed lags.
%   corr_ci: correlation confidence intervals for each run
%   corr_std: coorelation standard deviation from bootstrap.
%   tau: esimated decorrelation time for each run - note that this is
%   actually neglected in the bootstrap computation.  
%   runs_out: names of the runs for which correlations were computed.
%   C: cross-correlations for each run as a function of lags
%   L: vector of lags.
%   R2: pointers to the runs for which correlations were computed.
%
%
% INPUTS:
%   comp: 1, 2, or 3
%   tscale: which timescale to look at.
%	tscale = 1: 30-90 day (subseasonal) variations.
%	tscale = 2: 24-30 month variations (QBO).
%	tscale = 3: 2-7 year variation.
%	tscale = 4: 7-20 year variations.

%--------------------------------------------------------------------------

%---temp
clear all;
comp = 3;
tscale = 1;
alpha = .05;

%--how many bootstrap samples?
nboot = 1000;


%--------------------

% select the timescales to examine here,
% plus the runs we want to focus on for each timescale.
if tscale == 1	% subseasonal
  T = [30,90]; 
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
if tscale == 5      % annual
    T = [230,450];
    maxlag = 365;
end

date_start = [1960,1,1];
date_stop = [2008,12,31];


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




%---load EMAC runs to compare to GEO
% and select only the ones that are relevant.
[runs,names,ef,no_ib] = aam_paper_runs;
nf = length(runs);

    R = [1:8,10];	% runs to retrieve
    R2 = 1:8;		% runs to compute correlation for.

% set pointers to where OAM and HAM and CAM are.
% for consistency, use ERA-40 here.
%  oam = 7;
%  ham = 8;
%  cam = 9;
%  net = 11;
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
  [X,XF,~] = retrieve_AAM_filtered(runs(irun),T,filtype,fil_order,date_start,date_stop, ef(irun),no_ib(irun),'p');
  XX(irun,:) = squeeze(X(3,comp,:));
  XXF(irun,:) = squeeze(XF(3,comp,:));
end

%----single out the ERP variations ("GEO")

XFG = squeeze(XXF(geo,:));
%XG = squeeze(XX(geo,:));

%---compute correlations with the geodetic residual here---------------

% initialize output arrays
% these only hold statistics for the runs that aren't GEO
nR2 = length(R2);
corr_out = zeros(nR2,1);
corr_ci = zeros(nR2,2);
corr_std = zeros(nR2,1);
Nsamples = zeros(nR2,1);

nc = maxlag*2+1; 

C = zeros(nR2,nc);
L = zeros(nR2,nc);	4;

% clear some shit out so that matlab doesn't have a cow.
clear XX
clear XG
clear X
clear XF
clear MJD

for ii = 1:nR2
  irun = R2(ii);
  
  % find the points where each run overlaps with the observations and
  % neither one has a bad value.
  finite_test = isfinite(XXF(irun,:)) + isfinite(XFG);
  x = squeeze(XXF(irun,finite_test == 2));
  y = XFG(finite_test == 2);
  %[c,lag] = xcorr(x,y,maxlag,'coeff');
  %C(ii,:) = c;
  %L(ii,:) = lag; 
    
  % block-bootstrap confidence interval:
  %--old---[B,bootsam] = bootstrp(nboot,@corr,x',y');
  
  % Option 1: simple subsampling bootstrap.
  %[B] = bstrap(nboot,4*tau/length(x),'corr',[x',y']); 
  %B3 = squeeze(B(1,2,:));
  %Bsort = sort(B3);  %Rank the coefficients from low to high
  
  % Option 2: stationary bootstrap (blocks of varying sizes)
  [bsdata,~] = stationary_bootstrap((1:length(x))',nboot,4*tau);
  B = zeros(nboot,1);
  for iboot = 1:nboot
      B(iboot) = corr(x(bsdata(:,iboot))',y(bsdata(:,iboot))');
  end
  Bsort = sort(B);
  
  corr_ci(ii,:) = [Bsort(round((alpha/2)*nboot)); Bsort(round((1-alpha/2)*nboot))] ; %Report confidence intevals
  corr_out(ii) = mean(Bsort);
  corr_std(ii) = std(Bsort);
  
  % plot the histogram of estimated bootstrap correlations to see how good
  % the sample is. 
  %figure(2+ii)
  %histfit(B3,nboot/20);
  %title(['TimeScale ',num2str(tscale),'  Run ',names(R2(ii))])
  
  %  compute effective nr of DOF and estimate min corr. for significance
  %[c_auto,lag_auto] = xcorr(x,x,maxlag,'coeff');
  %[~,b2] = min(abs(c_auto));
  %tau(ii) = abs(lag_auto(b2));


end


runs_out = names(R2);



