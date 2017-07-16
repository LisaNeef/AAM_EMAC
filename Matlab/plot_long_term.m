function [h1,h2] = plot_long_term(comp, apply_lag,export)
%
% Make a plot showing the full evolution of ERPs, compared to the CAM data 
% supplied by Jan Hagedoorn.
% Then make another plot comparing the residual ("GEO-CAM") to atmospherie data
% from ERA-40, ERA-interim, and the CCMVal emac run.
%
% 13 March 2011
%
% MODS:
%  31 Mar 2011:  compute the residual for the obs with a lo-pass filter
%	also subtract net with the optimal lag.
%   4 Jun 2011: add the option not to apply a lag, based on discussion with
%   Jan H.
%
% INPUTS:
%   comp: vector component to plot.
%   apply_lag: set to 1 to apply a lag to the CAM timeseries when we compute the residual.
%------------------------------------------------------------


% temp inputs
%comp = 3;
%apply_lag = 0; 

%---define the runs that we want to look at here, and set up empty arrays

R = [4,5,6,9,10,11];
[runs,names_total,emac_flag,no_ib] = aam_paper_runs;
nR = length(R);

date_start = [1962,1,1];
date_stop = [2000,12,31];
mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;

XX = zeros(nR,nd);
XXF = zeros(nR,nd);

% filter settings - here do a low pass filter
fo = 2;           
ft = 1;                     
T = [5,0]*365; 


%---retrieve the (unfiltered) datasets

for irun = R
  % retrieve filtered timeseries
  [X,XF,MJD] =  retrieve_AAM_filtered(runs(irun),T,ft,fo,date_start,date_stop, emac_flag(irun),no_ib(irun),'p');
  XX(irun,:) = squeeze(X(3,comp,:));	% (total only, forget about wind and mass terms here)
  XXF(irun,:) = squeeze(XF(3,comp,:));	% (total only, forget about wind and mass terms here)
end

%---lag the cam timeseries by the optimal lag.
cam = 9;
net = 11;
geo = 10;

dum = squeeze(XX(cam,:));
g = find(isfinite(dum));
xcam = squeeze(XX(cam,g));
xgeo = squeeze(XX(geo,g));
[c,lag] = xcorr(xcam,xgeo,'coeff');
[maxcorr,b]=max(c);

xcam2 = xcam*0+NaN;
if apply_lag
    xcam_shift = circshift(xcam',-lag(b))';
    xcam2(1:length(xcam)-lag(b)-1) = xcam_shift(1:length(xcam_shift)-lag(b)-1);
else
    xcam2 = xcam;
end





%---compute the residual of GEO and CAMs4:


% compute the residual for the filtered geo signal.
XX(net,g) = squeeze(XX(geo,g)) - xcam2;
XXF(net,g) = squeeze(XXF(geo,g)) - xcam2;

% make the time axis
t = MJD*0+NaN;
[y, m, d] = mjd2date(MJD);
for ii=1:length(t)
  if isfinite(MJD(ii)) ==1, t(ii)=datenum([y(ii) m(ii) d(ii)]); end
end


%---plot settings

R1 = [10,9];	% pointer to GEO and CAMs4
R2 = [11,5,4];  % pointer to NET, ERA, CCMVal
nR1 = length(R1);
nR2 = length(R2);

switch comp
  case 1
    ax1 = 200*[-1,1];
    ax2 = 200*[-1,1];
  case 2
    ax1 = 200*[-1,1];
    ax2 = 200*[-1,1];
  case 3
    ax1 = 3.5*[-1,1];
    ax2 = [-0.4,1];
end


%---compute some correlations and print that shit out!

% correlation between CEF and DLOD
if apply_lag
    disp(['Correlation between CAM and DLOD:  ',num2str(maxcorr)])
    disp(['...with a lag of :  ',num2str(lag(b)/12), 'years'])
else
    disp(['Correlation between CAM and DLOD:  ',num2str(c(lag==0))])
    disp('...with no lag.')
end

% correlation between CEF and ERA40, CCMVal
era = 5;
ccmval = 4;
xera = squeeze(XXF(era,g));
xccmval = squeeze(XXF(ccmval,g));
xnet = squeeze(XXF(net,g));
[c_era,lag_era] = xcorr(xera,xnet,'coeff'); 
g2 = find(isfinite(xccmval));
[c_ccmval,lag_ccmval] = xcorr(xccmval(g2),xnet(g2),'coeff'); 

[a_era,b_era] = max(c_era);
[a_ccmval,b_ccmval] = max(c_ccmval);


if apply_lag
    disp(['Correlation between CAM and ERA-40:  ',num2str(a_era)])
    disp(['...with a lag of :  ',num2str(lag_era(b_era)/12), 'years'])
    disp(['Correlation between CAM and CCMVal:  ',num2str(a_ccmval)])
    disp(['...with a lag of :  ',num2str(lag_ccmval(b_ccmval)/12), 'years'])    
else
    disp(['Correlation between CAM and ERA-40:  ',num2str(c_era(lag_era==0))])
    disp(['Correlation between CAM and CCMVal:  ',num2str(c_ccmval(lag_ccmval==0))])
    disp('...with no lag.')
end

% define colors to be used here

col_obs = 0*ones(1,3);   
col_cam = 0.7*ones(1,3);
col_era = 0*ones(1,3);
col_emac = 0*ones(1,3);

col1 = [col_obs;col_cam];
col2 = [col_obs;col_era;col_emac];

% define linestyles to be used here.

LS_obs = '-';
%LS_cam = '-';
LS_era = '--';
LS_emac = '-.';

LS2 = {LS_obs;LS_era;LS_emac};

LW = 2;

% also define the names to be plotted here.

names1 = {'OBS';'CAM'};
names2 = {'OBS-CAM';'ERA-40';'CCMVal'};

%---plot 1: comparison of the timeseries: observed and CAM

figure(1),clf

h1 = subplot(2,1,1);
lh = zeros(nR1,1);
for irun = 1:nR1
  x = squeeze(XX(R1(irun),:));
  g3 = find(isfinite(x));
  tg = t(g3);
  lh(irun) = plot(tg,x(g3),'Color',col1(irun,:),'LineWidth',LW);
  hold on
end
legend(lh,names1,'Orientation','Horizontal')
  top = ax1(2);
  bot = ax1(1);
  dax = (top-bot)/10;
  set(gca,'YTick',bot:dax:top)
  legend('boxoff')
  xlabel('year')
  datetick('x','yyyy')
  ylabel('ms')
  axis([t(1),max(t), ax1(1), ax1(2)]);


%---plot 2: comparison of the residual to the AEFs


h2 = subplot(2,1,2);
for irun = 1:nR2
  x = squeeze(XXF(R2(irun),:));
  tg = t(g);
  lh(irun) = plot(tg,x(g),'Color',col2(irun,:),'LineStyle',char(LS2(irun)),'LineWidth',LW);
  hold on
  %names_here = {'CCMVal','ERA-40','Residual'};
end
legend(lh,names2,'Orientation','Horizontal')
  top = ax1(2);
  bot = ax1(1);
  dax = (top-bot)/10;
  set(gca,'YTick',bot:dax:top)
  legend('boxoff')
  xlabel('year')
  ylabel('ms')
  datetick('x','yyyy')
  if apply_lag
      axis([t(1),max(t), ax2(1), ax2(2)]);
  else
      axis([t(1),max(t), 0.5*ax1(1), 0.5*ax1(2)]);
  end


%---export:
if export
%rend = 'opengl';
LW = 2;
ph = 9;        % paper height
pw = 12;        % paper width
fs = 16;        % fontsize

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
fig_name = [plot_dir,'long_term_X',num2str(comp),'.png'];

exportfig(1,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
end

