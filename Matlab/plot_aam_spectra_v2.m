%function plot_aam_spectra_v2(comp,terms,plot_type,BW,plot_letter,annotate)
%
% Compare the frequency spectra of AAM in EMAC runs, reanalyses, and observations.
% Feb 2011
%
% Version 2:
% Instead of a choice of fft types, just use Christof's simple
% implementation.
% 
%
% INPUTS:
%   comp: 1, 2, or 3
%   terms: which terms in the AAM interval: wind (1), mass (2), total (3)
%   plot_type: set of runs to plot for specific purpose.
%	1: compare EMAC runs to ERA-40: do resolution / time-dependent BCs make a difference?
%	4: compare ERA-40/Interim to the raw observed signals, plus OAM and HAM.
%   BW: set to 1 for black and white plots.
%   xlabel_on: set to 1 to label the xaxis.
%   plot_letter: the label for the plot being made (for the paper)
%   annotate: add boxes indicating the important timescales.

%---------------------------------------------------------------

clear all;

% choose component and which paper plot type we want.
comp = 3;
BW = 1;

%pp = 'X1t';
%pp = 'X1m';
%pp = 'X1w';

%pp = 'X2t';
%pp = 'X2m';
%pp = 'X2w';

%pp = 'X3t';
%pp = 'X3m';
pp = 'X3w';

switch pp
    case 'X1t'
        plot_type = 4;
        plot_letter = '(a) \chi_1 Observation v. Excitation';
        terms = 3;
        xlabel_on = 0;
        comp = 1;
        annotate = 1;
    case 'X1m'
        plot_type = 1;
        plot_letter = '(c) \chi_1 AAM Mass Term';
        terms = 2;
        xlabel_on = 0;
        comp = 1;
        annotate = 0;
    case 'X1w'
        plot_type = 1;
        plot_letter = '(e) \chi_1 AAM Wind Term';
        terms = 1;
        xlabel_on = 1;
        comp = 1;
        annotate = 0;
    case 'X2t'
        plot_type = 4;
        plot_letter = '(b) \chi_2 Observation v. Excitation';
        terms = 3;
        xlabel_on = 0;
        comp = 2;
        annotate = 1;
    case 'X2m'
        plot_type = 1;
        plot_letter = '(d) \chi_2 AAM Mass Term';
        terms = 2;
        xlabel_on = 0;
        comp = 2;
        annotate = 0;
    case 'X2w'
        plot_type = 1;
        plot_letter = '(f) \chi_2 AAM Wind Term';
        terms = 1;
        xlabel_on = 1;
        comp = 2;
        annotate = 0;
    case 'X3t'
        plot_type = 4;
        plot_letter = '(a) \chi_3 Observation v. Excitation';
        terms = 3;
        xlabel_on = 0;
        comp = 3;
        annotate = 1;
    case 'X3m'
        plot_type = 1;
        plot_letter = '(b) \chi_3 AAM Mass Term';
        terms = 2;
        xlabel_on = 0;
        comp = 3;
        annotate = 0;
    case 'X3w'
        plot_type = 1;
        plot_letter = '(c) \chi_3 AAM Wind Term';
        terms = 1;
        xlabel_on = 1;
        comp = 3;
        annotate = 0;
end
   


BW = 1;
addpath('/home/ig/neef/MFiles/utilities/');

%---define the runs to be plotted here
switch plot_type
    case 1
        R = [2,4,5];
    case 2
        R = 4:5;
    case 3
        R = [5:8,10,11,19];
    case 4 
        R = [10,5,7,8];
    case 5
        R = [4,5,10];
    case 6
        R = [4,5,7,8,10,11];
end


[runs,names_total,run_flag,no_ib] = aam_paper_runs;
names = names_total(R);
nf = length(runs);
nR = length(R);

names = names_total(R);
nf = length(runs);

% default settings for which terms to plot.
if plot_type == 3|plot_type == 4|plot_type == 5| plot_type == 6, terms = 3; end

date_start = [1960,1,1];
date_stop = [1999,12,1];

% define filter settings
  fo = 2;               % greater than 2 seems to not work for ERA data
  ft = 1;                       % butterworth filter

%---initialize empty array

mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
XX = zeros(length(terms),nf,nd)+NaN;

%---cycle through datasets and extract timeseries for the total (wind plus pressure) of each component

% set some dummy filter options.
T = [1 20]*365;
fil_order = 1;	filtype = 1;
net1 = NaN;	net2 = NaN;
oam1 = NaN;	oam2 = NaN;
ham = NaN;
geo = NaN;
for irun = R
  [X,XF,MJD] =  retrieve_AAM_filtered(runs(irun),T,ft,fo,date_start,date_stop, run_flag(irun),no_ib(irun),'p');
  for iterm = terms
    XX(iterm,irun,:) = squeeze(X(iterm,comp,:));
  end

  % also note done where OAM, HAM, and CAM live so we can compute the residual AEF
  nn = char(runs(irun));
  if length(nn) == 5
    if nn(1:5) == 'OAM40', oam1 = irun; end
    if nn(1:5) == 'OAMec', oam2 = irun; end
    if nn(1:5) == 'HAM40', ham = irun; end
    if nn(1:5) == 'CAMs4', cam = irun; end
  end
  if length(nn) == 4
    if nn(1:4) == 'NET1', net1 = irun; end
    if nn(1:4) == 'NET2', net2 = irun; end
  end
  if nn(1:3) == 'GEO', geo = irun; end

end

%--compute AEF residual from GEO (axial terms only).

if isfinite(net1(1)) > 0
  XX(3,net1,:,:) = XX(3,geo,:,:)-XX(3,oam1,:,:)-XX(3,ham,:,:);
end
if isfinite(net2(1)) > 0
  XX(3,net2,:,:) = XX(3,geo,:,:)-XX(3,oam2,:,:)-XX(3,ham,:,:);
end

%---cycle through timeseries and compute spectra

%----new---------

FS = 365;                        % sample frequency in years^(-1)
L = size(XX,3);                  % Length of timeseries
t = (0:L-1)/FS;                  % Period vector in years.

NFFT = 2^nextpow2(L);               % Next power of 2 from length of y
F = FS/2*linspace(0,1,NFFT/2+1);    % frequency vector.
T = (F.^-1)*365.5;                          % period vector -- in days

P = zeros(3,nf,length(F))+NaN;	% empty array to hold spectra

for irun = R
  for iterm = terms
    Xrun = squeeze(XX(iterm,irun,:));
    gg = find(isfinite(Xrun)==1);
    Ydum = fft(Xrun(gg),NFFT)/L;
    Y = 2*abs(Ydum(1:NFFT/2+1,:));
    P(iterm,irun,:) = Y;
  end
end

%---plot settings
% for comparison to geo AEF residual, re-select the runs we show
switch plot_type 
  case 3
    R2 = [5,6];         % Plots ERAs only - add residual with shading later
  case 6
    R2 = [4,5,11];         % Plot only CCMVal, ERA-40, and NET1
  otherwise,
    nR2 = nR;
    R2 = R;
end
nR2 = length(R2);


col = aam_paper_colors;
transparency = 0.3;
LS = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};

%---black and white settings
if BW
  col(R,:) = zeros(length(R),3);
  gray = 0.5*ones(1,3);
  era = 5; oam = 7; ham = 8;
  if plot_type == 4
     LS(ham) = cellstr('-.'); 
     LS(oam) = cellstr('--'); 
     col(era,:) = gray;
     col(ham,:) = gray;
     col(oam,:) = gray;
  end	
  if plot_type == 1
    LS(era) = cellstr('--'); 
    dum = bone(6);%	(create a black & white colormap)
    col(1:4,:) = dum(1:4,:);
  end
else
  %shade_color = [0.8117 0.9019 0.1254];
  if isfinite(geo),  shade_color = col(geo,:); end
  % change the color of the CCMVal run to distinguish on the printer from GEO
  if plot_type == 5, col(4,:) = [0.2785    0.5469    0.9575]; end
  if plot_type == 6, col(4,:) = [0.2785    0.5469    0.9575]; end
  if plot_type == 2, col(4,:) = [0.2785    0.5469    0.9575]; end
end

lh = zeros(length(R2),1);

rend = 'opengl';
LW = 5;

% set axis limits
  xbot =  25/365;
  xtop =  10;

% specify the log limits of the y-axis
switch comp
  case 1           % X1, X2
    ybot = -3;
    ytop = 2.5;
  case 2           % X1, X2
    ybot = -3;
    ytop = 2.5;
  case 3        % X3
    ybot = -6;
    ytop = 1;
end 

axx = [xbot xtop 10^ybot 10^ytop];
y_log_steps = ybot:2:ytop;
yticks = 10.^y_log_steps;

TTw = ['\chi_' num2str(comp) ' Wind Term'];
TTm = ['\chi_' num2str(comp) ' Mass Term'];
TTt = ['\chi_' num2str(comp)];
TT = {TTw;TTm;TTt};

% annotation box settings
subs = find(T > 30 & T<90);
xsubs  = [30 90]/365;	
wsubs = xsubs(2)-xsubs(1);
dum = squeeze(P(terms,R,subs));	
ysubs = [0.9*min(min(dum)) 5*max(max(dum))];

inta = find(T > 2*365 & T < 7*365);
xinta  = [2 7];
winta = xinta(2)-xinta(1);
dum = squeeze(P(terms,R,inta));	
yinta = [0.8*min(min(dum)) 4*max(max(dum))];

ann = find(T> 360 & T < 370);
xann = [1 1];
if comp == 3
    yann = [500*axx(3) 0.0001*axx(4)];
end
if comp == 1|2
    yann = [50*axx(3) 0.001*axx(4)];
end

%---make plots
for iterm = terms
  figure(iterm),clf,hold off
      for irun = 1:nR2
      to_plot = squeeze(P(iterm,R2(irun),:)) ;
      lh(irun) = loglog(T/365.5,to_plot,'Linewidth',LW,'Color',col(R2(irun),:),'LineStyle',char(LS(R2(irun)))); 
      hold on
    end 
    % for plot_type 3, add shadded residual estimates.
    if plot_type == 3
      net = [net1,net2];
      p_net = squeeze(P(iterm,net,:));
      pnet_bot = squeeze(min(p_net));
      pnet_top = squeeze(max(p_net));
      %[h,msg] = jbfill(T/365.5,pnet_bot,pnet_top,1.4*col(net1,:),0.7*col(net1,:),1,transparency);
      [h,msg] = jbfill(T/365.5,pnet_bot,pnet_top,0.7*ones(1,3),0.7*col(net1,:),1,transparency);
      lh(irun+1) = h;
    end

    ylabel('PSD') 
    if xlabel_on, xlabel('Period (Years)'); end
    names_here = names_total(R2);
    if plot_type == 3, names_here = [names_total(R2);'Obs. AEF Residual']; end
    if plot_type == 5, names_here = {'CCMVal','ERA-40','Geodetic'}; end
    if plot_type == 6, names_here = {'Unconstrained Model AAM','ERA-40 AAM','Inferred from Obs.'}; end
    gl = find(isfinite(lh));
    legend(lh(gl),names_here(gl),'Orientation','Horizontal','Location','SouthEast')
    legend('boxoff')
    axis(axx);
    set(gca,'YTick',yticks)
    % add plot labels:
    if BW
      text(1.1*axx(1),0.3*axx(4),plot_letter,'FontSize',30)
    end
    % add a title to the top plot in the column
    if plot_type == 4 & BW == 1
     % title(['\chi_',num2str(comp)])
    end

    % option to annotate plot with boxes showing timescales of interest.
    if annotate
      % convert box limits to figure units
      [xsubsf,ysubsf] = ds2nfu(xsubs,ysubs);
      [xintaf,yintaf] = ds2nfu(xinta,yinta);
      [xannf,yannf] = ds2nfu(xann,yann);
      [xsannf,ysannf] = ds2nfu(0.5*xann,yann);
      wsubsf = xsubsf(2) - xsubsf(1);
      hsubsf = ysubsf(2) - ysubsf(1);
      wintaf = xintaf(2) - xintaf(1);
      hintaf = yintaf(2) - yintaf(1);
      subseasonal = [xsubsf(1) ysubsf(1) wsubsf hsubsf];
      interannual= [xintaf(1) yintaf(1) wintaf hintaf];
      % choose color, linestyles, etc.
      annotate_color = 0.7*ones(1,3);
      annotation('rectangle',subseasonal,'Color',annotate_color,'LineWidth',4)
      annotation('rectangle',interannual,'Color',annotate_color,'LineWidth',4)
      text(xsubs(1)+0.4*wsubs,3*ysubs(2),'Subseasonal','Color',annotate_color,'FontSize',30,'HorizontalAlignment','center')
      text(xinta(1)+0.4*winta,0.1*yinta(1),'Interannual','Color',annotate_color,'FontSize',30,'HorizontalAlignment','center')
      annotation('textarrow',xannf,yannf,'Color',zeros(1,3),'String','Annual')
      annotation('textarrow',xsannf,ysannf,'Color',zeros(1,3),'String','Semiannual')
    end
end

%---export plots

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
switch plot_type
    case 1
        suff1 = ['_emacruns']; 
    case 2
        suff1 = ['_CCMVal_vs_ERA']; 
    case 3
        suff1 = ['_ERA_vs_NET']; 
    case 4
        suff1 = ['_ERA_vs_GEO']; 
    case 5
        suff1 = ['_CCMVal_ERA_GEO']; 
    case 6
        suff1 = ['_CCMVal_ERA_NET']; 
end

fig_name_w = [plot_dir,'spectr_X',num2str(comp),'w',suff1,'_fft.png']; 
fig_name_m = [plot_dir,'spectr_X',num2str(comp),'m',suff1,'_fft.png']; 
fig_name_t = [plot_dir,'spectr_X',num2str(comp),'t',suff1,'_fft.png']; 
FN = {fig_name_w;fig_name_m;fig_name_t};

rend = 'opengl';
LW = 2;
ph = 6;        % paper height
pw = 20;        % paper width
fs = 25;        % fontsize

for iterm = terms
  exportfig(iterm,char(FN(iterm)),'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');
end






