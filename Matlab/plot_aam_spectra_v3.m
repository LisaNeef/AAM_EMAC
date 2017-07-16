%function plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate,export,extra_labels,label_x_axis,label_y_axis)
%
% Compare the frequency spectra of AAM in EMAC runs, reanalyses, and observations.
% Feb 2011
%
% V3: instead of plotting power, plot power*frequency -- that highlights
% the main signals more.
%
% INPUTS:
%   comp: 1, 2, or 3
%   terms: which terms in the AAM interval: wind (1), mass (2), total (3)
%   plot_type: set of runs to plot for specific purpose.
%	1: compare EMAC runs to ERA-40: do resolution / time-dependent BCs make a difference?
%	2: compare CCMVal run to ERA-40/ERA-interim
%	3: compare ERA-40/Interim to 2 estimates of the geodetic AEF residual.
%	4: compare ERA-40/Interim to the raw observed signals, plus OAM and HAM.
%	5: compare ERA-40 only to the raw observed signals AND CCMVal
%	6: like 5, but instead of GEO show residual - only 1 because ECCO is to short to properly 
%          compute the spectrum
%   spect_type: which type of spectrum to plot:
%	1: pwelch
%	2: periodogram
% 	3: pmtm
%   BW: set to 1 for black and white plots.
%   plot_letter: the label for the plot being made (for the paper)
%   annotate: add boxes indicating the important timescales.
%   export: set to 1 to export a png version of this plot, 0 to not do
%   that.
%   extra_labels: set to 0 to remove some of the annotation that is only
%   needed in the first plot column.
%
%---------------------------------------------------------------

% temp:
clear all;
annotate = 0;
terms = 3;
plot_letter = ' ';
plot_type = 4;
spect_type = 1;
BW = 0;
extra_labels = 1;
comp = 3;
export = 0;
figure(1),clf
label_x_axis = 1;
label_y_axis = 1;

addpath('/home/ig/neef/MFiles/utilities/');

%---define the runs to be plotted here
switch plot_type
    case 1
        R = 1:5;
    case 2
        R = 4:5;
    case 3
        R = [5:8,10,11,19];
    case 4 
        R = [5,7,8,10];
    case 5
        R = [4,5,10];
    case 6
        R = [4,5,7,8,10,11];
end


[runs,names_total,emac_flag,no_ib] = aam_paper_runs;
nf = length(runs);

% default settings for which terms to plot.
if plot_type == 3||plot_type == 4||plot_type == 5|| plot_type == 6, terms = 3; end

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
net1 = NaN;	net2 = NaN;
oam1 = NaN;	oam2 = NaN;
ham = NaN;
geo = NaN;
for irun = R
  [X,~,~] =  retrieve_AAM_filtered(runs(irun),T,ft,fo,date_start,date_stop, emac_flag(irun),no_ib(irun),'p');
  for iterm = terms
    XX(iterm,irun,:) = squeeze(X(iterm,comp,:));
  end

  % also note done where OAM, HAM, and CAM live so we can compute the residual AEF
  nn = char(runs(irun));
    if strcmp(nn, 'OAM40'), oam1 = irun; end
    if strcmp(nn, 'OAMec'), oam2 = irun; end
    if strcmp(nn, 'HAM40'), ham = irun; end
    if strcmp(nn,'NET1'), net1 = irun; end
    if strcmp(nn,'NET2'), net2 = irun; end
    if strcmp(nn,'GEO'), geo = irun; end
    if strcmp(nn,'ERA40'),era = irun; end
end

%--compute AEF residual from GEO (axial terms only).

if isfinite(net1(1)) > 0
  XX(3,net1,:,:) = XX(3,geo,:,:)-XX(3,oam1,:,:)-XX(3,ham,:,:);
end
if isfinite(net2(1)) > 0
  XX(3,net2,:,:) = XX(3,geo,:,:)-XX(3,oam2,:,:)-XX(3,ham,:,:);
end

%---cycle through timeseries and compute spectra

Tmin = 5;	% minimum period in days
%Tmax = nd/3;	% max period in days.
Tmax = 15*365.5;
T = Tmin:1:Tmax;
F = 1./T;	% frequency range in cycles/d
Fs = 1;	% sampling frequency in 1/d



P = zeros(3,nf,length(F))+NaN;	% empty array to hold spectra

for irun = R
  for iterm = terms
    Xrun = squeeze(XX(iterm,irun,:));
    gg = find(isfinite(Xrun)==1);
    switch spect_type
      case 1 % pwelch with default settings, freq in cpd
        p = pwelch(Xrun(gg),[],[],F,Fs); 
      case 2 % periodogram
        p = periodogram(Xrun(gg),[],F,Fs);
      case 3 % multi-taper method
        p = pmtm(Xrun(gg),[],F,Fs);
    end
    P(iterm,irun,:) = p;
  end
end

%---plot settings
% for comparison to geo AEF residual, re-select the runs we show
switch plot_type 
  case 3,
    R2 = [5,6];         % Plots ERAs only - add residual with shading later
  case 6,
    R2 = [4,5,11];         % Plot only CCMVal, ERA-40, and NET1
  case 4,
    R2 = [5,7,8];
  otherwise,
    R2 = R;
end
nR2 = length(R2);


%---color & linstyle settings

LW = 2;
FS = 14;

col = aam_paper_colors;
transparency = 0.3;
LS = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
gray = 0.7*ones(1,3);

if BW
  col(R,:) = zeros(length(R),3);
  gray = 0.5*ones(1,3);
  era = 5; oam = 7; ham = 8;
  if plot_type == 4
     LS(era) = cellstr('--');
     LS(ham) = cellstr('-.'); 
     LS(oam) = cellstr('-'); 
     col(era,:) = zeros(1,3);
     col(ham,:) = gray;
     col(oam,:) = gray;
  end	
  if plot_type == 1
    LS(era) = cellstr('--'); 
    dum = bone(6);%	(create a black & white colormap)
    col(1:4,:) = dum(1:4,:);
  end
else
  if isfinite(geo),  shade_color = col(geo,:); end
  % change the color of the CCMVal run to distinguish on the printer from GEO
  if plot_type == 5, col(4,:) = [0.2785    0.5469    0.9575]; end
  if plot_type == 6, col(4,:) = [0.2785    0.5469    0.9575]; end
  if plot_type == 2, col(4,:) = [0.2785    0.5469    0.9575]; end
  if plot_type == 4
      col(R2,:) = winter(length(R2));
      col(era,:) = zeros(1,3);
  end
end

lh = zeros(length(R2),1);

% set axis limits
  %xtop =  log10(max(T)/365);
  xtop =  log10(max(T));
  xbot = -xtop;
  xbot = log10(min(T));

% specify the log limits of the y-axis
% axes for plots that compare different excitation sources

if plot_type == 4
  switch comp
    case {1,2}
      ybot = -2;
      ytop = 7;
    case 3
      ybot = -6;
      ytop = 1;  
  end
  y_steps = ybot:1:ytop;
end

if plot_type == 1
  switch comp
    case {1,2}
      ybot = 2;
      ytop = 7;
    case 3
      ybot = -3;
      ytop = 3.4;  
  end
  y_steps = ybot:1:ytop;
end

axx = [xbot xtop ybot ytop];

yticks = y_steps;
yticks(length(yticks)) = NaN;
if ~label_x_axis
  yticks(1) = NaN;
end

% annotation box settings
subs = T > 30 & T<90;
xsubs  = log10([30 90]);	
wsubs = xsubs(2)-xsubs(1);
dum = squeeze(P(terms,R,subs));	
ysubs = 1.1*log10([min(min(dum)) max(max(dum))]);

inta = T > 2*365 & T < 7*365;
xinta  = log10([2 7]*365);
winta = xinta(2)-xinta(1);
dum = squeeze(P(terms,R,inta));	
yinta = log10([0.9*min(min(dum)) 1.1*max(max(dum))]);

ann = find(T> 360 & T < 370);
xann = log10([1 1]*365);
dum = squeeze(P(terms,R,ann));
dum2 = min(min(dum));
if dum2 > 0
    yann = [log10(dum2)-5, log10(dum2)-2];
else
    yann = [.2 .8]*log10(min(min(dum)));
end

sann = round(ann/2);
xsann = log10([0.5,0.5]*365);
dum = squeeze(P(terms,R,sann));
if dum2 > 0
    ysann = [log10(dum2)-5, log10(dum2)-3];
else
    ysann = [.2 .8]*log10(min(min(dum)));
end


%---make plots
for iterm = terms
  if export, figure(iterm),clf; end
    axis(axx);
    yy = get(gca,'YLim');
    dylim = yy(2)-yy(1);
    %xlim = get(gca,'XLim');

    if annotate
        boundedline(xsubs,[yy(2),yy(2)],[dylim dylim],'cmap',0.7*flipud(bone)); hold on
        boundedline(xinta,[yy(2),yy(2)],[dylim dylim],'cmap',0.7*flipud(bone)); hold on    
    end
    
  if (plot_type == 4) 
    
    to_plot = squeeze(P(iterm,geo,:)).*F' ; 
    if BW
        lh(nR2+1) = boundedline(log10(T),log10(to_plot),[-ybot+log10(to_plot),zeros(length(F),1)],'alpha','cmap',bone); hold on
    else
        lh(nR2+1) = boundedline(log10(T),log10(to_plot),[-ybot+log10(to_plot),zeros(length(F),1)],'cmap',0.7*flipud(summer)); hold on
    end
  end

    for irun = 1:nR2
      to_plot = squeeze(P(iterm,R2(irun),:)).*F' ;
      lh(irun) = plot(log10(T),log10(to_plot),'Color',col(R2(irun),:),'LineStyle',char(LS(R2(irun))),'LineWidth',LW); 
      hold on
    end
    
    % reapply axis in case it didn't get applied above
    axis(axx)
    yy = get(gca,'YLim');
    dylim = yy(2)-yy(1);
    xlim = get(gca,'XLim');
    
    % for plot_type 3, add shadded residual estimates.
    if plot_type == 3
      net = [net1,net2];
      p_net = squeeze(P(iterm,net,:));
      pnet_bot = squeeze(min(p_net));
      pnet_top = squeeze(max(p_net));
      %[h,msg] = jbfill(T/365.5,pnet_bot,pnet_top,1.4*col(net1,:),0.7*col(net1,:),1,transparency);
      [h,~] = jbfill(log10(T),log10(pnet_bot),log10(pnet_top),shade_color,0.7*col(net1,:),1,transparency);
      lh(irun+1) = h;
    end
    
    
    % add vertical lines denoting interannual and
    % subseasonal bands
    if annotate
        plot(xann,yy,'--','Color',gray,'LineWidth',2)
        plot(xsann,yy,'--','Color',gray,'LineWidth',2)
    end

    if label_y_axis
       ylabel('Log PSD * frequency')
    else
       set(gca,'YTickLabel','') 
    end
    if label_x_axis
      xlabel('days')
      set(gca,'XTick',log10([10,50,100,300,500,1000,2500,5000]));
      set(gca,'XTickLabel',{'10','50','100','300','500','1000','2500','5000'})
    else
      set(gca,'XTickLabel','')  
    end
    gl = find(isfinite(lh));
    switch plot_type
        case 1
          names_here = names_total(R2);
          if extra_labels
             legend(lh(gl),names_here(gl),'Orientation','Vertical','Location','NorthEast')
          end
        case 3
           names_here = [names_total(R2);'Obs. AEF Residual'];
        case 4
            names_here = {'AAM','OAM','HAM','OBS'}; 
            legend(lh(gl),names_here(gl),'Orientation','Vertical','Location','SouthEast')
        case 5
            names_here = {'CCMVal','ERA-40','Geodetic'};
        case 6
            names_here = {'Unconstrained Model AAM','ERA-40 AAM','Inferred from Obs.'};
        case default
            names_here = names_total(R2);
            legend(lh(gl),names_here(gl),'Orientation','Vertical','Location','SouthEast')
    end
    legend('boxoff')
    set(gca,'YTick',yticks)
    % add plot labels:
    if plot_type ==4
        text(xlim(1),yy(2),plot_letter,'HorizontalAlignment','Right','FontSize',FS)
    end
    if plot_type == 1
        text(xlim(1),yy(2)-0.1*dylim,plot_letter,'HorizontalAlignment','Left','FontSize',FS)
    end        
  
   
end



