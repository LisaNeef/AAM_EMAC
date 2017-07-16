function plot_aam_acvar(comp,term,lonplot,var_to_plot,textlabel)
% plot the annual variance in regional contributions to EAFs
% The goal is to see how different regions affect the annual cycles in EAFs
% This code makes 3 plots: for the m term, the w term, and the total.  
%
% 17 Jan 2011
%
% Inputs:
%   comp: which vector component of AAM to plot.
%   term: 1 for wind, 2 for mass, 3 for total.
%   lonplot: set to 1 to get a plot over longitudes, 0 to get over latitudes
%   var_to_plot: put 'X' for AAM, 'u' for u300, or 'p' for ps
%   text: the text with which to annotate the plot.
%-------------------------------------------------------------------------------------------

if var_to_plot == 'u', comp = 1; end
if var_to_plot == 'p', comp = 1; end

% read in lat / lon slices from previously stored AAM mat files.
[lat31, lon31, varX31, varU31, varPS31] = aam_acvar_slice(comp,'ref2_T31L39');
[lat42, lon42, varX42, varU42, varPS42] = aam_acvar_slice(comp,'t7_T42L39');
[lat63, lon63, varX63, varU63, varPS63] = aam_acvar_slice(comp,'ref_T63L39');
[latCC, lonCC, varXCC, varUCC, varPSCC] = aam_acvar_slice(comp,'CCMval');
[latEI, lonEI, varXEI, varUEI, varPSEI] = aam_acvar_slice(comp,'ERAint');

names = {'EMAC31','EMAC42', 'EMAC63','CCMVal','ERA-Int'};

% average over all the years of each run, depending on variable.

if var_to_plot == 'X'
  var31 =  (nanmean(varX31,2));
  var42 =  (nanmean(varX42,2));
  var63 =  (nanmean(varX63,2));
  varCC =  (nanmean(varXCC,2));
  varEI =  (nanmean(varXEI,2));
end

if var_to_plot == 'u'
  var31 =  (nanmean(varU31,2));
  var42 =  (nanmean(varU42,2));
  var63 =  (nanmean(varU63,2));
  varCC =  (nanmean(varUCC,2));
  varEI =  (nanmean(varUEI,2));
end

if var_to_plot == 'p'
  var31 =  (nanmean(varPS31,2));
  var42 =  (nanmean(varPS42,2));
  var63 =  (nanmean(varPS63,2));
  varCC =  (nanmean(varPSCC,2));
  varEI =  (nanmean(varPSEI,2));
end


% shift lons for ERA
%lonEI = lonEI-180;

% add the weighting functions (for u) to plots of  AAM in lat/lon belts.
d2r = pi/180;
rlat63 = d2r*lat63;
rlon63 = d2r*lon63;

Wlon = zeros(3,length(lon63));
Wlat = zeros(3,length(lat63));

Wlon(1,:) = cos(rlon63); 		% X1 u-weighting with constant lat
Wlon(2,:)= sin(rlon63); 		% X2 weighting with constant lat
Wlon(3,:) = zeros(size(rlon63));
Wlat(1,:) = sin(rlat63).*cos(rlat63);	% X1 weighting with constant lon
Wlat(2,:) = sin(rlat63).*cos(rlat63);	% X2 weighting with constant lon
Wlat(3,:) = cos(rlat63).^2;		% X3 weighting

% --plots: lat profiles of variances

% define plot settings
col = aam_paper_colors;

  % pointers to where colors are:
  e31 = 1;
  e42 = 2;
  e63 = 3;
  cc =  4;
  ei =  6;

LW = 3;
LS = '-';
LSweights = {'--','--','--'};
lh = zeros(1,5);

ax = zeros(3,3,2);
if var_to_plot == 'X'
  ax(1,1,:) = 25000*[0,1];    % range for X1t in mas^2
  ax(1,2,:) = 25000*[0,1];    % range for X2t in mas^2
  ax(1,3,:) = 0.12*[0,1];       % range for DLODt in ms^2
  ax(2,1,:) = 25000*[0,1];    % range for X1w in mas^2
  ax(2,2,:) = 25000*[0,1];    % range for X2w in mas^2
  ax(2,3,:) = 0.12*[0,1];       % range for DLODw in ms^2
  ax(3,1,:) = 300*[0,1];    % range for X1m in mas^2
  ax(3,2,:) = 300*[0,1];    % range for X2m in mas^2
  ax(3,3,:) = 0.002*[0,1];       % range for DLODm in ms^2
end
if strcmp(var_to_plot,'u')
   ax(1,1,:) = 250*[0,1];	% range for u in (m/s)^2
end
if strcmp(var_to_plot,'ps')
   ax(1,1,:) = 250*[0,1];	% range for ps in (Pa)^2
end
if lonplot, ax = ax*0.6; end


% longitude profiles: mass, wind, total
ii = term;
  figure(ii),clf
  hold on
    if lonplot 
      if var_to_plot == 'X'
        plot(lon63,0.3*ax(ii,comp,2)*(1+Wlon(comp,:)),'LineStyle',char(LSweights(comp)),'Color',0.5*ones(1,3),'LineWidth',5)
      end
      lh(1) = plot(lon31, squeeze(nanmean(var31(ii,:,:,:),3)),'Color',col(e31,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(2) = plot(lon42, squeeze(nanmean(var42(ii,:,:,:),3)),'Color',col(e42,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(3) = plot(lon63, squeeze(nanmean(var63(ii,:,:,:),3)),'Color',col(e63,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(4) = plot(lonCC, squeeze(nanmean(varCC(ii,:,:,:),3)),'Color',col(cc,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(5) = plot(lonEI, squeeze(nanmean(varEI(ii,:,:,:),3)),'Color',col(ei,:),'LineWidth',LW,'LineStyle',char(LS));
      %axis([-180,180,ax(ii,comp,1), ax(ii,comp,2)])
      xlabel('Longitude')
      if comp < 3, ylabel('mas^2'), else ylabel('ms^2'), end
      if var_to_plot == 'u', ylabel('(m/s)^2'); end
      if var_to_plot == 'p', ylabel('Pa^2'); end
      legend(lh, names,'Orientation','Vertical','Location','NorthEast')
      set(gca,'XTick',-180:45:180);
      set(gca,'XLim',[-180 180]);
   else
      if var_to_plot == 'X'
        plot(0.3*ax(ii,comp,2)*(1+Wlat(comp,:)),lat63,'LineStyle',char(LSweights(comp)),'Color',0.5*ones(1,3),'LineWidth',5)
      end
      lh(1) = plot(squeeze(nanmean(var31(ii,:,:,:),4)),lat31,'Color',col(e31,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(2) = plot(squeeze(nanmean(var42(ii,:,:,:),4)),lat42,'Color',col(e42,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(3) = plot(squeeze(nanmean(var63(ii,:,:,:),4)),lat63,'Color',col(e63,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(4) = plot(squeeze(nanmean(varCC(ii,:,:,:),4)),latCC,'Color',col(cc,:),'LineWidth',LW,'LineStyle',char(LS));
      lh(5) = plot(squeeze(nanmean(varEI(ii,:,:,:),4)),latEI,'Color',col(ei,:),'LineWidth',LW,'LineStyle',char(LS));
      ylabel('Latitude')
      if comp < 3, xlabel('mas^2'), else xlabel('ms^2'), end
      if var_to_plot == 'u', xlabel('(m/s)^2'); end
      if var_to_plot == 'p', xlabel('Pa^2'); end
      legend(lh, names,'Orientation','Vertical','Location','SouthEast')
      set(gca,'YTick',-90:20:90);
      set(gca,'YLim',[-90 90]);
    end
  xx = get(gca,'XLim');
  dx = xx(2)-xx(1);
  yy = get(gca,'Ylim');
  dy = yy(2)-yy(1);
  if lonplot
      text(xx(1)+0.01*dx,yy(1)+0.95*dy,textlabel,'HorizontalAlignment','left')
  else
      text(xx(2),yy(1)+0.9*dy,textlabel,'HorizontalAlignment','right')
  end
  legend('BoxOff')
  grid on
% export plots here:

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

if lonplot 
  fig_name_w = [plot_dir,'lonslice_acvar_',var_to_plot,num2str(comp),'w.png'];
  fig_name_m = [plot_dir,'lonslice_acvar_',var_to_plot,num2str(comp),'m.png'];
  fig_name_t = [plot_dir,'lonslice_acvar_',var_to_plot,num2str(comp),'t.png'];
else
  fig_name_w = [plot_dir,'latslice_acvar_',var_to_plot,num2str(comp),'w.png'];
  fig_name_m = [plot_dir,'latslice_acvar_',var_to_plot,num2str(comp),'m.png'];
  fig_name_t = [plot_dir,'latslice_acvar_',var_to_plot,num2str(comp),'t.png'];
end


LW = 4;
ph = 10;        % paper height
pw = 10;        % paper width
fs = 30;        % fontsize

switch term
    case 1
        fig_name = fig_name_w;
    case 2
        fig_name = fig_name_m;
    case 3
        fig_name = fig_name_t;
end
exportfig(term,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');


