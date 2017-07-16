function plot_aam_mm_slice(comp,lonplot)
% Plot the differences between monthly contributions to AAM terms between 
% different data sets or runs (here: EMAC runs with different resolutions)
% The motivation is that hemispheric asymmetry is really what makes the various terms wobble.
%
% This code makes 3 plots: for the m term, the w term, and the total.  
%
% 1 Nov 2010 - see Notes vol 4, p. 2  
% Mods:
%  14 Nov 2010
%
% Inputs:
% comp: which vector component of AAM to plot.
% longplot: set to 1 to get a plot over longitudes, 0 to get over latitudes
%-------------------------------------------------------------------------------------------


% read in lat / lon slices from previously stored AAM mat files.
last20 = 0;	% flag for last 20 years.
[lat42, lon42, dum42] = aam_mm_slice(comp,'t7_T42L39',last20);
[lat31, lon31, dum31] = aam_mm_slice(comp,'ref2_T31L39',last20);
[lat63, lon63, dum63] = aam_mm_slice(comp,'ref_T63L39',last20);
[latEI, lonEI, dumEI] = aam_mm_slice(comp,'ERAinterim',last20);
[latCC, lonCC, dumCC] = aam_mm_slice(comp,'CCMval',last20);

names = {'T31L39','T42L39', 'T63L39','T42L90','ERAint'};

X42 = [dum42;sum(dum42,1)];
X31 = [dum31;sum(dum31,1)];
X63 = [dum63;sum(dum63,1)];
XEI = [dumEI;sum(dumEI,1)];
XCC = [dumCC;sum(dumCC,1)];

% shift lons for ERA
lonEI = lonEI-180;

% these are the months to plot.
months = [1,7];
nm = length(months);

% also compute jan - july differences
D42 = X42(:,1,:,:) - X42(:,7,:,:) ;
D31 = X31(:,1,:,:) - X31(:,7,:,:) ;
D63 = X63(:,1,:,:) - X63(:,7,:,:) ;
DEI = XEI(:,1,:,:) - XEI(:,7,:,:) ;
DCC = XCC(:,1,:,:) - XCC(:,7,:,:) ;

% --for each month make a plot comparing latitudinal profiles

% define plot settings
col = [colormap(summer(4))];
col31 = col(1,:);
col42 = col(2,:);
col63 = col(3,:);
colCC = col(4,:);
colEI = zeros(1,3);

LW = 3;
LS = {'-','-.'};
lh = zeros(1,3);
TTw = {'\chi_1 Wind Term','\chi_2 Wind Term','\Delta LOD Wind Term'};
TTm = {'\chi_1 Mass Term','\chi_2 Mass Term','\Delta LOD Mass Term'};
TTt = {'\chi_1 Total','\chi_2 Total','\Delta LOD Total'};
TT = [TTw; TTm; TTt];
%TT = {'\chi_1 AAM Contributions','\chi_2 AAM Contributions','\Delta LOD AAM Contributions'};
pl = {'Wind Term','Mass Term','Total'};

ax = zeros(2,3,2);
ax(1,1,:) = 80*[-1,1];    % range for X1w in mas
ax(1,2,:) = 80*[-1,1];    % range for X2w in mas
ax(1,3,:) = 0.4*[-1,1];       % range for DLODw in ms
ax(2,1,:) = 10*[-1,1];    % range for X1m in mas
ax(2,2,:) = 10*[-1,1];    % range for X2m in mas
ax(2,3,:) = 0.05*[-1,1];       % range for DLODm in ms

line1 = -180:1:180;
line2 = -90:90;

% longitude profiles: mass, wind, total
for ii = 1:2
  figure(ii),clf
  plot(line1,line1*0,'k')
  hold on
  for im = 1:2
    m = months(im);
    if lonplot 
      lh(1,im) = plot(lon31, squeeze(mean(X31(ii,m,:,:),3)),'Color',col31,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(2,im) = plot(lon42, squeeze(mean(X42(ii,m,:,:),3)),'Color',col42,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(3,im) = plot(lon63, squeeze(mean(X63(ii,m,:,:),3)),'Color',col63,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(4,im) = plot(lonCC, squeeze(mean(XCC(ii,m,:,:),3)),'Color',colCC,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(5,im) = plot(lonEI, squeeze(mean(XEI(ii,m,:,:),3)),'Color',colEI,'LineWidth',LW,'LineStyle',char(LS(im)));
      axis([-180 180 ax(ii,comp,1) ax(ii,comp,2)])
      %text(-170,0.8*ax(ii,comp,2),TT(ii,comp))
      xlabel('Longitude')
      if ii < 3, ylabel('mas'), else ylabel('ms'), end
   else
      lh(1,im) = plot(squeeze(mean(X31(ii,m,:,:),4)),lat31,'Color',col31,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(2,im) = plot(squeeze(mean(X42(ii,m,:,:),4)),lat42,'Color',col42,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(3,im) = plot(squeeze(mean(X63(ii,m,:,:),4)),lat63,'Color',col63,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(4,im) = plot(squeeze(mean(XCC(ii,m,:,:),4)),latCC,'Color',colCC,'LineWidth',LW,'LineStyle',char(LS(im)));
      lh(5,im) = plot(squeeze(mean(XEI(ii,m,:,:),4)),latEI,'Color',colEI,'LineWidth',LW,'LineStyle',char(LS(im)));
      axis([ax(ii,comp,1) ax(ii,comp,2) -90 90])
      %text(.8*ax(ii,comp,1),75,pl(ii))
      ylabel('Latitude')
      if ii < 3, xlabel('mas'), else xlabel('ms'), end
   end
  end
  legend(lh(:,1), names,'Orientation','Horizontal','Location','SouthOutside')
  legend('BoxOff')
  title(TT(ii,comp));
end


% export plots here:

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

if lonplot 
  fig_name_w = [plot_dir,'lonslice_X',num2str(comp),'w.png']
  fig_name_m = [plot_dir,'lonslice_X',num2str(comp),'m.png']
else
  fig_name_w = [plot_dir,'latslice_X',num2str(comp),'w.png']
  fig_name_m = [plot_dir,'latslice_X',num2str(comp),'m.png']
end


LW = 4;
ph = 10;        % paper height
pw = 10;        % paper width
fs = 30;        % fontsize


exportfig(1,fig_name_w,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
exportfig(2,fig_name_m,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');


