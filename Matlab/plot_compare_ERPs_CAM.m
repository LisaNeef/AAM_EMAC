% plot_compare_ERPs_CAM.m
%
% read in the ERPs and compare the correspondig core-mantle excitation functions
% (of which there are 4 sets)
% to get an idea of what they look like.
%
%  29 Nov 2010 - Notes vol. 4, p. 33
%---------------------------------------------------------------------------

% get the ERP timeseries.

[x1,x2,dlod,mjd_geo,ex1,ex2,edlod] = read_eops('IERS',1);
GEO_trend = [x1,x2,dlod]';

% take out its long-term mean.
GEO = GEO_trend*0;
for ii = 1:3
  GEO(ii,:) = detrend(GEO_trend(ii,:),'constant');
end 


% read in the four CAM series

CAM = zeros(4,3,457);

[x1_cam,x2_cam,dlod_cam,mjd_cam,date_cam] = read_CAM_JH('t',1);
CAM(1,:,:) = [x1_cam,x2_cam,dlod_cam]';

[x1_cam,x2_cam,dlod_cam,mjd_cam,date_cam] = read_CAM_JH('t',4);
CAM(2,:,:) = [x1_cam,x2_cam,dlod_cam]';

[x1_cam,x2_cam,dlod_cam,mjd_cam,date_cam] = read_CAM_JH('s',1);
CAM(3,:,:) = [x1_cam,x2_cam,dlod_cam]';

[x1_cam,x2_cam,dlod_cam,mjd_cam,date_cam] = read_CAM_JH('s',4);
CAM(4,:,:) = [x1_cam,x2_cam,dlod_cam]';

names = {'Obs.','T1','T4','S1','S4'};
TT = {'\chi_1','\chi_2','\Delta LOD'};


% turn the mjd series into matlab datenum series
  tgeo = mjd_geo*0+NaN;
  tcam = mjd_cam*0+NaN;
  [y, m, d] = mjd2date(mjd_geo);
  for ii=1:length(tgeo),  tgeo(ii)=datenum([y(ii) m(ii) d(ii)]); end
  [y, m, d] = mjd2date(mjd_cam);
  for ii=1:length(tcam),  tcam(ii)=datenum([y(ii) m(ii) d(ii)]); end


% plot it!

LW = 3;
col = [colormap(jet(4))];

for ii = 1:3  
  lh = zeros(1,5);
  figure(ii),clf
    lh(1) = plot(tgeo,GEO(ii,:),'k','LineWidth',LW);
    hold on
    for icam = 1:4, lh(icam+1) = plot(tcam,squeeze(CAM(icam,ii,:)),'Color',col(icam,:),'LineWidth',LW); end
    legend(lh,names,'Location','South','Orientation','Horizontal')
    title([char(TT(ii)),' Comparison of Observed and CAM'])
    xlabel('year')
    if ii < 3, ylabel('mas'), else ylabel('ms'), end
    datetick('x','yyyy')
end

% export it!

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
rend = 'opengl';
LW = 3;
ph = 8;        % paper height
pw = 20;        % paper width
fs = 26;        % fontsize

for ii = 1:3
  fig_name = [plot_dir,'comp_ERPs_CAM_X',num2str(ii),'.png'];
  exportfig(ii,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');
end






