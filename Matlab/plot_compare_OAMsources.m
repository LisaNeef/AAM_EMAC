%  plot_compare_OAMsouruces
%
%  Make a plot comparing different available estimates for OAM.
%  31 Jan 2011
%  see notes vol. 4, p. 82
%
%
%
%---------------------------------------------------------------



% load the different runs

[Xw_ECCO,Xm_ECCO,mjd_ECCO] = read_EFs('oam','ECCO',0);
[Xw_ERA40,Xm_ERA40,mjd_ERA40] = read_EFs('oam','ERA40',0);
[Xw_ERAint,Xm_ERAint,mjd_ERAint] = read_EFs('oam','ERAinterim',0);
[Xw_ECMWF,Xm_ECMWF,mjd_ECMWF] = read_EFs('oam','ECMWF',0);

aam_constants_gross

% remove long term mean form all
Xw_ECCO2 = detrend(Xw_ECCO,'constant'); 	Xm_ECCO2 = detrend(Xm_ECCO,'constant');
Xw_ERA402 = detrend(Xw_ERA40,'constant'); 	Xm_ERA402 = detrend(Xm_ERA40,'constant');
Xw_ERAint2 = detrend(Xw_ERAint,'constant'); 	Xm_ERAint2 = detrend(Xm_ERAint,'constant');
Xw_ECMWF2 = detrend(Xw_ECMWF,'constant'); 	Xm_ECMWF2 = detrend(Xm_ECMWF,'constant');

%---plots!

random_colors

LW = 2;
col_ERA40 = rand(1,3);
col_ERAint= rand(1,3);
col_ECMWF = rand(1,3);
col_ECCO= acid;

%---wind terms

for ii = 1:3
  figure(ii),clf
  lh = zeros(4,1);
  if ii < 3, c = rad2mas; else, c = LOD0_ms; end
  lh(1) = plot(mjd_ERA40,c*Xw_ERA402(ii,:),'Color',col_ERA40,'LineWidth',2);
  hold on
  lh(2) = plot(mjd_ERAint,c*Xw_ERAint2(ii,:),'Color',col_ERAint,'LineWidth',2);
  lh(3) = plot(mjd_ECMWF,c*Xw_ECMWF2(ii,:),'Color',col_ECMWF,'LineWidth',2);
  lh(4) = plot(mjd_ECCO,c*Xw_ECCO2(ii,:),'Color',col_ECCO,'LineWidth',2);
  if ii == 1, title('OAM Comparison: \chi_1 Current Term'); end
  if ii == 2, title('OAM Comparison: \chi_2 Current Term'); end
  if ii == 3, title('OAM Comparison: \chi_3 Current Term'); end
  xlabel('MJD')
  ylabel('mas')
  legend(lh,{'OMCT ERA40','OMCT ERA-Interim','OMCT-ECMWF','ECCO'},'Orientation','Horizontal','Location','South')
end


%---mass terms

for ii = 1:3
  figure(ii+3),clf
  lh = zeros(4,1);
  lh(1) = plot(mjd_ERA40,rad2mas*Xm_ERA402(ii,:),'Color',col_ERA40,'LineWidth',2);
  hold on
  lh(2) = plot(mjd_ERAint,rad2mas*Xm_ERAint2(ii,:),'Color',col_ERAint,'LineWidth',2);
  lh(3) = plot(mjd_ECMWF,rad2mas*Xm_ECMWF2(ii,:),'Color',col_ECMWF,'LineWidth',2);
  lh(4) = plot(mjd_ECCO,rad2mas*Xm_ECCO2(ii,:),'Color',col_ECCO,'LineWidth',2);
  if ii == 1, title('OAM Comparison: \chi_1 Mass Term'); end
  if ii == 2, title('OAM Comparison: \chi_2 Mass Term'); end
  if ii == 3, title('OAM Comparison: \chi_3 Mass Term'); end
  xlabel('MJD')
  ylabel('mas')
  legend(lh,{'OMCT ERA40','OMCT ERA-Interim','OMCT-ECMWF','ECCO'},'Orientation','Horizontal','Location','South')
end

%---export plots!

LW = 4;
ph = 5;        % paper height
pw = 15;        % paper width
FS = 20;        % fontsize

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

exportfig(1,[plot_dir,'comp_OAMsources_X1w.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
exportfig(2,[plot_dir,'comp_OAMsources_X2w.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
exportfig(3,[plot_dir,'comp_OAMsources_X3w.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

exportfig(4,[plot_dir,'comp_OAMsources_X1m.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
exportfig(5,[plot_dir,'comp_OAMsources_X2m.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
exportfig(6,[plot_dir,'comp_OAMsources_X3m.png'],'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');



