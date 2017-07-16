% quick plots of the CAM timeseries.
% 25 nov 2010


[X1,X2,dlod_s1,mjd,date] = read_CAM_JH('s',1);
[X1,X2,dlod_t1,mjd,date] = read_CAM_JH('t',1);
[X1,X2,dlod_s4,mjd,date] = read_CAM_JH('s',4);
[X1,X2,dlod_t4,mjd,date] = read_CAM_JH('t',4);

lh = zeros(1,4);

figure(1),clf
  lh(1) = plot(date,dlod_s1,'Color',rand(3,1),'LineWidth',3);
  hold on
  lh(2) = plot(date,dlod_t1,'Color',rand(3,1),'LineWidth',3);
  lh(3) = plot(date,dlod_s4,'Color',rand(3,1),'LineWidth',3);
  lh(4) = plot(date,dlod_t4,'Color',rand(3,1),'LineWidth',3);

legend(lh,'S1','T1','S4','T4')
title('\Delta LOD Due to CAM Variations')
xlabel('date')
ylabel('ms')

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
fig_name = [plot_dir,'CAM_JH.png'];

LW = 4;
ph = 10;        % paper height
pw = 10;        % paper width
fs = 26;        % fontsize

exportfig(1,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

