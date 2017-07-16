% make a plot examining the correlations between 4 different CAM datasets and the observed ERPs.
% 14 Dec 2010
% notes vol4, p. 50
% Mods
%  17 Dec 2010 ...fix conversion from days on time axis to years

%---plot settings

LW = 3;

T1 = 'Correlation Between CAM Timeseries and Observed ERPs: \chi_1';
T2 = 'Correlation Between CAM Timeseries and Observed ERPs: \chi_2';
T3 = 'Correlation Between CAM Timeseries and Observed ERPs: \chi_3';
TT = {T1;T2;T3};


%---cycle components, retrieve data, make plot
exp_type = 2;
tscale = 4;
for comp = 1:3

  [mc,lag,cc,runs,tt,C,L] =  correlations_aam_filtered_v2(comp,tscale,exp_type);

  figure(comp),clf
    nR = length(runs);
    col = colormap(hsv(nR));
    lh = zeros(1,nR);
    for ii = 1:nR
      lh(ii) = plot(L(ii,:)/365.5,C(ii,:),'LineWidth',LW,'Color',col(ii,:));
      hold on
    end
    legend(lh,runs)
    axis([-40,40,-0.7,0.7]);
    xlabel('Lag (Years)')
    ylabel('Correlation')
    title(TT(comp))
end



%---export plots

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
rend = 'opengl';
LW = 3;
ph = 8;        % paper height
pw = 20;        % paper width
fs = 26;        % fontsize

for ii = 1:3
  fig_name = [plot_dir,'corr_ERPs_CAM_X',num2str(ii),'.png'];
  exportfig(ii,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');
end


