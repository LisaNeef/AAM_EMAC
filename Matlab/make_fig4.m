% make plots of the correlations between models and the geodetic residual

figure(1),clf
clear all;
plot_type = 1;
alpha = 0.05;       % 95% confidence limit.            


  plot_correlations_bar_bootstrap(1,alpha,plot_type,0)

set(gca,'LineWidth',1)

% make it look less shitty

set(gcf,'Position',[461 508 1200 463])
ytick = -0.5:0.25:1.0;


    set( gca                       , ...
      'FontName'   , 'Helvetica' ,...
      'FontSize'   , 16         , ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'YTickLabel'  ,  {'-0.5','-0.25','0','0.25','0.5','0.75','1.0'},...
      'YTick'       , ytick,...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );


% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig4.eps



% less shitty lines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig4.eps')
