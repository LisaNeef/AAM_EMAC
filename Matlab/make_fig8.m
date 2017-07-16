% make plots of the correlations between models and the geodetic residual

figure(1),clf
clear all;
plot_type = 1;
alpha = 0.05;       % 95% confidence limit.            
h = [0,0];

tscale = [3,2];

for ii = 1:2
    h(ii) = subplot(2,1,ii);
    plot_correlations_bar_bootstrap(tscale(ii),alpha,plot_type,0)
end

set(gca,'LineWidth',1)

% make it look less shitty

set(gcf,'Position',[10 30 1200 2*463])

x0 = 0.09;                  % left position
dy = .05;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .08;
w2 = (w-2*dw-x0);               % width per figure
ht = (y0-3*dy)/2;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])

y3 = y0-2*dy-2*ht;
set(h(2),'Position',[x0 y3 w2 ht])

ytick = -0.5:.25:1;

for k = 1:2
    if k ==1
        set(h(k),'XTick',[])
    end
    set( h(k)                       , ...
      'FontName'   , 'Helvetica' ,...
      'FontSize'   , 16         ,...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'YTick'       , ytick,...
      'YTickLabel'  ,  {'-0.5','-0.25','0','0.25','0.5','0.75','1.0'},...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );
end

% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig8.eps



% less shitty lines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig8.eps')


