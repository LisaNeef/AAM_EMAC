% figure 9 for the aam paper

% initial settings

figure(1),clf
clear all;


%make the plot

comp = 3;
apply_lag = 0;
export = 0;

[h1,h2] = plot_long_term(comp, apply_lag,export);
h(1) = h1;
h(2) = h2;

% make it not suck

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
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );
end

% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig9.eps



% less shitty lines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig9.eps')


