% commands to make the annual cycle plots for the AAM diagnosis paper.

%---plots with comparison to obs.
figure(1),clf
comp_geo = 1;
label_x_axis = [0,0,1];
h = zeros(1,3);
for ii = 1:3
    switch ii
        case 1
            plot_label = '(a) \chi_1';
        case 2
            plot_label = '(b) \chi_2';
        case 3
            plot_label = '(c) \chi_3';
    end
    
    h(ii) = subplot(3,1,ii);
    plot_aam_ann_cycle_compare(ii,3,comp_geo,plot_label,0,label_x_axis(ii),1);
    
end


% now go through panesl and make them look right

set(gcf,'Position',[911 551 400 800])

x0 = 0.09;                  % left position
dy = .05;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .05;
w2 = w-dw-x0;               % width per figure
ht = (1-4*dy)/3;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])
set(h(1),'XTick',[])


y3 = y0-2*dy-2*ht;
set(h(2),'Position',[x0 y3 w2 ht])
set(h(2),'XTick',[])

y3 = y0-3*dy-3*ht;
set(h(3),'Position',[x0 y3 w2 ht])

for k = 1:3
set( h(k)                       , ...
  'FontName'   , 'Helvetica' ,...
  'FontSize'   , 16 ,...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );
end

% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig5.eps

% less shitty gridlines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig5.eps')



