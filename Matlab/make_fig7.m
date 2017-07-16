% figure 3 for the AAM paper

clear all;
figure(1),clf
plot_type = 2;
terms = 3;
error_case = 1;
tscale = 3;


h = zeros(1,2);

% X2 plot for the paper
plot_label = '(a) \chi_2';
add_SOI = 0;
comp = 2;
h(1) = subplot(2,1,1);
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label,0,0)


% X3 plot for the paper
plot_label = '(b) \chi_3';
comp = 3;
add_SOI = 1;
subplot(2,1,2);
[hdum1,hdum2,hlegend]=  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label,0,1);
h(2) = hdum1;
h(3) = hdum2;
  
% make the spaces and stuff look good

set(gcf,'Position',[911 551 999 659])

x0 = 0.09;                  % left position
dy = .05;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .08;
w2 = (w-dw-x0);               % width per figure
ht = (y0-3*dy)/2;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])

y3 = y0-2*dy-2*ht;
set(h(2),'Position',[x0 y3 w2 ht])
set(h(3),'Position',[x0 y3 w2 ht])


for k = 1:3
    if k ==1
        set(h(k),'XTick',[])
    end
    set( h(k)                       , ...
      'FontName'   , 'Helvetica' ,...
            'FontSize'   , 16 ,...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.01 .01] , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );
end

% also turn off the box on the SOI plot's legend.
set(hlegend,'Box','off')

% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig7.eps

% less shitty lines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig7.eps')


