
%AAM Paper plot comparing annual cycle in 3 AEFs between ERA and various
%EMAC simulations.


figure(1),clf
comp_geo = 0;
extra_labels = zeros(3,2);
extra_labels(1,2) = 1;
label_x_axis = zeros(3,2);
label_x_axis(3,:) = [1,1];
label_y_axis = zeros(3,2);
label_y_axis(:,2) = ones(3,1);
h = zeros(1,6);

PLm = {'(a) \chi_1 Mass Term', '(c) \chi_2 Mass Term','(e) \chi_3 Mass Term'};
PLw = {'(b) \chi_1 Wind Term', '(d) \chi_2 Wind Term','(f) \chi_3 Wind Term'};

jj = 1;

for ii = 1:3
    for term = [2,1]
        switch term
            case 1
                plot_label = PLw(ii);
            case 2
                plot_label = PLm(ii);
        end
        h(jj) = subplot(3,2,jj);
        plot_aam_ann_cycle_compare(ii,term,comp_geo,plot_label,0,label_x_axis(ii,term),label_y_axis(ii,term));
        jj = jj+1; 

    end
end


% now go through panesl and make them look right

set(gcf,'Position',[911 551 845 719])

x0 = 0.09;                  % left position
dy = .05;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .05;
w2 = (w-3*dw-x0)/2;               % width per figure
ht = (y0-4*dy)/3;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])
set(h(2),'Position',[x0+w2+dw y3 w2 ht])

y3 = y0-2*dy-2*ht;
set(h(3),'Position',[x0 y3 w2 ht])
set(h(4),'Position',[x0+w2+dw y3 w2 ht])

y3 = y0-3*dy-3*ht;
set(h(5),'Position',[x0 y3 w2 ht])
set(h(6),'Position',[x0+w2+dw y3 w2 ht])


for k = 1:6
    if k < 5
        set(h(k),'XTick',[])
    end
    set( h(k)                       , ...
      'FontName'   , 'Helvetica' ,...
      'FontSize'   , 12 ,...
     'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );
end


% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig6.eps

% less shitty gridlines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig6.eps')



