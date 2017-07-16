% code and export settings Fig. 1 of my aam paper.
% last update: 12 Sept 2011./

clear all;

plot_type = 1;
spect_type = 1;
BW = 0;
terms = [1,2];

PLm = {'(a) \chi_1 Mass'; '(c) \chi_2 Mass'; '(e) \chi_3 Mass'};
PLw = {'(b) \chi_1 Wind'; '(d) \chi_2 Wind'; '(f) \chi_3 Wind'};


h = zeros(1,6);
extra_labels = zeros(3,2);
extra_labels(1,2) = 1;
label_x_axis = zeros(3,2);
label_x_axis(3,:) = [1,1];
label_y_axis = zeros(3,2);
label_y_axis(:,2) = ones(3,1);

figure(1),clf
jj = 1;
for ii = 1:3
    for term = [2,1]
      switch term
          case 1
              label = PLw(ii);
          case 2
              label = PLm(ii);
      end
      h(jj) = subplot(3,2,jj);
      plot_aam_spectra(ii,term,plot_type,spect_type,BW,label,0,0,extra_labels(ii,term),label_x_axis(ii,term),label_y_axis(ii,term));
      grid on
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
      'FontSize'   , 16         ,...
      'Box'         , 'off'     , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'LineWidth'   , 1         );
end


% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig2.eps

% less shitty gridlines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig2.eps')
