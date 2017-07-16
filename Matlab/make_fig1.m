% code and export settings Fig. 1 of my aam paper.
% last update: 12 Sept 2011./

clear all;

plot_type = 4;
spect_type = 1;
BW = 0;
terms = [1,2];

PL = {'(a) \chi_1'; '(b) \chi_2'; '(c) \chi_3'};

h = zeros(1,3);
annotate = [1,0,0];
label_x_axis = [0,0,1];

figure(1),clf
for ii = 1:3
    
    h(ii) = subplot(3,1,ii);
    plot_aam_spectra(ii,terms,plot_type,spect_type,BW,PL(ii),annotate(ii),0,0,label_x_axis(ii),1);
        
end


% now go through panesl and make them look right

set(gcf,'Position',[911 551 572 719])

x0 = 0.09;                  % left position
dy = .045;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .05;
w2 = w-dw-x0;               % width per figure
ht = (y0-4*dy)/3;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])
set(h(1),'XTick',[])

y3 = y0-2*dy-2*ht;
set(h(2),'Position',[x0 y3 w2 ht])
set(h(1),'XTick',[])

y3 = y0-3*dy-3*ht;
set(h(3),'Position',[x0 y3 w2 ht])


% export to eps

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig1.eps

% less shitty gridlines

fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig1.eps')
