clear all;
runid = 'CCMval';
fil_order = 2;
abs_val =  1;
tscale = 1;

figure(1),clf

h = zeros(1,3);

PL = {'(a) \chi_1 Mass Term', '(b) \chi_2 Mass Term','(c) \chi_3 Wind Term'};

jj = 1;
terms = [2,2,1];
comp = [1,2,3];
extras = zeros(1,3);        % plot extras: set to 1 to get colorbar, etc
extras(3) = 1;

% make the plots

for jj = 1:3
        
        plot_label = PL(jj);
        h(jj) = subplot(3,1,jj);
        plot_covariance_maps_bootstrap(runid,tscale,comp(jj),terms(jj),abs_val,fil_order,plot_label,0,extras(jj))
  
end


% make them look good.



set(gcf,'Position',[0.545 0.05 0.355 0.47])

x0 = 0.09;                  % left size margin
dy = 0;
y0 = 0.96;                  % top of plot
w = 1; 
dw = .1;                   % space on right side of plots
w2 = (w-dw-x0);               % width per figure
dy = 0.05;
ht = (y0-2*dy)/3-.01;          % height per figure

y3 = y0-ht;
set(h(1),'Position',[x0 y3 w2 ht])

y3 = y0-2*ht-dy;
set(h(2),'Position',[x0 y3 w2 ht])

y3 = y0-2*dy-3*ht;
set(h(3),'Position',[x0 y3 w2 ht])


    for k = 1:3
        if k < 5
            set(h(k),'XTick',[])
        end
        set( h(k)                       , ...
          'FontName'   , 'Helvetica' ,...
          'FontSize'    ,12         , ...
          'Box'         , 'off'     , ...
          'TickDir'     , 'out'     , ...
          'TickLength'  , [.02 .02] , ...
          'XColor'      , [.3 .3 .3], ...
          'YColor'      , [.3 .3 .3], ...
          'LineWidth'   , 1         );
    end


    % export to eps

    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig10.eps

    % less shitty gridlines

    fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig10.eps')




