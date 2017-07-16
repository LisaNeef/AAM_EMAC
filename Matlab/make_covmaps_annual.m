clear all;
runid = 'CCMval';
fil_order = 2;
abs_val =  1;
tscale = 2;

figure(1),clf

h = zeros(1,6);

PL = {'(a) \chi_1 Mass Term', '(b) \chi_2 Mass Term','(c) \chi_3 Wind Term', '(d) \chi_3 Mass Term'};

jj = 1;
terms = [2,2,1,2];
comp = [1,2,3,3];
extras = zeros(1,4);        % plot extras: set to 1 to get colorbar, etc
extras(3:4) = 1;

% make the plots

for jj = 1:4
        
        plot_label = PL(jj);
        h(jj) = subplot(2,2,jj);
        plot_covariance_maps_bootstrap(runid,tscale,comp(jj),terms(jj),abs_val,fil_order,plot_label,0,extras(jj))
  
end


% make them look good.



set(gcf,'Position',[911 551 845 719])

x0 = 0.09;                  % left position
dy = 0;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .03;
w2 = (w-1*dw-x0)/2;               % width per figure
ht = (y0-.05)/2;          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])
set(h(2),'Position',[x0+w2+dw y3 w2 ht])

y3 = y0-dy-2*ht;
set(h(3),'Position',[x0 y3 w2 ht])
set(h(4),'Position',[x0+w2+dw y3 w2 ht])


    for k = 1:4
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
    print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig11_covariance.eps

    % less shitty gridlines

    fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig11_covariance.eps')




