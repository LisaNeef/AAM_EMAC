clear all;
runid = 'CCMval';
fil_order = 2;
abs_val =  1;
tscale = 3;

figure(1),clf

h = zeros(1,2);

PL = {'(a) \chi_2 Mass Term','(b) \chi_3 Wind Term'};
terms = [2,1];
comp = [2,3];


jj = 1;

extras = zeros(1,5);        % plot extras: set to 1 to get colorbar, etc


% make the plots

for jj = 1:2    
        h(jj) = subplot(1,2,jj);
        plot_covariance_maps_bootstrap(runid,tscale,comp(jj),terms(jj),abs_val,fil_order,PL(jj),0,extras(jj))
        jj = jj+1; 
end



% make them look good.



set(gcf,'Position',[834 595 845 378])

x0 = 0.09;                  % left position
dy = .03;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .07;
w2 = (w-1*dw-2*x0)/2;               % width per figure
ht = (y0-6*dy);          % height per figure

y3 = y0-dy-ht+0;
set(h(1),'Position',[x0 y3 w2 ht])
set(h(2),'Position',[x0+w2+dw y3 w2 ht])



    for k = 1:2
        %if k < 5
       %     set(h(k),'XTick',[])
       % end
        set( h(k)                       , ...
          'FontName'   , 'Helvetica' ,...
          'Box'         , 'off'     , ...
          'TickDir'     , 'out'     , ...
          'TickLength'  , [.02 .02] , ...
          'XColor'      , [.3 .3 .3], ...
          'YColor'      , [.3 .3 .3], ...
          'LineWidth'   , 1         );
    end


    % export to eps

    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 -r600 -painters /home/ig/neef/Documents/AAM_Paper/Fig12.eps

    % less shitty gridlines

    fixPSlinestyle('/home/ig/neef/Documents/AAM_Paper/Fig12.eps')




