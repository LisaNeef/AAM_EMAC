%% make cov_corr_maps_together.m
% Different drafts of the correlation (or covariance) maps to present in
% in the AAM paper, for each timescale. --> which format is best?
% 3 Nov 2011

%% main settings to change

clear all;
tscale = 2;
to_show = 2;       % 1 = correlation, 2 = covariance

%% code that doesn't change

runid = 'CCMval';
fil_order = 2;
abs_val =  1;

figure(1),clf

h = zeros(1,6);

PLm = {'(a) \chi_1 Mass Term', '(c) \chi_2 Mass Term','(e) \chi_3 Mass Term'};
PLw = {'(b) \chi_1 Wind Term', '(d) \chi_2 Wind Term','(f) \chi_3 Wind Term'};

jj = 1;

extras = zeros(3,2);        % plot extras: set to 1 to get colorbar, etc
extras(3,:) = 1;

% make the plots

for comp = 1:3
    for term = [2,1]
        switch term
            case 1
                plot_label = PLw(comp);
            case 2
                plot_label = PLm(comp);
        end
        h(jj) = subplot(3,2,jj);
        switch to_show
            case 1
                plot_correlation_maps_bootstrap(runid,tscale,comp,term,abs_val,fil_order,plot_label,0,extras(comp,term))
            case 2
                plot_covariance_maps_bootstrap(runid,tscale,comp,term,abs_val,fil_order,plot_label,0,extras(comp,term))
        end
        jj = jj+1; 
      end 
end



%% make them look good.



set(gcf,'Position',[911 551 845 719])

x0 = 0.09;                  % left position
dy = .03;
y0 = 0.99;                  % top of plot
w = 1; 
dw = .03;
w2 = (w-1*dw-2*x0)/2;               % width per figure
ht = (y0-6*dy)/3;          % height per figure

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
          'Box'         , 'off'     , ...
          'TickDir'     , 'out'     , ...
          'TickLength'  , [.02 .02] , ...
          'XColor'      , [.3 .3 .3], ...
          'YColor'      , [.3 .3 .3], ...
          'LineWidth'   , 1         );
    end


%% export to eps

    set(gcf, 'PaperPositionMode', 'auto');
    
switch to_show
    case 1
        switch tscale
            case 1
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/subseasonal_correlation.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/subseasonal_correlation.eps')  
            case 2
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/annual_correlation.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/annual_correlation.eps')  
            case 3
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/IAV_correlation.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/IAV_correlation.eps')  
        end
    case 2
        switch tscale
            case 1
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/subseasonal_covariance.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/subseasonal_covariance.eps')  
            case 2
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/annual_covariance.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/annual_covariance.eps')  
            case 3
                print -depsc2 -r600 -painters /home/ig/neef/Documents/Plots/aam_run_comparison/IAV_covariance.eps
                fixPSlinestyle('/home/ig/neef/Documents/Plots/aam_run_comparison/IAV_covariance.eps')  
        end
end



