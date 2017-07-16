%% make_TF_polar_stereogr_plots.m
% Make grouped figures showing the AM transfer (weighting) functions in a
% polar stereographic projection.  Plots are grouped according to variable,
% and for each variable, each AM component is shown.
%
% Lisa Neef, 13 April 2012
%



%% User inputs

variable = 'U';
hostname = 'ig48';
clc;

%% cycle AM components and produce plots

if strcmp(variable,'V')
    C = {'X1','X2'};
    nc = 2;
else
    C = {'X1','X2','X3'};
    nc = 3;
end

h = zeros(1,nc);

switch hostname
    case 'ig48'
        figH = figure(1);clf
    case 'blizzard'
        figH = figure('visible','off') ; 
end

for ii = 1:nc
    h(ii) = subplot(1,nc,ii);
    plot_TF_polar_stereogr(char(C(ii)),variable)
    freezeColors
    box off
    axis off
end


%% adjust the axes, make everything look right.
set(gcf,'Position',[60 576 1082 341])

x0 = 0.01;                  % left position
dy = 0;
y0 = 0.97;                  % top of plot
w = 1; 
dw = .07;
ht = (y0-dy);          % height per figure
y = y0-dy-ht+0;

    pw = 12;        % paper width
    w2 = (w-3*dw-x0)/3;               % width per figure
    set(h(1),'Position',[x0 y w2 ht])
    set(h(2),'Position',[x0+w2+dw y w2 ht])
if length(h) == 3
    set(h(3),'Position',[x0+2*w2+2*dw y w2 ht])
end


    

   for k = 1:length(h)
    
       set( h(k)                       , ...
       'FontName'   , 'Helvetica' ,...
       'FontSize'   , 16         ,...
       'Box'         , 'off'     , ...
       'YGrid'       , 'on'      , ...
       'XGrid'       , 'on'      , ...
       'TickDir'     , 'out'     , ...
       'TickLength'  , [.02 .02] , ...
       'XColor'      , [.3 .3 .3], ...
       'YColor'      , [.3 .3 .3], ...
       'LineWidth'   , 1         );
   end



%% some plot settings
LW = 1;
ph = 4;        % paper height
fs = 16;        % fontsize

switch hostname
    case 'ig48'
        plot_dir = '/home/ig/neef/Documents/Plots/AAM/';
    case 'blizzard'
        plot_dir = '/work/bb0519/b325004/DART/ex/Comparisons/';
end

fig_name = ['TF_',variable,'_polarstereogr.png'];


%% export

    exportfig(figH,[plot_dir,fig_name],'width',pw,'height',ph,...
              'fontmode','fixed', 'fontsize',fs,'color','cmyk',...
              'LineMode','fixed','LineWidth',LW,'format','png');
    if strcmp(hostname,'blizzard')
         close(figH) ;
    end




