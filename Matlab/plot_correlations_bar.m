% function plot_correlations_bar(comp,bw)
% plot_correlations_bar.m
%
% Make a bar graph of the correlations between different data sets and 
% the observed ERPs, at different timescales.
% started 22 Nov 2010
%
% Mods:
%    24 Nov 2010 - instead of correlation, let's plot its square -- the
%          percent of explained variance.
%    27 Jun 2011 - back to v1, and back to plotting correlations instead of
%    explained variance.
%
%  INPUTS:
%	comp: vector component of AAM
%	bw: set to 1 for black and white plots
%-------------------------------------------------------------------------


%--temp inputs
clear all;
comp =3;
bw = 0;

%--initialize arrays holding all the data.
nf = 8;
RHO = zeros(nf,3);
LAG = zeros(nf,3);


%---loop over timescales and compute the correlations.
tscales = [1,3,2];
for itscale = 1:3
  [pmax,lag,Neff, tau,runs_out,C,L,R,corr_sig] = correlations_aam_filtered_v2(comp,tscales(itscale),1);
  LAG(:,itscale) = lag; 
  RHO(:,itscale) = pmax;
  % mask out correlations that are not significant at the 95% level.
  RHO(find(corr_sig == -1),itscale) = NaN;
  RHO(find(pmax < corr_sig'),itscale) = NaN;
  LAG(find(corr_sig == -1),itscale) = NaN;
  LAG(find(pmax < corr_sig'),itscale) = NaN;
end


% also retrieve the names of the runs to be plotted here
[runs,names,ef,no_ib] = aam_paper_runs;



%---some plot settings.

col = aam_paper_colors;

switch comp
    case 1
        TT = '\chi_1';
    case 2
        TT = '\chi_2';
    case 3
        TT = '\chi_3';
end

Tnames = {'Subseasonal (30-90 d)', '2-7 y', 'QBO'};


%--- make a bar graph.


figure(1),clf
  bar(RHO')
  ylabel('correlation')
  title(TT)
  set(gca,'XTickLabel',Tnames)
  legend(names(R))
  legend('Location','EastOutside')
  colormap(col(R,:))
  
figure(2),clf
    bar(LAG')
    ylabel('Optimal Lag (Days)')
    title(TT)
    set(gca,'XTickLabel',Tnames)
    legend(names(R))
    legend('Location','EastOutside')
    colormap(col(R,:))

disp(TT(comp))


%---export these plots.

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

fig_name_1 = [plot_dir,'correlations_bar_X',num2str(comp),'_','.png'];
fig_name_2 = [plot_dir,'correlations_bar_X',num2str(comp),'_lag.png']; 

LW = 4;
ph = 10;        % paper height
pw = 20;        % paper width
fs = 30;        % fontsize

exportfig(1,fig_name_1,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

exportfig(2,fig_name_2,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

