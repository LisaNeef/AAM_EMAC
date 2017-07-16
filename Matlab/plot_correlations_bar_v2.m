% function plot_correlations_bar_v2(tscale,plot_type)
% plot_correlations_bar.m
%
% Make a bar graph of the correlations between different data sets and 
% the observed ERPs, at different timescales.
% started 22 Nov 2010
% Version 2: instead of lumping plots by component, do it by timescale.
% Mods:
%  8 Dec 2010: switch back to plotting correlation instead of variance explained.
% 17 Dec 2010: streamline the code and update colors.
% 11 Mar 2011: eliminate the choice of central correlation, for simplicity (notes vol. 4, p. 118)
%		Also add plot_type as input
%  INPUTS:
%	tscale : which timescale.
%         tscale = 1: 30-90 day (subseasonal) variations.
%         tscale = 2: 24-30 month variations (QBO).
%         tscale = 3: 2-7 year variation.
%         tscale = 4: 7-20 year variations.
%	plot_type: which plot to make?
%	  1: compare all EMAC runs, ERAs, OAM and HAM 
%	  2: as in 1, but also add the best-fitting CAM run (Need to work on this)
%-------------------------------------------------------------------------

% temp inputs
tscale = 3;
plot_type = 1;

switch plot_type
  case 1
    nf = 8;
    exp_type = 1;
  case 2
    nf = 9;
    exp_type = 3;
end

%--initialize arrays holding all the data.
RHO = zeros(nf,3);
LAG = zeros(nf,3);

%---loop over timescales and compute the correlations.
for comp = 1:3
  [pmax,lag,N,tau,runs,C,L,R,corr_sig] = correlations_aam_filtered_v2(comp,tscale,exp_type);
  LAG(:,comp) = lag; 
  RHO(:,comp) = pmax;
  % mask out correlations that are not significant at the 95% level.
  RHO(find(corr_sig == -1),comp) = NaN;
  RHO(find(pmax < corr_sig'),comp) = NaN;
  LAG(find(corr_sig == -1),comp) = NaN;
  LAG(find(pmax < corr_sig'),comp) = NaN;

end

% also retrieve the names of the runs to be plotted here
[runs,names,ef,no_ib] = aam_paper_runs;

%---some plot settings.

col = aam_paper_colors;

if tscale == 1, tscalename = 'Subseasonal'; end
if tscale == 2, tscalename = 'Quasi-Biennial'; end
if tscale == 3, tscalename = 'Interannual'; end
if tscale == 4, tscalename = 'Long-Term'; end

TT = [tscalename, ' Variations'];
Xnames = {'X1', 'X2', 'X3'};


% Print out the number of samples, decorrelation times, and effective DOF
disp('Run Name   # Samples  Decorr Time  correlation       sig. limit')
for ii = 1:nf
  disp([char(names(R(ii))),'  ',num2str(N(ii)),'  ',num2str(tau(ii)),'  ',num2str(pmax(ii)),'   ',char(corr_sig(ii))])
end 
disp(['Mean decorreation time:   ',num2str(mean(tau))])

% compute mean Neff and, from there, significance limit.
% here don't count ERA-interim, OAM, or HAM
for ii = 1:nf
  name = char(names(R(ii)));
  if length(name) == 6
    if name(1:3) == 'HAM'; ham = ii; end
    if name(1:3) == 'OAM'; oam = ii; end
    if name(5:6) == 'In'
      erainterim = ii;
    end
  end
end
Rmod = R*0+1;
Rmod(erainterim) = 0;
Rmod(ham) = 0;
Rmod(oam) = 0;



Neff_mean = mean(Rmod'.*Neff);
disp(['Mean N_eff:   ',num2str(Neff_mean)])

sig_cutoff = erfinv(0.95)*sqrt(2/Neff_mean);

%--- make a bar graph.

figure(1),clf
  bar(RHO')
  ylabel('correlation')
  text(0.6,0.9,TT)
  set(gca,'XTickLabel',Xnames)
  legend(names(R))
  legend('Location','EastOutside')
  colormap(col(R,:))
  axis([0.5 3.5 0 1])
 
   figure(2),clf
    bar((LAG)')
    ylabel('Optimal Lag (Days)')
    title(TT)
    set(gca,'XTickLabel',Xnames)
    legend(names(R))
    legend('Location','EastOutside')
    colormap(col(R,:))

disp(TT(comp))
%---export these plots.

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

fig_name_1 = [plot_dir,'correlations_bar_X',tscalename,'_','.png'];
fig_name_2 = [plot_dir,'correlations_bar_X',tscalename,'_lag.png']; 

LW = 4;
ph = 10;        % paper height
pw = 20;        % paper width
fs = 30;        % fontsize

exportfig(1,fig_name_1,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

exportfig(2,fig_name_2,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

