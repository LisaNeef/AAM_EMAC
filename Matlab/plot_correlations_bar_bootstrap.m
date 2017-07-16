 function plot_correlations_bar_v2(tscale,alpha,plot_type,export)
% plot_correlations_bar.m
%
% Make a bar graph of the correlations between different data sets and 
% the observed ERPs, at different timescales.
% _bootstrap: in this version we plot to correlation confidence estimates
% estimated via the bootstrap method.
%
% Mods: 
%   12 Sept 2011: make more customizable with external call.
%
%  INPUTS:
%	tscale : which timescale.
%         tscale = 1: 30-90 day (subseasonal) variations.
%         tscale = 2: 24-30 month variations (QBO).
%         tscale = 3: 2-7 year variation.
%         tscale = 4: 7-20 year variations.
%-------------------------------------------------------------------------

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
SIG = zeros(nf,3);
CI  = zeros(nf,2,3);
LAG = zeros(nf,3);
TAU = zeros(nf,3);

%---loop over timescales and compute the correlations.
for comp = 1:3
  [corr_out,corr_ci,corr_std,Nsamples, tau,runs_out,C,L,R] = correlations_aam_filtered_bootstrap(comp,tscale,alpha) ;
  RHO(:,comp) = corr_out;
  CI(:,:,comp) = corr_ci;
  SIG(:,comp) = corr_std;
end

% also retrieve the names of the runs to be plotted here
[runs,names,ef,no_ib] = aam_paper_runs;

%---compute the upper and lower error bounds based on the confidence
%interval given.
error_lo = RHO-squeeze(CI(:,1,:));
error_hi = RHO-squeeze(CI(:,2,:));

%---some plot settings.

col = aam_paper_colors;

col(5,:) = 0.7*ones(1,3);

if tscale == 1, tscalename = 'Subseasonal'; end
if tscale == 2, tscalename = 'Quasi-Biennial'; end
if tscale == 3, tscalename = 'Interannual'; end
if tscale == 4, tscalename = 'Long-Term'; end

TT = [tscalename, ' Variations'];
Xnames = {'X1', 'X2', 'X3'};


FS = 18;

% Print out the number of samples, decorrelation times, and effective DOF
%disp('Run Name   # Samples  Decorr Time  correlation       sig. limit')
%for ii = 1:nf
%  disp([char(names(R(ii))),'  ',num2str(N(ii)),'  ',num2str(tau(ii)),'  ',num2str(pmax(ii)),'   ',char(corr_sig(ii))])
%end 
%disp(['Mean decorreation time:   ',num2str(mean(tau))])



%---plots!

barweb(RHO', error_lo',error_hi', 1, ['X1';'X2';'X3'], [], [], [], col(R,:), 'y', names(R), 2, 'plot','EastOutSide')
  ylabel('correlation','FontSize',FS)
  text(0.6,0.9,TT,'FontSize',FS)
  axis([0.5 3.5 -0.5 1])
  grid off
disp(TT)


%---export these plots.
if export
    plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

    fig_name_1 = [plot_dir,'correlations_bar_X',tscalename,'_bootstrap_CI.png'];

       LW = 2;
    ph = 10;        % paper height
    pw = 20;        % paper width
    fs = 30;        % fontsize

    exportfig(1,fig_name_1,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
end


