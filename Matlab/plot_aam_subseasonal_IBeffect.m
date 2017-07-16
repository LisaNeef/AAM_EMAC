function plot_aam_subseasonal_IBeffect(comp)
%
% Comparing the effect of having the IP approximation in 
% subseasonal AAM Variations
% 6 Dec 2010
% Notes vol. 4, p. 41
%
% INPUTS:
%  comp: which vector component to look at (1,2, or 3)
%-----------------------------------------------------

% define the run names and other settings

nf = 3;
runid = {'CCMval';'CCMval';'ERA40'};
names = {'T42L90';'T42L90-noIB';'ERA-40'};
emac_flag = [1,1,0];
no_ib = [0,1,0];

T = [30,90]; Ty = (1/365)*T; 
date_start = [1980,1,1];
date_stop = [1985,1,1];

fil_order = 2;		% greater than 2 seems to not work for ERA data
filtype = 1;			% butterworth filter
%filtype = 2;			% chebyeshev-1 filter
%filtype = 3;			% chebyeshev-2 filter


% initialize empty array

mjd0 = date2mjd(date_start(1),date_start(2),date_start(3));
mjdf = date2mjd(date_stop(1),date_stop(2),date_stop(3));
nd = mjdf-mjd0+1;
XX = zeros(nf,3,nd);
XXF = zeros(nf,3,nd);

% retrieve the filtered runs

for irun = 1:nf

  [X,XF,MJD] = retrieve_AAM_filtered(runid(irun),T,filtype,fil_order,date_start,date_stop, emac_flag(irun),no_ib(irun));
  XX(irun,:,:) = X(:,comp,:);
  XXF(irun,:,:) = XF(:,comp,:);

end

% craft a time series.
t = MJD*0+NaN;
[y, m, d] = mjd2date(MJD);
for ii=1:length(t)
  if isfinite(MJD(ii)) ==1, t(ii)=datenum([y(ii) m(ii) d(ii)]); end
end

% set plot settings

YL = {'\chi_1 (mas)', '\chi_2 (mas)', '\Delta LOD (ms)'};
TT1 = {'\chi_1 Wind Term 30-90 Day Variation','\chi_1 Mass Term 30-90 Day Variation','\chi_1 30-90 Day Variation'};
TT2 = {'\chi_2 Wind Term 30-90 Day Variation','\chi_2 Mass Term 30-90 Day Variation','\chi_2 30-90 Day Variation'};
TT3 = {'\Delta LOD Wind Term 30-90 Day Variation','\Delta LOD Mass Term 30-90 Day Variation','\Delta LOD 30-90 Day Variation'};
TT = [TT1;TT2;TT3];

aam_paper_colors
transparency = 0.3;

rend = 'opengl';
LW = 4;
ph = 8;        % paper height
pw = 25;        % paper width
fs = 26;        % fontsize

ax(1,:) = 45*[-1,1];
ax(2,:) = 50*[-1,1];
ax(3,:) = 0.5*[-1,1];

col2 = col;
col2(1,:) = col(4,:);
col2(2,:) = col(11,:);
col2(3,:) = col_era;

% make the plots

for iterm = 1:3
  lh = zeros(2,1);
  figure(iterm), clf
  for irun =1:nf
    x = squeeze(XXF(irun,iterm,:));
    gg = find(isfinite(x) == 1);
    if length(gg) > 0
      lh(irun) = plot(t(gg),x(gg),'LineWidth',LW,'Color',col2(irun,:));
      hold on
    end
  end
  datetick('x','mm/yy')
  ylabel(YL(comp))
  xlabel('time')
  title(TT(comp,iterm))
  set(gca,'LineWidth',3.0)
  lg = find(isfinite(lh)); 
  legend(lh(lg),char(names(lg)),'Orientation','Horizontal','Location','SouthOutside')
  legend('boxoff')
  axis([t(1),max(t), ax(comp,1), ax(comp,2)]);

end


% export the plots

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
fig_name_w = [plot_dir,'comp_subseasonal_IBeffect_X',num2str(comp),'w.png'];
fig_name_m = [plot_dir,'comp_subseasonal_IBeffect_X',num2str(comp),'m.png'];
fig_name_t = [plot_dir,'comp_subseasonal_IBeffect_X',num2str(comp),'t.png'];

rend = 'opengl';
LW = 4;
ph = 10;        % paper height
pw = 25;        % paper width
fs = 30;        % fontsize

exportfig(1,fig_name_w,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');
exportfig(2,fig_name_m,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');
exportfig(3,fig_name_t,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','renderer',rend,'LineMode','fixed','LineWidth',LW,'format','png');



