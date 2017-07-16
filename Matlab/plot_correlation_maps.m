%function plot_correlation_maps(runid,tscale,comp,corr_type,abs_val,fil_order)
%plot_correlation_maps(runid,tscale,comp)
%
% load pre-computed maps of correlation and covariance and plot on a map.
% 13 Dec 2010, Notes vol. 4, p. 51.
%
% MODS:
%  change the colormap to something easier to read.
%	for correlations: diverging scheme with white center
%	for covariances (which are positive-definite): single-hue progression 
%  14 Mar 2011: update paths to CCMVal data
%  14 Mar 2011: make plotting correlations vs covariances optional
%  15 Mar 2011: for mass term cov/corr, mask out the oceans.
%  15 Mar 2011: add option for showing the absolute values.
%  16 Jun 2011: change the filename call to make sure it chooses
%  correlation maps with AAM.
%
% INPUTS
%  runid: which run to do the map for
%  tscale: which timescale to plot
%  comp: which vector component to plot
%
% OPTIONAL INPUTS:
%  corr_type: 1 for correlation, 2 for covariance
%  abs_val: set to 1 to have absolute value
%  fil_order: which filter order file to load

%----------------------------------------------------

% temp input
runid = 'CCMval';
tscale = 1;
comp = 3;
corr_type = 1;
fil_order = 2;
abs_val = 0;

% set default parameters
switch nargin
  case 3
    corr_type = 1;
    abs_val = 1;
    fil_order = 2;
  case 2
    abs_val = 1;
    fil_order = 2;
  case 1
    fil_order = 2;
end  
  

% load preprocessed data.
if tscale == 1, tscalename = 'Subseasonal'; end
if tscale == 2, tscalename = 'Annual'; end
if tscale == 3, tscalename = 'Interannual'; end
if tscale == 4, tscalename = 'Long-Term'; end


datadir = ['/dsk/nathan/lisa/EMAC_ex/',runid,'/mat/'];
ff = dir([datadir,'corr_aam_X',num2str(comp),'*',tscalename,'_filorder',num2str(fil_order),'*.mat']);
filename = [datadir,ff.name]
if length(ff) > 0
  disp(['loading file ',filename])
  load(filename)
else
  disp(['Cant find the input file ',filename])
  return
end

%---show abosolute values instead? 

if abs_val
  RHO2 = abs(RHO);
  SIG2 = abs(SIG);
else
  RHO2 = RHO;
  SIG2 = SIG;
end

%--- figure settings
LW = 1;
ph = 10;        % paper height
pw = 10;        % paper width
FS = 20;        % fontsize

YL = {'mas^2';'mas^2';'ms^2'};

% set axis limits -- based on maxima?
rho_max = zeros(1,3);
sig_max = zeros(1,3);
for ii = 1:3
  rho_max(ii) = max(max(abs(RHO(ii,:,:))));
  sig_max(ii) = max(max(abs(SIG(ii,:,:))));
end

% for covariances it's good to adjust the limits by timescale.
%if tscale == 1,  sig_max = [0.5e6, 0.5e6, 20]; end
sig_max = [0.5e6, 0.5e6, 30];
if tscale == 2,  sig_max = [5e6, 5e6, 30]; end

% set axis limits by tscale?
% default values first
sig_max = [5e5, 5e5, 30];
rho_max = [1,1,1];
switch tscale
    case 1
        rho_max = [0.2,0.2,0.2];
    case 2
        sig_max = [5e6,5e6,30];
        rho_max = [0.5,0.5,0.5];
    case 3
        rho_max = [0.5,0.5,0.5];
end

        

if abs_val
  cor_limits = rho_max'*[0,1];
  cov_limits = sig_max'*[0,1];
else
  cor_limits = rho_max'*[-1,1];
  cov_limits = sig_max'*[-1,1];
end

% set colors for the plots.

gray = 0.7*ones(1,3);

col_neg = [0.2622    0.6028    0.7112]; % dark teal
%col_neg = [0.030   0.353   1.000];	% blue
%col_neg = [0, 136, 55]/356;		% green

col_pos = [0.9961    0.0782    0.4427]; % fuschia
%col_pos = [ 0.9134    0.6324    0.0975]; % mustard
%col_pos = [123, 50, 148]/356;		% purple
%col_pos = [0.9649    0.1576    0.9706]; % magenta

col_zero = ones(1,3);
col_high = [0.3500    0.1966    0.2511];

nlevels_corr = 15;
nlevels_cov = 7;

if abs_val
  col_map_corr  = makeColorMap(col_zero,col_neg,col_pos,10); 
  col_map_cov  = col_map_corr;
else
  col_map_corr  = makeColorMap(col_neg,col_zero,col_pos,nlevels_corr); 
  col_map_cov  = makeColorMap(col_neg,col_zero,col_pos,nlevels_corr); 
end

%-- make plots


for ii = 1:3		% loop over mass, wind, and total
  figure(0+ii),clf
  ax = axesm('MapProjection','eqdcylin','grid','on',...
      'MeridianLabel','on','ParallelLabel','on',...
      'PLabelLocation',30,'MLabelLocation',90)
   
% NOTE: sometgimes the plot disappears when we add a grid.
% this prob seems to only be present in matlab2010b (not a)
      
  switch corr_type
    % correlation plots
    case 1
      contourfm(lat,lon,squeeze(RHO2(ii,:,:)),nlevels_corr,'LineStyle','none');
      caxis(cor_limits(comp,:));
      title_pref = 'Correlation with Global \chi_';

    % covariance plots
    case 2
      contourfm(lat,lon,squeeze(SIG2(ii,:,:)),nlevels_cov,'LineStyle','none');
      caxis(cov_limits(comp,:));
      title_pref = 'Covariance with Global \chi_';

  end
  
  colormap(col_map_corr);
  setm(ax,'Fontsize',16);
  setm(ax,'Fontcolor','Black');
  c = load('coast');
  plotm(c.lat,c.long,'Color','black','LineWidth',1);

  % for mass terms, mask out the sea.  This is done with a weird trick...
  if ii == 2
    % geoshow(flipud(c.lat),flipud(c.long),'DisplayType','Line')
    % geoshow('seaareas.shp', 'FaceColor', 'black');
     geoshow(flipud(c.lat), flipud(c.long),'DisplayType', 'polygon', 'FaceColor', gray)

  end

  % clean up a little bit
  box off
  axis off
  grid on
  h = colorbar; % Make colorbar and save handle.
  switch ii
    case 1
      title([title_pref,num2str(comp),', Wind Term,  ',num2str(tscalename)]);
    case 2
      title([title_pref,num2str(comp),', Mass Term,  ',num2str(tscalename)]);
    case 3
      title([title_pref,num2str(comp),'  ',num2str(tscalename)]);
  end

end


%---export plots

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

  for ii = 1:3
    if ii == 1, suff = [num2str(comp),'w_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png']; end
    if ii == 2, suff = [num2str(comp),'m_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png']; end
    if ii == 3, suff = [num2str(comp),'t_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png']; end
    if corr_type == 1
      fig_name = [plot_dir,'corrmap_chi',suff];
    else
      fig_name = [plot_dir,'covmap_chi',suff];
    end
    exportfig(ii,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
  end



