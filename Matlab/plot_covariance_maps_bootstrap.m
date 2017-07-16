function plot_covariance_maps_bootstrap(runid,tscale,comp,term,abs_val,fil_order,plot_label,export,plot_extras)
%plot_correlation_maps(runid,tscale,comp)
%
% load pre-computed maps of correlation and covariance and plot on a map.
% in this version, plot covariances (rather than correlations)
%
% MODS:
%
% INPUTS
%  runid: which run to do the map for
%  tscale: which timescale to plot
%  comp: which vector component to plot
%  abs_val: set to 1 to have absolute value
%  fil_order: which filter order file to load

%----------------------------------------------------

% temp stuff
%clear all;
%runid = 'CCMval';
%fil_order = 2;
%abs_val =  1;
%comp = 2;
%tscale = 2;
%term = 2;
%plot_label = '(d)';
% export = 0;

% set default parameters
switch nargin
  case 3
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
ff = dir([datadir,'cov_EF_X',num2str(comp),'*',tscalename,'_filorder',num2str(fil_order),'*_bootstrap500.mat']);
filename = [datadir,ff.name];
if ~isempty(ff)
  disp(['loading file ',filename])
  load(filename)
else
  disp(['Cant find the input file ',filename])
  return
end

%---show abosolute values instead? 

if abs_val
  SIG2 = abs(SIG);
else
  SIG2 = SIG;
end

% mask out everything that comes out as not significant according to
% bootstrap
% the way we test for this is: if there is a sign-switch between the
% 98-percentile limes of each estimate (S_LO and S_HI), then 0 is included
% in the 98-percentile correlation estimates.
% And how to we see a sign change? A: if the product of the two numbers
% makes a negative.  (29 Aug 2011)

insignificant = (S_LO.*S_HI) <= 0;
SIG2(insignificant) = 0;

%--- figure settings
LW = 1;
ph = 10;        % paper height
pw = 10;        % paper width
FS = 20;        % fontsize


%make sig_max a function of the timescale
switch tscale
    case 1
        switch comp
            case {1,2}
                sig_max = 50;
            case 3
                sig_max = 4E-3;
        end
    case 2
        switch comp
            case {1,2}
                sig_max = 600;
            case 3
                sig_max = 0.1;
        end
    case 3
        switch comp
            case {1,2}
                sig_max = 50;
            case 3
                sig_max = 4E-3;
        end
end

if abs_val
  cor_limits = sig_max*[0,1];
else
  cor_limits = sig_max*[-1,1];
end

% set colors for the plots.

gray = 0.7*ones(1,3);
nlevels_corr = 10;
if ~abs_val
    nlevels_corr = 11;
end

% pointer to nonzero longitudes -- needed in plot to avoid a vertical
% stipe.
gg = find(abs(lon) > 0);



%-- make plots


  ii = term;
 % figure(1)
 % subplot(spv)
  ax = axesm('MapProjection','eqdcylin','grid','on',...
      'MeridianLabel','on','ParallelLabel','on',...
      'PLabelLocation',30,'MLabelLocation',90);
   
% NOTE: sometgimes the plot disappears when we add a grid.
% this prob seems to only be present in matlab2010b (not a)
      
  contourfm(lat,lon(gg),squeeze(SIG2(ii,:,gg)),nlevels_corr,'LineStyle','none');
  caxis(cor_limits);
 
  if abs_val
      col_map_corr = jet(nlevels_corr);
      col_map_corr(1,:) = ones(1,3);
  else
%      col_map_corr  = makeColorMap(col_neg,col_zero,col_pos,nlevels_corr);
      col_map_corr = jet(nlevels_corr);
      col_map_corr(round(nlevels_corr/2),:) = ones(1,3);
  end  

  colormap(col_map_corr)
  setm(ax,'Fontsize',5);
  setm(ax,'Fontcolor','Black');
  c = load('coast');
  plotm(c.lat,c.long,'Color','black','LineWidth',1);

  % for mass terms, mask out the sea.  This is done with a weird trick...
  if ii == 2
    geoshow(flipud(c.lat), flipud(c.long),'DisplayType', 'polygon', 'FaceColor', gray)
  end

  % clean up a little bit
  box off
  axis off
  grid on
%  if plot_extras
     colorbar('location','EastOutside'); 
%  end
  
 title(plot_label)



%---export plots
if export
    
plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

      if ~abs_val
          switch ii
              case 1
                suff = [num2str(comp),'w_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png'];
              case 2
                suff = [num2str(comp),'m_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png'];
              case 3
                suff = [num2str(comp),'t_',tscalename,'_',runid,'_filorder',num2str(fil_order),'.png'];
          end
      else
          switch ii
              case 1
                suff = [num2str(comp),'w_',tscalename,'_',runid,'_filorder',num2str(fil_order),'_absval.png'];
              case 2
                suff = [num2str(comp),'m_',tscalename,'_',runid,'_filorder',num2str(fil_order),'_absval.png'];
              case 3
                suff = [num2str(comp),'t_',tscalename,'_',runid,'_filorder',num2str(fil_order),'_absval.png'];
          end
      end
    fig_name = [plot_dir,'corrmap_chi',suff];
    exportfig(gcf,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',FS,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');

end

