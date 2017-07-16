function plot_aam_mm_regions_ERAdiffs(comp)
% plot_aam_mm_regions_ERAdiffs.m--------------
%
% 8 Oct 2010
% Mods
%  14 Nov 2010: add option to look at diffs in era winds and ps fields.
%  19 Nov 2010: make the plot exporting automatic
% INPUTS:
%    comp: choose which component to plot 
%	1,2,3: components of AAM (X1, X2, or DLOD)
%	4: u300 
%	5: ps
%
%-----------------------------------------------------------------

%---paths...
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/bipolar_colormap');

% load the data 
ccmval20 = 1;	% flag to only use 20 years of ccmval run.
if comp > 3, comp2 = 1; else comp2 = comp; end
[lat,lon,X_emac,u_emac,ps_emac] = aam_mm_slice(comp2,'CCMval',ccmval20);
[lat,lon,X_era,u_era,ps_era] = aam_mm_slice(comp2,'ERAinterim',0);

  [nm,nlat,nlon] = size(u_emac);
  D = zeros(3,nm,nlat,nlon);

if comp < 4
   x_era = sum(X_era,1);
   x_emac = sum(X_emac,1);
   D(1,:,:,:) = x_era;
   D(2,:,:,:) = x_emac;
   D(3,:,:,:) = x_era-x_emac;
end
if comp == 4
  D(1,:,:,:) = u_era;
  D(2,:,:,:) = u_emac;
  D(3,:,:,:) = u_era-u_emac;
end
if comp == 5
  D(1,:,:,:) = ps_era;
  D(2,:,:,:) = ps_emac;
  D(3,:,:,:) = ps_era-ps_emac;
end

%--plot specifications:
if comp == 1
  T4jan = {'\chi_1: JAN, ERA','\chi_1: JAN, EMAC', '\chi_1: JAN, ERA-EMAC'};
  T5jan = {'\chi_1: JAN, ERA','\chi_1: JAN, EMAC', '\chi_1: JAN, ERA-EMAC'};
  Tjan = [T4jan;T4jan;T4jan;T4jan;T5jan];
  T4jul = {'\chi_1: JUL, ERA','\chi_1: JUL, EMAC', '\chi_1: JUL, ERA-EMAC'};
  T5jul = {'\chi_1: JUL, ERA','\chi_1: JUL, EMAC', '\chi_1: JUL, ERA-EMAC'};
  Tjul = [T4jul;T4jul;T4jul;T4jul;T5jul];
end
if comp == 4
  T4jan = {'u(300hPa): JAN, ERA','u(300hPa): JAN, EMAC', 'u(300hPa): JAN, ERA-EMAC'};
  T5jan = {'u(300hPa): JAN, ERA','u(300hPa): JAN, EMAC', 'u(300hPa): JAN, ERA-EMAC'};
  Tjan = [T4jan;T4jan;T4jan;T4jan;T5jan];
  T4jul = {'u(300hPa): JUL, ERA','u(300hPa): JUL, EMAC', 'u(300hPa): JUL, ERA-EMAC'};
  T5jul = {'u(300hPa): JUL, ERA','u(300hPa): JUL, EMAC', 'u(300hPa): JUL, ERA-EMAC'};
  Tjul = [T4jul;T4jul;T4jul;T4jul;T5jul];
end
if comp == 5
  T4jan = {'p_s: JAN, ERA','p_s: JAN, EMAC', 'p_s: JAN, ERA-EMAC'};
  T5jan = {'p_s: JAN, ERA','p_s: JAN, EMAC', 'p_s: JAN, ERA-EMAC'};
  Tjan = [T4jan;T4jan;T4jan;T4jan;T5jan];
  T4jul = {'p_s: JUL, ERA','p_s: JUL, EMAC', 'p_s: JUL, ERA-EMAC'};
  T5jul = {'p_s: JUL, ERA','p_s: JUL, EMAC', 'p_s: JUL, ERA-EMAC'};
  Tjul = [T4jul;T4jul;T4jul;T4jul;T5jul];
end

nt = 3;

%---set the axes for the different terms here
cax = zeros(6,2);
cax(4,:) = 50*[-1,1];	% axis for u300 plot


iplot=1;
for im = [1,7];
  for iterm = 1:nt
    figure(iplot),clf
    axesm robinson; framem;gridm
    c = load('coast');
    geoshow(c.lat,c.long,'Color',0*ones(3,1))
    contourfm(lat,lon,squeeze(D(iterm,im,:,:)),30); 
    caxis(cax(comp,:))
    if iterm == nt, caxis(0.4*cax(comp,:)), end	% (zoom for difference plot)
    if im == 1, title(Tjan(comp,iterm)), end
    if im == 7, title(Tjul(comp,iterm)), end
    %colormap(bluewhitered(256)), colorbar('SouthOutside')
    colormap(jet(256)), colorbar('West')
    iplot = iplot + 1;
  end
end


%----export plots!

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';
if comp == 4, varname = 'u300'; end
month = {'jan','jan','jan','jul','jul','jul'};
dataset = {'ERA','EMAC','ERA-EMAC','ERA','EMAC','ERA-EMAC'};
LW = 4;
ph = 10;        % paper height
pw = 10;        % paper width
fs = 30;        % fontsize

for ifig = 1:6
  fig_name = [plot_dir,char(varname),'_',char(month(ifig)),'_',char(dataset(ifig)),'.png']
  exportfig(ifig,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
end


