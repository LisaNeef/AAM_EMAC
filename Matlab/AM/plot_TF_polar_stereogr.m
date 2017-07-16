%% plot_TF_polar_stereogr.m
function plot_TF_polar_stereogr(comp,variable)
%
% plot_TF_polar_stereogr.m
% make plots of the AAM excitation functions in polar-stereographic form.
% 
% Lisa Neef, 13 April 2012
%
% INPUTS:
%   comp: angular momentum component (X1,X2,X3)
%   variable: 'U', 'V', or 'PS'
%
% MODS:
%  1 Feb 2013: flip the sign on the TFs for X1 and X2, in order to include that 
%      negative prefactor in these integrals.


%% temp inputs
%comp = 'X1';
%variable = 'V';

%% retrieve the desired excitation function

lat = -90:.5:90;
lon = -180:.5:180;

W = eam_weights(lat,lon,comp,variable);

[X,Y] = meshgrid(lon,lat);

%% multiply the transfer functions for X1 and X2 by -1, to account for the negative prefactor.
switch comp
case {'X1','X2'}
	W2 = -W;
case 'X3'
	W2 = W;
end


%% make the polar stereographic plot.

axesm ('MapProjection','stereo', 'Frame', 'on', 'Grid', 'on','MapLatLimit',[0 90]);
contourfm(lat,lon,W2',10)
switch comp
    case {'X1','X2'}
        colormap(flipud(div_red_yellow_blue11))
    case 'X3'
        colormap(seq_yellow_green_blue9)
end
coast = load('coast');
plotm(coast.lat,coast.long,'Color',zeros(1,3))
box off
colorbar
cbfreeze
grid on

% plot annotation
ylim = get(gca,'YLim');
xlim = get(gca,'XLim');
dxlim = xlim(2)-xlim(1);
dylim = ylim(2)-ylim(1);

switch comp
    case 'X1'
        textstring = '\chi_1';
    case 'X2'
        textstring = '\chi_2';
    case 'X3'
        textstring = '\chi_3';
end

text(xlim(1)+0.1*dxlim,ylim(2),textstring,'FontSize',16,'Color',zeros(1,3))
