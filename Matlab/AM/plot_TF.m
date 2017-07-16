%% plot_TF_polar_stereogr.m
function plot_TF_polar_stereogr(comp,variable)
%
% plot_TF_polar_stereogr.m
% make plots of the AAM excitation functions in mercator projection
% Lisa Neef, 13 April 2012
%
% INPUTS:
%   comp: angular momentum component (X1,X2,X3)
%   variable: 'U', 'V', or 'PS'


%% temp inputs
%comp = 'X1';
%variable = 'V';

%% retrieve the desired excitation function

lat = -90:.5:90;
lon = -180:.5:180;

W = eam_weights(lat,lon,comp,variable);

[X,Y] = meshgrid(lon,lat);


%% make the polar stereographic plot.

% select the colormap - diverging or progression

switch comp
    case {'X1','X2'}
       cmap = flipud(div_red_yellow_blue11);
    case 'X3'
       cmap1 = flipud(div_red_yellow_blue11);
       cmap = cmap1(6:11,:);
end


%axesm ('MapProjection','stereo', 'Frame', 'on', 'Grid', 'on','MapLatLimit',[0 90]);
contourf(X,Y,W',10)
colormap(cmap)
coast = load('coast');
hold on
plot(coast.long,coast.lat,'Color',zeros(1,3))
box off
colorbar
grid on
freezeColors

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

title(textstring)
%text(xlim(1)+0.1*dxlim,ylim(2),textstring,'FontSize',16,'Color',zeros(1,3))
