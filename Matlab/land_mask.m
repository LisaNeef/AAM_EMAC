%% this function builds a land mask for a given grid

function [out] = land_mask(lon,lat,detail)
% variables:
%           lon:    vector containing all longitudes 
%           lat:    vector containing all latitudes from -90[S]:90[N]
%           detail: detail level [1:537] --> [worst:best]

% By Christof Petrick, 31.03.2010

% lon = [-15:0.1:30];
% lat = [30:0.1:70];
% detail = 500

lon(lon>180) = lon(lon>180)-360;


[LONGITUDE,LATITUDE] = meshgrid(lon,lat);

s=shaperead('landareas.SHP','UseGeoCoords',true);

% add for antarctic 
    i = 1
    X = s(i,1).Lon;
    Y = s(i,1).Lat;
    X = [X,180,-180,-180];
    Y = [Y,-90.5,-90.5,Y(1,1)];
    out = inpolygon(LONGITUDE,LATITUDE,X,Y);
    
for i = 2:detail
    X = s(i,1).Lon;
    Y = s(i,1).Lat;
    in = inpolygon(LONGITUDE,LATITUDE,X,Y+0.01);
    out = out + in;
    i
end



figure
axesm('gstereo','Frame', 'on','Grid','on',...
    'MapLonLimit',[-180 180],'MapLatLimit',[-90 90],'MeridianLabel','on','ParallelLabel','on')
geoshow(LATITUDE,LONGITUDE,out,'DisplayType', 'Texturemap')
axis tight
hold on 
c =load('coast')
plotm(c.lat,c.long,'k')


