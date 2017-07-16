function W = eam_weights(lat,lon,comp,variable)
% 
% Compute the approproate geometric weights for the dynamical variables U,
% V, and PS, needed to compute the effective angular momentum (EAM)
% functions X1, X2, and X3, given a certain lat/lon grid.
%
%
% INPUTS:
%   lat: latitude grid in degrees.
%   lon: longitude grid in degrees
%   comp: angular momentum component: 'X1', 'X2', or 'X3'
%   var: variable for which weighting is needed: 'U','V', or 'PS'
%
% MODS:
%   6 Mar 2011: for X3, variabe 'V' should be weighted with zeros, since it
%   doesn't appear in the integral.  (Previously had it weighted as 1
%   ...woops!)
%
%  Lisa Neef / 29 Nov 2011.

%% convert lat and lon to radians & make a meshgrid.
rlon = lon*pi/180;
rlat = lat*pi/180;

[LAT,LON] = meshgrid(rlat,rlon);

%% compute the appropriate weighting functions

switch comp
    case 'X1'
        switch variable
            case 'U'
                W = sin(LAT).*cos(LAT).*cos(LON);
            case 'V'
                W = -cos(LAT).*sin(LON);                          
            case 'PS'
                W = sin(LAT).*cos(LAT).*cos(LAT).*cos(LON);
        end
    case 'X2'
        switch variable
            case 'U'
                W = sin(LAT).*cos(LAT).*sin(LON);
            case 'V'
                W = cos(LAT).*cos(LON);
            case 'PS'
                W = sin(LAT).*cos(LAT).*cos(LAT).*sin(LON);
        end        
    case 'X3'
        switch variable
            case 'U'
                W = cos(LAT).^2;
            case 'V'
                W = zeros(size(LAT));                
            case 'PS'
                W = cos(LAT).^3;
        end
    case 'none'
        W = ones(size(LAT));
end




