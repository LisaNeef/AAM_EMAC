% Definitions of the various constants used in AAM computations,
% all taken from the Gross, 2009, "Earth Variations - Long Period"
%
% This script should be run in MATLAB to get all the right
% constants onto the desktop.
%
% THIS NEEDS TO BE FINISHED!
% Updates:
%  6 Jul 2011: correction to LOD0 -- should be sidereal day, not solar day.
%---------------------------------------------------------------------------------



M_atm	=  5.1441e18;		% mass of the atmosphere in kg
M_oc	 =  1.4e21;  		% mass of the ocena in kg

% observed whole-earth parameters

Q	=  7.292115e-5;   	% rot. rate of the earth in rad/s
M	=  5.9737e24;     	% mass of the earth in kg
C	=  8.0365e37;     	% axial principal moment of inertia  (kg m2)
B	=  8.0103e37;		% next-largest principal MOI (kg m2)
A	=  8.0101e37;		% next-largest principal MOI (kg m2)
CminusA	=  2.6398e35;		% kg m2
CminusB	=  2.6221e35;		% kg m2
BminusA	=  1.763e33;		% kg m2

% modeled whole-earth parameters

Re_m	= 6.371e6;		% radius of earth (m)
Re_km	= 6371.0;		% radius of earth (km)
g 	= 9.81;               	% grav constant in m/s2


% crust and mantle parameters
Mm	= 4.0337e24;		% mass of mantle (km)
Cm 	= 7.1236e37;		% principal MOI of mantle (kgm^2)
Am	= 7.0999e37;		% next-largest MOI (kgm^2)

% other conversions, etc.

d2r = pi/180.;        		 %  degrees to radian
rad2microas = (180/pi)*60*60*1e6;%  radians to micro arcseconds
rad2mas = (180/pi)*60*60*1e3;	 %  radians to milli arcseconds
LOD0_ms = double(86160*1e3);     % sidereal LOD in milliseconds.

% these are the EAF prefactors.


