% ----read_eops.m----------------
function [X1_out,X2_out,dlod_out,mjd_out,eX1_out,eX2_out,edlod_out] = read_eops
%
%  Read the earth orientation parameters from IERS (eopc04.62) into matlab, 
%    compute the implied "geodetic" excitation functions and delta-LOD, and
%    return a time series
%
%  Mods:
%   19 March 2010 to deal with files from http://hpiers.obspm.fr/eop-pc/
%    which don't have tides.
%   8 April 2010: Instead of my 1/2-assed digital filter (which doesn't really work) use canned filtering 
%	thing to take out the chandler wobble.  Also make the output in mas instead of radians.
%   14 April 2010: fixed the error that p=x+iy (actually it's p = x-iy), see notes vol 2 p73 and Nathan blog)
%       so now don't need the canned filter anymore.
%   15 April 2010: add option to use COMB2000 data set from JPL (see, eg Gross et al 2003)
%   10 Nov 2010: make dlod output in ms, and add option to only have daily values
%   19 Nov 2010: also return the uncertainties here.
%   11 Mar 2011: fixed small bug in the reading of IERS2
%   21 Mar 2011: downloaded new data set encompassing 1962-2008, for complete continuity (notes v4, p123)
%
%
%  OUTPUTS:
%	X1, X2: equatorial excitation functions in mas
%	dlod: LOD fluctations in ms
%	mjd: mod julian day.



%--file paths and stuff
  dir = '/dsk/nathan/lisa/IERS-ERP/';
  d = 1;
  ff = 'C04_1962_2008_notides.txt'
  disp('Reading IERS dataset 1962-2008')

fname=[dir,ff];

% set default of dailyvalues
dailyvalues = 1;

%--constants!

sigc = 2*pi/433;	% chandler frequency in rad/days
Tc = 431.2;		% chandler period in days
Q = 179;		% chandler resonance quality factor
sigc = 2*pi/Tc*(1+i/(2*Q));	% chandler frequency in radians/day

%--read in the file
eop = importdata(fname,' ',2);

%--select the important arrays
if d == 1
  mjd	= eop.data(:,1);
  x 	= eop.data(:,2);	% PM-x (mas)
  ex 	= eop.data(:,3);	% error in PM-x (mas)
  y 	= eop.data(:,4);	% PM-y (mas)
  ey 	= eop.data(:,5);	% error in PM-y (mas)
  dlod_in	= eop.data(:,10);	% LOD change (ms)
  edlod 	= eop.data(:,11);	% error in LOD change 
end


if d == 2
  mjd	= eop.data(:,1);
  x 	= eop.data(:,2)*1e-3;	% PM-x (mas)
  y 	= eop.data(:,3)*1e-3;	% PM-y (mas)
  dlod_in	= eop.data(:,5)*1e-4;	% LOD change (ms)
end


%--convert polar motion terms into X1 and X2:
% where X(t) = p(t)+(i/sigma)*deriv(p(t))
p = x-i*y;		% note the negative here -- see Gross et al 1996a
ep = ex-i*ey;		% note the negative here -- see Gross et al 1996a
pdot = p*0;
epdot = p*0;
h = 10+0*i;
for k = real(h)+1:length(p)-real(h)
  pdot(k) = (p(k+h)-p(k-h))/(2*h);
  epdot(k) = (ep(k+h)-ep(k-h))/(2*h);
end

X = p+(i/sigc)*pdot;
%***NEED TO MAKE SURE THIS ERROR MAPPING IS CORRECT
eX = ep+(i/sigc)*epdot;

X1 = real(X);
X2 = imag(X);
eX1 = real(eX);
eX2 = imag(eX);

if dailyvalues
  mjd_test = mjd- round(mjd);
  mjd_out = mjd(find(mjd_test == 0));
  X1_out = X1(find(mjd_test == 0));
  X2_out = X2(find(mjd_test == 0));
  dlod_out = dlod_in(find(mjd_test == 0));
  eX1_out = eX1(find(mjd_test == 0));
  eX2_out = eX2(find(mjd_test == 0));
  edlod_out = edlod(find(mjd_test == 0));
else
  mjd_out = mjd;
  X1_out = X1;
  X2_out = X2;
  dlod_out = dlod_in;
  eX1_out = eX1;
  eX2_out = eX2;
  edlod_out = edlod;
end

