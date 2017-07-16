function [X1,X2,dlod,mjd,date] = read_CAM_JH(topog,cond)
% function [X1,X2,dlod,mjd,date] = read_CAM_JH(topog,cond)
%
% Read in the CAM timeseries produced by Jan Hagedoorn.  
% 25 Nov 2010 (see notes vol. 4 p. 33)
% 
%  INPUTS:
%	topog: select topographic model used: 's' or 't'
%       cond: select conductivity profile used: 1 or 4.
%	
%
%  OUTPUTS:
%	X1, X2: equatorial excitation functions in mas
%	dlod: LOD fluctations in ms
%	mjd: mod julian day.

%--file paths and stuff
datadir = '/dsk/nathan/lisa/CAM_JanH/';
ff = ['exf_',topog,num2str(cond),'.dat'];
fname=[datadir,ff];
exist(fname);

%--read in the file
% convert chi1 and chi2 from radians to mas
aam_constants_gross
cam = importdata(fname,' ',0);

date = squeeze(cam(:,1));
X1 = rad2mas*squeeze(cam(:,2));
X2 = rad2mas*squeeze(cam(:,3));
dlod = LOD0_ms*squeeze(cam(:,4));

% convert date to MJD
year = floor(date);
extra_days = (date-year)*365;
mjd = date2mjd(year,1,1,0,0,0)+extra_days;

