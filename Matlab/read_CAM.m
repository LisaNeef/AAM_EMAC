%----read_CAM.m----------------
%
% This function reads in the CAM data from an IERS file so that
% it can be subtracted from the EOPs.
% syntax:
% [mjd,dLOD] = read_CAM(inventory)
%
% INPUT: choose which dataset to look at
%  choose either 'pais' or 'jack'
%  (there are also others but they have weird formats that I haven't dealt with yet)
% 
%  Started 8 March 2010.
function [mjd,dLOD] = read_CAM(inventory)

%--file paths and stuff

dir = '/dsk/nathan/lisa/Data/IERS-SBC/';
if inventory(1:4) == 'pais' 
  ff = 'pais_table.txt';
  nh = 1;
end
if inventory(1:4) == 'jack'
  ff = 'jackson_table.txt'; 
  nh = 1;
end
fname=[dir,ff];



%--read in the file
cam = importdata(fname,' ',nh);

year = cam.data(:,1);
dLOD = cam.data(:,2);


%--convert the output year into MJD for easier comparison to other data.

mjd = date2mjd(year);

mjd(find(year == 0)) = NaN;
dLOD(find(year == 0)) = NaN;

end
