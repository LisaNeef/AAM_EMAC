% ----read_EFs.m----------------
%
%  Read the OAM, AAM, and HAM produced by models (OMCT and LSDM) coupled to reanalysis 
%    data for the atmosphere
%  MODS:
%    28 Jan 2010: enable reading OAM from ECCO model (notes v4, p 81)
%
%  function [Xw,Xm,mjd] = read_EFs(AM,RA)
%
%  Inputs:
%	AM: AM type - choose oam, aam, or ham
%	RA: reanalysis set (or model) name.  Choices:
%	  ERA40, ERAinterim, ECMWF
%	  ECCO (Ocean only)
%	dailvalues - set to 1 to only choose once a day values
%		...this one only works if supdaily data are available,
%		but not for obs that are less frequent than daily (e.g. ECCO)
%  Outputs:
%	[Xw_out,Xm_out,mjd_out]
%	terms 1:2 in Xw_out and Xm_out are in radians.
%	term 3 is dimensionless.
function [Xw_out,Xm_out,mjd_out] = read_EFs(AM,RA,dailyvalues)


%--file paths and stuff

% flags.
era40 = 0;
ecmwf = 0;
eraint= 0;
ecco= 0;

%--specify start / stop years available for each scenario
if RA(1:4) == 'ECMW' 
  y0 = 2000;
  yf = 2009;
  RAdir = 'IERS_opECMWF'; 
  RA2 = 'opECMWF';
  ecmwf = 1;
  obs_per_day = 4;
end
if RA(1:4) == 'ERA4' 
  y0 = 1958;
  yf = 2001;
  RAdir = 'IERS_ERA40'; 
  RA2 = RA;
  era40 = 1;
  obs_per_day = 4;
end
if RA(1:4) == 'ERAi' 
  y0 = 1989;
  yf = 2008;
  RAdir = 'IERS_ERAinterim'; 
  RA2 = RA;
  eraint = 1;
  obs_per_day = 4;
end
if RA(1:4) == 'ECCO' 
  y0 = 1993;
  yf = 2003;
  RAdir = 'IERS-ECCO'; 
  RA2 = RA;
  ecco = 1;
  obs_per_day = 1;
end


datadir = ['/dsk/nathan/lisa/' , RAdir];



%--initialize arrays
nyears = yf-y0;
nt1 = nyears*366*obs_per_day;

year = zeros(1,nt1)+NaN;
month = zeros(1,nt1)+NaN;
day = zeros(1,nt1)+NaN;
hour = zeros(1,nt1)+NaN;
mjd = zeros(1,nt1)+NaN;
Xm  = zeros(3,nt1)+NaN;
Xw = zeros(3,nt1)+NaN;


%--cycle through years and read in the files
if eraint | era40 | ecmwf
  k1 = 1;
  for iy = y0:yf
    fname = [datadir '/'  RA2 '.' num2str(iy) '.' AM];
   
    if exist(fname) == 0, disp(['Cant find reanalysis file  ', fname]), end
    % specify nr of header lines.
    if AM == 'oam', nh = 43; end 
    if AM == 'aam', nh = 36; end
    if AM == 'ham', nh = 48; end
    if exist(fname) == 2
      dum = importdata(fname,' ',nh);
      k2 = k1-1+size(dum.data,1);
      year(k1:k2) 	= dum.data(:,1); 
      month(k1:k2) 	= dum.data(:,2); 
      day(k1:k2) 		= dum.data(:,3); 
      hour(k1:k2) 	= dum.data(:,4); 
      mjd(k1:k2) 		= dum.data(:,5); 
      Xm(1,k1:k2) 	= dum.data(:,6); 
      Xm(2,k1:k2) 	= dum.data(:,7); 
      Xm(3,k1:k2) 	= dum.data(:,8); 
      Xw(1,k1:k2) 	= dum.data(:,9); 
      Xw(2,k1:k2) 	= dum.data(:,10); 
      Xw(3,k1:k2) 	= dum.data(:,11); 
      k1 = k2+1;
    end
  end
end

if ecco
  % for ecco, if daily values is not selected, use the 10 year averages instead.
  if dailyvalues
    fname = [datadir,'/ECCO_kf049f.oam']; 
    disp('Using ECCO daily values, short data set')
  else
    fname = [datadir,'/ECCO_50yr.oam']; 
    disp('Using ECCO Weekyl Averages, 50 year data set')
  end
  
  if exist(fname) == 0, disp(['Cant find reanalysis file  ', fname]), return, end
  nh = 42;
  dum = importdata(fname,' ',nh);

  % for not daily values, ECCO OAM data are 10-d averages, so need to reinitialize arrays
  % note also that for ECCO OAM, unit is kg-m**2/s: nondimensionalize (to rad, rad/s)
  mjd			= dum.data(:,1);
  Xm  = zeros(3,length(mjd));
  Xw  = zeros(3,length(mjd));

  % prefactors (see Gross 09 and notes, vol. 4 p. 82)
  aam_constants_gross
  Re = Re_m;      % (use earth radius in meters)
  pw12 = 1.608/(CminusA*Q)
  pm12 = 1.608*0.684/(CminusA*Q);
  pw3  = 0.997/(Cm*Q);
  pm3  = 0.997*0.750/(Cm*Q);
  
  Xm(1,:)		= pm12*dum.data(:,2);
  Xm(2,:)		= pm12*dum.data(:,3);
  Xm(3,:)		= pm3*dum.data(:,4);
  Xw(1,:)		= pw12*dum.data(:,5);
  Xw(2,:)		= pw12*dum.data(:,6);
  Xw(3,:)		= pw3*dum.data(:,7);

end


test_mjd = mjd - round(mjd);
round_mjd = find(test_mjd == 0);

if dailyvalues 
  mjd_out = mjd(round_mjd);
  Xm_out = Xm(:,round_mjd);
  Xw_out = Xw(:,round_mjd);
else
  mjd_out = mjd;
  Xm_out = Xm;
  Xw_out = Xw;
end
