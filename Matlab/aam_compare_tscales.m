function [t,M,O,R,G,H,MJD] = aam_compare_tscales(Tmin,Tmax,eopset,RAset)
%---plot_aam_compare.m-------------
%
% new code to compare AAM from observations, reanalysis, and emac
% started 8 april 2010
% last update 19 april 2010
% Inputs:
%       Tmin and Tmax -- the periods
%	** if Tmax is set to zero, we do no filtering at all
%       eopset -- the dataset we use for EOPs
%       RAset -- the dateset we use for reanalysis.  Choose 'ERA' or 'ECMWF'
% Outputs:
%	t = time vector in datenum values for matlab to read
%	M = filtered 3-d model vector
%	R = filtered 3-d reanalysis vector
%	G = filtered 3-d geodetic vector
%	O = filtered 3-d oceanic excitation function
%	H = filtered 3-d hydrosphere excitation function
%	* the LOD output is in seconds, and the equatorial EF output in mas

%----------------------------------------------



%---set paths
addpath('/home/ig/neef/MFiles/utilities/')
addpath('/home/ig/neef/MFiles/utilities/m_map/')


%---load or read in the relevant data
LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate
rad2mas = (180/pi)*60*60*1000;	% conversion of radians to mas

% emac output, processed.  Convert from radians to milli-arcseconds
load '/home/ig/neef/MFiles/AM/aam_mass_CCMval_1960_2000_withIB.mat'
%load  '/home/ig/neef/MFiles/AM/aam_emac_1960_2000_withIB_set3.mat'
%load '/home/ig/neef/MFiles/AM/test_aamfile.mat'
Xtemp_M = X;
MJD_M = MJD;
g_M = find(isfinite(Xtemp_M(3,:)) == 1);
Xtemp_M(3,:) = LOD0*(Xtemp_M(3,:) - mean(Xtemp_M(3,g_M)));	% (conversion to LOD)
Xtemp_M(1:2,:) = Xtemp_M(1:2,:)*rad2mas;			% (conversion to mas)

% reanalyses - read in and convert to milliarcseconds
% here there's the option option to compare the ECMWFop AAM to geodetic (2000-2010 only)
if RAset(1:3) == 'ERA'
  [Xw_EI,Xm_EI,MJD_EI] = read_EFs('aam','ERAinterim');
  [Xw_E4,Xm_E4,MJD_E4] = read_EFs('aam','ERA40');
  s = find(MJD_E4 == min(MJD_EI));
  Xtemp_E4 = Xw_E4(:,1:s-1)+Xm_E4(:,1:s-1);
  Xtemp_EI = Xw_EI+Xm_EI;
  Xtemp_R = [Xtemp_EI,Xtemp_E4];
  MJD_R = [MJD_E4(1:s-1),MJD_EI];
end
if RAset(1:3) == 'ECM'
  [Xw_EC,Xm_EC,MJD_R] = read_EFs('aam','ECMWF');
  Xtemp_R = Xw_EC+Xm_EC;
end
g_R = find(isfinite(Xtemp_R(3,:)) == 1);
Xtemp_R(1:2,:) = Xtemp_R(1:2,:)*rad2mas;			% (conversion to mas)
Xtemp_R(3,:) = LOD0*(Xtemp_R(3,:) - mean(Xtemp_R(3,g_R)));	% (conversion to LOD)

% other earth system components: unit is radians
[Xw_OAM1,Xm_OAM1,MJD_OAM1] = read_EFs('oam','ERA40');
[Xw_HAM1,Xm_HAM1,MJD_HAM1] = read_EFs('ham','ERA40');
[Xw_OAM2,Xm_OAM2,MJD_OAM2] = read_EFs('oam','ERAinterim');
[Xw_HAM2,Xm_HAM2,MJD_HAM2] = read_EFs('ham','ERAinterim');
s_O = find(MJD_OAM1 == min(MJD_OAM2));
s_H = find(MJD_HAM1 == min(MJD_HAM2));

% slap timeseries together and also convert equatorial terms to milli-arcseonds
Xtemp_O = [Xw_OAM1(:,1:s_O-1)+Xm_OAM1(:,1:s_O-1),Xw_OAM2+Xm_OAM2];
Xtemp_H = [Xw_HAM1(:,1:s_H-1)+Xm_HAM1(:,1:s_H-1),Xw_HAM2+Xm_HAM2];
MJD_O = [MJD_OAM1(1:s_O-1), MJD_OAM2];
MJD_H = [MJD_HAM1(1:s_H-1), MJD_HAM2];
Xtemp_O(1:2,:) = Xtemp_O(1:2,:)*rad2mas;			% (conversion to mas)
Xtemp_H(1:2,:) = Xtemp_H(1:2,:)*rad2mas;			% (conversion to mas)
g_O = find(isfinite(Xtemp_O(3,:)) == 1);
g_H = find(isfinite(Xtemp_H(3,:)) == 1);
Xtemp_O(3,:) = LOD0*(Xtemp_O(3,:) - mean(Xtemp_O(3,g_O)));	% (conversion to LOD)
Xtemp_H(3,:) = LOD0*(Xtemp_H(3,:) - mean(Xtemp_H(3,g_H)));	% (conversion to LOD)

% EOP-implied AAM terms: X1 and X2 in mas, dlodg in seconds
[X1g,X2g,dlodg,MJD_G] = read_eops(eopset);
dlodg2 = dlodg-mean(dlodg);
Xtemp_G = transpose([X1g, X2g, dlodg2]);


%---unify the time axes

d0 = date2mjd(1962,1,1,0,0,0);
df = date2mjd(1999,12,31,0,0,0);
if eopset(1:3) == 'JPL', d0 = date2mjd(1962,1,22,0,0,0); end
if RAset(1:3) == 'ECM' 
  d0 = date2mjd(2001,1,1,0,0,0);
  df = date2mjd(2010,2,26,0,0,0);
end
nt = df-d0+1;
MJD = d0:df;

X_M   = zeros(3,nt)+NaN;
X_R   = zeros(3,nt)+NaN;
X_G   = zeros(3,nt)+NaN;
X_O   = zeros(3,nt)+NaN;
X_H   = zeros(3,nt)+NaN;



for ii = 1:nt
  iday = MJD(ii);
  iM = find(MJD_M == iday+0.5);
  iR = find(MJD_R == iday);
  iG = find(MJD_G == iday);
  iO = find(MJD_O == iday);
  iH = find(MJD_H == iday);
  % if the days can be matched, fill in the corresponging X vectors
  if isempty(iM) == 0, X_M(:,ii) = Xtemp_M(:,iM); end
  if isempty(iR) == 0, X_R(:,ii) = Xtemp_R(:,iR);end
  if isempty(iG) == 0, X_G(:,ii) = Xtemp_G(:,iG);end
  if isempty(iO) == 0, X_O(:,ii) = Xtemp_O(:,iO);end
  if isempty(iH) == 0, X_H(:,ii) = Xtemp_H(:,iH);end
end


%---convert MJD array to a t axis that MATLAB likes
t = MJD*0;
nt = size(t,2);
[y, m, d] = mjd2date(MJD);
for ii=1:nt, t(ii)=datenum([y(ii) m(ii) d(ii)]); end



% convert to matlab timeseries

T1_M = timeseries(X_M(1,:),t,'isDatenum',true);
T2_M = timeseries(X_M(2,:),t,'isDatenum',true);
T3_M = timeseries(X_M(3,:),t,'isDatenum',true);

T1_O = timeseries(X_O(1,:),t,'isDatenum',true);
T2_O = timeseries(X_O(2,:),t,'isDatenum',true);
T3_O = timeseries(X_O(3,:),t,'isDatenum',true);

T1_R = timeseries(X_R(1,:),t,'isDatenum',true);
T2_R = timeseries(X_R(2,:),t,'isDatenum',true);
T3_R = timeseries(X_R(3,:),t,'isDatenum',true);

T1_G = timeseries(X_G(1,:),t,'isDatenum',true);
T2_G = timeseries(X_G(2,:),t,'isDatenum',true);
T3_G = timeseries(X_G(3,:),t,'isDatenum',true);

T1_H = timeseries(X_H(1,:),t,'isDatenum',true);
T2_H = timeseries(X_H(2,:),t,'isDatenum',true);
T3_H = timeseries(X_H(3,:),t,'isDatenum',true);

%if the option is set, filter out the periods of interest
if Tmax == 0
  F1_M = T1_M;	F2_M = T2_M;	F3_M = T3_M;
  F1_O = T1_O;	F2_O = T2_O;	F3_O = T3_O;
  F1_R = T1_R;	F2_R = T2_R;	F3_R = T3_R;
  F1_G = T1_G;	F2_G = T2_G;	F3_G = T3_G;
  F1_H = T1_H;	F2_H = T2_H;	F3_H = T3_H;
else
  fmax = 1/Tmin;
  fmin = 1/Tmax;
  F1_M = idealfilter(T1_M,[fmin,fmax],'pass');
  F2_M = idealfilter(T2_M,[fmin,fmax],'pass');
  F3_M = idealfilter(T3_M,[fmin,fmax],'pass');
  F1_O = idealfilter(T1_O,[fmin,fmax],'pass');
  F2_O = idealfilter(T2_O,[fmin,fmax],'pass');
  F3_O = idealfilter(T3_O,[fmin,fmax],'pass');
  F1_R = idealfilter(T1_R,[fmin,fmax],'pass');
  F2_R = idealfilter(T2_R,[fmin,fmax],'pass');
  F3_R = idealfilter(T3_R,[fmin,fmax],'pass');
  F1_G = idealfilter(T1_G,[fmin,fmax],'pass');
  F2_G = idealfilter(T2_G,[fmin,fmax],'pass');
  F3_G = idealfilter(T3_G,[fmin,fmax],'pass');
  F1_H = idealfilter(T1_H,[fmin,fmax],'pass');
  F2_H = idealfilter(T2_H,[fmin,fmax],'pass');
  F3_H = idealfilter(T3_H,[fmin,fmax],'pass');
end

M = zeros(3,nt);
O = zeros(3,nt);
R = zeros(3,nt);
G = zeros(3,nt);
H = zeros(3,nt);

gm = find(isfinite(X_M(1,:)) == 1);
go = find(isfinite(X_O(1,:)) == 1);
gr = find(isfinite(X_R(1,:)) == 1);
gg = find(isfinite(X_G(1,:)) == 1);
gh = find(isfinite(X_H(1,:)) == 1);

M(1,:) = F1_M.data;
M(2,:) = F2_M.data;
M(3,:) = F3_M.data;

O(1,:) = F1_O.data;
O(2,:) = F2_O.data;
O(3,:) = F3_O.data;

R(1,:) = F1_R.data;
R(2,:) = F2_R.data;
R(3,:) = F3_R.data;

G(1,:) = F1_G.data;
G(2,:) = F2_G.data;
G(3,:) = F3_G.data;

H(1,:) = F1_H.data;
H(2,:) = F2_H.data;
H(3,:) = F3_H.data;

end
