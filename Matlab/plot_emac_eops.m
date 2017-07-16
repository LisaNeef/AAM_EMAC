% plot a comparison between EMAC AAM components and reanalysis / EOPS over certain latitude bands
% mod 8 April 2010: change units to mas (or ms for LOD) and change selection of period band
%
% Inputs:----
%	Tmin = lowest period considered
%	Tmax= highest period considered
function plot_emac_eops(Tmin,Tmax)


%---set paths
addpath('/home/ig/neef/MFiles/utilities/')

% conversion constants
LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate
rad2mas = (180/pi)*60*60*1000;	% conversion of radians to mas


%---load saved EMAC AAM computations and convert eq components to mas
load  '/home/ig/neef/MFiles/AM/aam_emac_1960_2000_withIB_set3.mat'
% rename the equatorial excitation functions:
Xtemp_M = X;
mjd_AAM = MJD;

%---read in EOPS 

[X1g,X2g,dlodg,MJDg] = read_eops;
Xg = transpose([X1g, X2g, dlodg]);

%---read in AM from other sources
[Xw_O1,Xm_O1,mjd_O1] = read_EFs('oam','ERA40');
[Xw_H1,Xm_H1,mjd_H1] = read_EFs('ham','ERA40');
[Xw_R1,Xm_R1,mjd_R1] = read_EFs('aam','ERA40');
[Xw_O2,Xm_O2,mjd_O2] = read_EFs('oam','ERAinterim');
[Xw_H2,Xm_H2,mjd_H2] = read_EFs('ham','ERAinterim');
[Xw_R2,Xm_R2,mjd_R2] = read_EFs('aam','ERA40');


%--unify them all on the same t axis
d0 = date2mjd(1962,1,2,0,0,0);
df = date2mjd(1999,12,31,0,0,0);
nt = df-d0+1;
MJD = d0:df;

Xo = zeros(3,nt)+NaN;
Xh = zeros(3,nt)+NaN;
Xg = zeros(3,nt)+NaN;
Xm = zeros(3,nt)+NaN;
Xr = zeros(3,nt)+NaN;

for ii=1:nt
  fo1= find(MJD(ii) == mjd_O1);
  fo2 = find(MJD(ii) == mjd_O2);
  fh1= find(MJD(ii) == mjd_H1);
  fh2 = find(MJD(ii) == mjd_H2);
  fr1= find(MJD(ii) == mjd_R1);
  fr2 = find(MJD(ii) == mjd_R2);
  fg = find(MJD(ii) == mjd_g);
  fm = find(MJD(ii)+0.5 == mjd_M);
  if isempty(fo1) == 0, Xo(:,ii) = Xw_O1(:,f1)+Xm_O1(:,f1); end
  if isempty(fo2) == 0, Xo(:,ii) = Xw_O2(:,f1)+Xm_O2(:,f1); end
  if isempty(f1b) == 0, Xoam(:,ii) = Xw_OAM2(:,f1b)+Xm_OAM2(:,f1b); end
  if isempty(f2) == 0, Xham(:,ii) = Xw_HAM1(:,f2)+Xm_HAM1(:,f2); end
  if isempty(f2b) == 0, Xham(:,ii) = Xw_HAM2(:,f2b)+Xm_HAM2(:,f2b); end
  if isempty(f3) == 0, Xgeo(:,ii) = Xg(:,f3); end
  if isempty(f4) == 0, Xaam(:,ii) = Xtemp_M(:,f4); end
end

% subtract other stuff from geodetic timeseries
% FOR NOW: filtering of tides and cam happens below.
Xnet = Xgeo - Xham - Xoam;


% compute implied dlod in each case
LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate
g1 = find(isfinite(Xaam(3,:)) == 1); 
g2 = find(isfinite(Xoam(3,:)) == 1); 
g3 = find(isfinite(Xham(3,:)) == 1); 
g4 = find(isfinite(Xgeo(3,:)) == 1); 

dlod_AAM = (Xaam(3,:)-mean(Xaam(3,g1)))*LOD0;
dlod_OAM = (Xoam(3,:)-mean(Xoam(3,g2)))*LOD0;
dlod_HAM = (Xham(3,:)-mean(Xham(3,g3)))*LOD0;
dlod_GEO = Xgeo(3,:);
dlod_NET = Xnet(3,:);

%---also compute deviations from the mean in t
dXaam = Xaam - mean(Xaam(:,g1),2)*ones(1,nt);
dXoam = Xoam - mean(Xoam(:,g2),2)*ones(1,nt);
dXham = Xham - mean(Xham(:,g3),2)*ones(1,nt);
dXgeo = Xgeo - mean(Xgeo(:,g4),2)*ones(1,nt);
dXnet = Xnet - mean(Xnet(:,g4),2)*ones(1,nt);


%---convert MJD array to a t axis that MATLAB likes
t = MJD*0+NaN;
[y, m, d] = mjd2date(MJD);
for ii=1:nt, t(ii)=datenum([y(ii) m(ii) d(ii)]); end

%---select timescales on which to focus and filter from tseries.
nb = 1;
Tband = zeros(7,2);
Tband(1,:) = [1,14];            % daily to weekly
Tband(1,:) = [7,14];            % daily to weekly
Tband(2,:) = [20,40];           % Monthly Variations
Tband(3,:) = [60,120];          % Intraseasonal
Tband(4,:) = [120,300];         % Interseasonal
Tband(5,:) = [360 370];         % Annual cycle
Tband(6,:) = [600,900];         % QBO
Tband(7,:) = [1, 10288];        % everything

Tb = Tband(bandno,:);

fband = ones(7,2)./Tband;       % frequency bands in 1/day

% for each component, create a timeseries to filter around
ts_X1aam =timeseries(Xaam(1,:),t,'Name','X1aam','isDatenum',true);
ts_X1geo =timeseries(Xgeo(1,:),t,'Name','X1geo','isDatenum',true);
ts_X1netA =timeseries(Xnet(1,:),t,'Name','X1net','isDatenum',true);

ts_X2aam =timeseries(Xaam(2,:),t,'Name','X2aam','isDatenum',true);
ts_X2geo =timeseries(Xgeo(2,:),t,'Name','X2geo','isDatenum',true);
ts_X2netA =timeseries(Xnet(2,:),t,'Name','X2net','isDatenum',true);

ts_DLaam =timeseries(dlod_AAM,t,'Name','X3aam','isDatenum',true);
ts_DLgeo =timeseries(dlod_GEO,t,'Name','X3geo','isDatenum',true);
ts_DLnetA =timeseries(dlod_NET,t,'Name','X3net','isDatenum',true);

% the "net" timeseries is additionally freed of tides and CAM
tide1 = [1/9, 1/10];
tide2 = [1/13, 1/14];
tide3 = [1/182, 1/183];

ts_X1netB = idealfilter(ts_X1netA,tide1,'notch');
ts_X1netC = idealfilter(ts_X1netB,tide2,'notch');
ts_X1netD = idealfilter(ts_X1netC,tide3,'notch');
ts_X1net = ts_X1netD;
ts_X2netB = idealfilter(ts_X2netA,tide1,'notch');
ts_X2netC = idealfilter(ts_X2netB,tide2,'notch');
ts_X2netD = idealfilter(ts_X2netC,tide3,'notch');
ts_X2net = ts_X2netD;
ts_DLnetB = idealfilter(ts_DLnetA,tide1,'notch');
ts_DLnetC = idealfilter(ts_DLnetB,tide2,'notch');
ts_DLnetD = idealfilter(ts_DLnetC,tide3,'notch');
ts_DLnet = ts_DLnetD;

% set up the matrix that will hold the correlations.
% p = [X component, lag time]
p = zeros(3,2*nt-1);

% axes for each time window
ax = zeros(nb,2);
dT = Tb(2);
ax= [1,20*dT];

% filter out different frequency bands, compute correlations, and plot.

col_GEO = [0.392,0.6555,0.1712];
col_AAM = 0*ones(3,1);
col_NET = [0.917,0.2858,0.7572];
col_res = [1,0,0];


ib = bandno;
  fX1aam = idealfilter(ts_X1aam,fband(ib,:),'pass');
  fX1geo = idealfilter(ts_X1geo,fband(ib,:),'pass');
  fX1net = idealfilter(ts_X1net,fband(ib,:),'pass');
  fX2aam = idealfilter(ts_X2aam,fband(ib,:),'pass');
  fX2geo = idealfilter(ts_X2geo,fband(ib,:),'pass');
  fX2net = idealfilter(ts_X2net,fband(ib,:),'pass');
  fDLaam = idealfilter(ts_DLaam,fband(ib,:),'pass');
  fDLgeo = idealfilter(ts_DLgeo,fband(ib,:),'pass');
  fDLnet = idealfilter(ts_DLnet,fband(ib,:),'pass');

  % compute time-lagged correlations for each one.
  [c1,lags] = xcorr(fX1aam.data,fX1net.data,'coeff');
  [c2,lags] = xcorr(fX2aam.data,fX2net.data,'coeff');
  [c3,lags] = xcorr(fDLaam.data,fDLnet.data,'coeff');

  p(1,:) = c1;
  p(2,:) = c2;
  p(3,:) = c3;

  % plot!
  figure(2),clf
    if ib < 7
      subplot(4,1,1)
    else
      subplot(3,1,1)
    end
    AAMplot= plot(fX1aam,'Color',col_AAM);
    hold on
    GEOplot= plot(fX1geo,'Color',col_GEO);
    NETplot= plot(fX1net,'Color',col_NET);
    lhandle = [AAMplot(1)  GEOplot(1) NETplot(1)];
    legend(lhandle,'EMAC+IB','Geod.','Geod.-HAM-OAM' )
    xlabel('Time')
    ylabel(['X_1'])
    [a,b] = max(squeeze(p(1,:)));
    title(['Period Band  ' num2str(Tband(ib,1)) ' - ' num2str(Tband(ib,2)) ' days'])
    bot = min(min(fX1aam.data, fX1geo.data));
    top = max(max(fX1aam.data,fX1geo.data));
    if ib < 7
      axis([ax bot top])
    end
    disp(['Period Band  ' num2str(Tband(ib,1)) ' - ' num2str(Tband(ib,2)) ' days'])
    disp(['    X1 \rho = ' num2str(a) '  lag = ' num2str(lags(b))])


    if ib < 7
      subplot(4,1,2)
    else
      subplot(3,1,2)
    end
    AAMplot= plot(fX2aam,'Color',col_AAM);
    hold on
    GEOplot= plot(fX2geo,'Color',col_GEO);
    NETplot= plot(fX2net,'Color',col_NET);
    xlabel('Time')
    ylabel(['X_2'])
    [a,b] = max(squeeze(p(2,:)));
    bot = min(min(fX2aam.data, fX2geo.data));
    top = max(max(fX2aam.data,fX2geo.data));
    if ib <7 
      axis([ax bot top])
    end
    disp(['    X2 \rho = ' num2str(a) '  lag = ' num2str(lags(b))])

    if ib < 7
      subplot(4,1,3)
    else
      subplot(3,1,3)
    end
    AAMplot= plot(fDLaam,'Color',col_AAM);
    hold on
    GEOplot= plot(fDLgeo,'Color',col_GEO);
    NETplot= plot(fDLnet,'Color',col_NET);
    xlabel('Time')
    ylabel(['\Delta LOD'])
    [a,b] = max(squeeze(p(3,:)));
    bot = min(min(fDLaam.data, fDLgeo.data));
    top = max(max(fDLaam.data,fDLgeo.data));
    if ib < 7
      axis([ax bot top])
    end
    disp(['    dLOD \rho = ' num2str(a) '  lag = ' num2str(lags(b))])

    if ib < 7
    subplot(4,1,4)
    plot(lags,squeeze(c1),lags,squeeze(c2),lags,squeeze(c3))
    xlabel('lag (days)')
    ylabel('correlation')
    axis([-10*dT 10*dT -1 1])
    end


end










