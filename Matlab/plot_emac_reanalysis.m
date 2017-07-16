function plot_emac_reanalyis(Tmin,Tmax)
%
%  Compare EMAC and reanalysis at different frequency bands
%  mods: 
%   24 March 2010 ...take out the IB data sets
%   7 April 2010: do it so that we can specify the period band beforehand
%
%  TO DO:
%	- make it so that we can put in any frequency instead of selecting the specific band.


%---set paths
addpath('/home/ig/neef/MFiles/utilities/')
addpath('/home/ig/neef/MFiles/utilities/m_map/')

LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate

% emac output, processed.
load  '/home/ig/neef/MFiles/AM/aam_emac_1960_2000_withIB_set3.mat'
Xtemp_M = X;
MJD_M = MJD;

% reanalyses
[Xw_EI,Xm_EI,MJD_EI] = read_EFs('aam','ERAinterim');
[Xw_E4,Xm_E4,MJD_E4] = read_EFs('aam','ERA40');
s = find(MJD_E4 == min(MJD_EI));
Xtemp_E4 = Xw_E4(:,1:s-1)+Xm_E4(:,1:s-1);
Xtemp_EI = Xw_EI+Xm_EI;
Xtemp_R = [Xtemp_EI,Xtemp_E4];
MJD_R = [MJD_E4(1:s-1),MJD_EI];
g_R = find(isfinite(Xtemp_R(3,:)) == 1);

%---unify the datasets on a single time axis, putting ERAinterim and ERA40 together.
% (set up such that if ERAinterim data exists, that's the one that is used)

d0 = date2mjd(1960,1,1,0,0,0);
df = date2mjd(1999,12,31,0,0,0);
nt = df-d0+1;
MJD = d0:df;
X_M = zeros(3,nt)+NaN;
X_R = zeros(3,nt)+NaN;

for ii = 1:nt
  iday = MJD(ii);
  iM = find(MJD_M == iday+0.5);
  iR = find(MJD_R == iday);
  % if the days can be matched, fill in the corresponging X vectors
  if isempty(iM) == 0, X_M(:,ii) = Xtemp_M(:,iM); end
  if isempty(iR) == 0, X_R(:,ii) = Xtemp_R(:,iR);end
end

g = find(isfinite(X_M(1,:)) == 1);
dl_M = LOD0*(X_M(3,:) - mean(X_M(3,g)));      % (conversion to LOD)
dl_R = LOD0*(X_R(3,:) - mean(X_R(3,g)));      % (conversion to LOD)

%---convert MJD array to a time axis that MATLAB likes
t = MJD*0+NaN;
[y, m, d] = mjd2date(MJD);
for ii=1:nt, t(ii)=datenum([y(ii) m(ii) d(ii)]); end

%---select timescales on which to focus and filter from tseries.
nb = 7;		% number of frequency bands to look at

Tband = [Tmin, Tmax]
fband = ones(1,2)./Tband;	% frequency bands in 1/day


% for each component, create a timeseries to filter around
ts_X1_R =timeseries(X_R(1,:),t,'Name','X1R','isDatenum',true);
ts_X1_M =timeseries(X_M(1,:),t,'Name','X1M','isDatenum',true);

ts_X2_M =timeseries(X_M(2,:),t,'Name','X2M','isDatenum',true);
ts_X2_R =timeseries(X_R(2,:),t,'Name','X2R','isDatenum',true);

ts_DL_M =timeseries(dl_M,t,'Name','dLOD EMAC +  IB','isDatenum',true);
ts_DL_R =timeseries(dl_R,t,'Name','dLOD ERA-40/Interim','isDatenum',true);

% filter out different frequency bands, compute correlations, and plot.

col_R = [0.917,0.2858,0.7572];
col_M = 0*ones(3,1);
col_M2 = 0.7*ones(3,1);
col_res = [1,0,0];

% set up the matrix that will hold the correlations.
% p = [X component,  lag time]
p = zeros(3,2*nt-1);

% axes for each time window
dT = Tband(2)-Tband(1);
ax = [1,20*dT]

  fX1_M = idealfilter(ts_X1_M,fband,'pass');  
  fX2_M = idealfilter(ts_X2_M,fband,'pass');  
  fDL_M = idealfilter(ts_DL_M,fband,'pass');  
  fX1_R = idealfilter(ts_X1_R,fband,'pass'); 
  fX2_R = idealfilter(ts_X2_R,fband,'pass');
  fDL_R = idealfilter(ts_DL_R,fband,'pass'); 


  % compute correlations for this timescale
  [c1,lags] = xcorr(fX1_M.data,fX1_R.data,'coeff');
  [c2,lags] = xcorr(fX2_M.data,fX2_R.data,'coeff');
  [c3,lags] = xcorr(fDL_M.data,fDL_R.data,'coeff');
  
  p(1,:) = c1;
  p(2,:) = c2;
  p(3,:) = c3;


  figure(1),clf
      subplot(4,1,1)
    Mplot= plot(fX1_M,'Color',col_M);
    hold on
    ERplot = plot(fX1_R,'Color',col_R);
    lhandle = [Mplot(1)  ERplot(1) ];
    legend(lhandle,'EMAC+IB','ERAinterim / ERA-40')
    xlabel('Time')
    ylabel(['X_1'])
    [a,b] = max(squeeze(p(1,:)));
    title(['Period Band  ' num2str(Tband(1)) ' - ' num2str(Tband(2)) ' days'])
    bot = min(min(fX1_M.data, fX1_R.data));
    top = max(max(fX1_M.data,fX1_R.data));
    %axis([ax(1) ax(2) bot top])
    disp(['Period Band  ' num2str(Tband(1)) ' - ' num2str(Tband(2)) ' days'])
    disp(['    X1 \rho = ' num2str(a) '  lag = ' num2str(lags(b))])

      subplot(4,1,2)
    Mplot= plot(fX2_M,'Color',col_M);
    hold on
    ERplot = plot(fX2_R,'Color',col_R);
    xlabel('Time')
    ylabel(['X_2'])
    [a,b] = max(squeeze(p(2,:)));
    bot = min(min(fX2_M.data, fX2_R.data));
    top = max(max(fX2_M.data,fX2_R.data));
    %axis([ax(1) ax(2) bot top])
    text(0.2*top,0.2*ax(2),['\rho = ' num2str(a) '  lag = ' num2str(lags(b))])
    disp(['    X2 \rho = ' num2str(a) '  lag = ' num2str(lags(b))])
 
      subplot(4,1,3)
    Mplot = plot(fDL_M,'Color',col_M);
    hold on
    EIplot = plot(fDL_R,'Color',col_R);
    xlabel('Time')
    ylabel('\Delta LOD')
    [a,b] = max(squeeze(p(3,:)));
    bot = min(min(fDL_M.data, fDL_R.data));
    top = max(max(fDL_M.data,fDL_R.data));
    %axis([ax(1) ax(2) bot top])
    text(0.2*top,0.2*ax(2),['\rho = ' num2str(a) '  lag = ' num2str(lags(b))])
    disp(['   dLOD \rho = ' num2str(a) '  lag = ' num2str(lags(b))])

    subplot(4,1,4)
    plot(lags,squeeze(c1),lags,squeeze(c2),lags,squeeze(c3))
    xlabel('lag (days)')
    ylabel('correlation')
    %axis([-10*dT 10*dT -1 1])

figure(2),clf
   Mplot = plot(fDL_M,'Color',col_M);
    hold on
    EIplot = plot(fDL_R,'Color',col_R);
    xlabel('Time')
    ylabel('\Delta LOD')
    [a,b] = max(squeeze(p(3,:)));
    bot = min(min(fDL_M.data, fDL_R.data));
    top = max(max(fDL_M.data,fDL_R.data));
    %axis([ax(1) ax(2) bot top])
    text(0.2*top,0.2*ax(2),['\rho = ' num2str(a) '  lag = ' num2str(lags(b))])
    disp(['   dLOD \rho = ' num2str(a) '  lag = ' num2str(lags(b))])

end
