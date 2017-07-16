% Compare the AAM produced by EMAC runs at various resolutions
% 19 Oct 2010
%--------------------------------------------------------------------------


clear all;

%---specify filtering periods

%T = [400,700];		% interannual variations
T = [30,90];		% annual cycle
fil_order = 3;

% choose type of filter to plot
filtype = 1;		% butterworth filter
%filtype = 2;		% chebyeshev-1 filter
%filtype = 3;		% chebyeshev-2 filter


%---specify the mat files for the runs to compare

d1 = '/dsk/nathan/lisa/CCMval/Mat/';
n1 = 'CCMVal';

d2 = '/dsk/nathan/lisa/ex/t7_T42L39/mat/';
n2 = 'T42L39';

d3 = '/dsk/nathan/lisa/ex/ref2_T31L39/mat/';
n3 = 'T31L39';

nf = 3;

names = {n1,n2, n3};
dd = 1960:10:1970;
years = 1960:1:1980;
if T(2) < 100
  dd = 1960;
  years = 1960:1970
end
lys = leapyear(years);
n_normal = length(find(lys == 0));
n_leap = length(find(lys == 1));
ndays = n_normal*365+n_leap*366;

fp1 = [d1,'aam_vol_ccmval_'];
fp2 = [d2,'aam_vol_t7_T42L39_'];
fp3 = [d3,'aam_vol_ref2_T31L39_'];
fp = {fp1, fp2, fp3};
  
% load constants
aam_constants_gross
Re = Re_m;      % (use earth radius in meters)
pm12 = -1.608*0.684*Re^4/(g*(C-A));
pw12 = -1.608*Re^3/(Q*g*(C-A));

%---PLOT SETTINGS HERE--------------------
col = rand(nf,3);
col = rand(nf,3);
  col(1,:) = 0.5*ones(1,3);
  col(2,:) = [0.8797    0.5944    0.3127];
  col(3,:) = [0.6753    0.3868    0.4624];
lh = zeros(nf,1);
LW = 2;
YL = {'\chi_1 (mas)', '\chi_2 (mas)', '\Delta LOD (ms)'};
TT = {'Wind Term','Mass Term'};

ax = zeros(3,2);
ax(1,:) = 20*[-1,1];
ax(2,:) = 40*[-1,1];
ax(3,:) = 0.5*[-1,1];



%---------------------------------------

% cycle through available decades and put together long timeseries

figure(1),clf
figure(2),clf
figure(3),clf

for ifile = 1:nf
  xm = zeros(3,ndays)+NaN;
  xw = zeros(3,ndays)+NaN;
  mjd = zeros(1,ndays)+NaN;
  k0 = 1;
  X = zeros(2,3,ndays);
  XF = zeros(2,3,ndays)+NaN;

  for idec = 1:length(dd)
    dec = dd(idec)
    f = [char(fp(ifile)),num2str(dec),'s.mat'];
    load(f)
  
    kf = k0+length(MJD)-1;
    xm(:,k0:kf) = Xm;
    xw(:,k0:kf) = Xw;
    mjd(k0:kf) = MJD;
    k0 = k0+kf;      
  end

  % now take out LT mean and convert units, ignoring NaNs.  
  gw = find(isfinite(xw(1,:)) == 1);
  gm = find(isfinite(xm(1,:)) == 1);
  X(1,1,gw) = rad2mas*(detrend(xw(1,gw)));
  X(1,2,gw) = rad2mas*(detrend(xw(2,gw)));
  X(1,3,gw) = LOD0_ms*(detrend(xw(3,gw)));
  X(2,1,gm) = rad2mas*(detrend(xm(1,gm)));
  X(2,2,gm) = rad2mas*(detrend(xm(2,gm)));
  X(2,3,gm) = LOD0_ms*(detrend(xm(3,gm)));

  % clear some space sothat MATLAB moves better
  clear Xm Xw MJD lat lon ww

  % craft a time array that MATLAB can deal with
  nt = size(mjd,2);
  t = mjd*0+NaN;
  [y, m, d] = mjd2date(mjd);
  for ii=1:nt
    if isfinite(mjd(ii)) ==1, t(ii)=datenum([y(ii) m(ii) d(ii)]); end
    % weird correction for CCMVal data 
    if mjd(ii) < mjd(1), t(ii) = NaN; end
  end

  % now cycle through all 6 AAM components and apply the time filter of choice.
  for ii = 1:3
    figure(ii)
    for iterm = 1:2
      gg = find(isfinite(X(iterm,ii,:)) == 1);
      Xgood = squeeze(X(iterm,ii,gg));
      xf = cfilter(Xgood,1,T(1),T(2),fil_order,'days');
      XF(iterm,ii,gg) = xf(:,filtype); 

  % throw the filtered series onto the appropriate plot
      subplot(1,2,iterm) 
      plot(t,squeeze(XF(iterm,ii,:)),'LineWidth',LW,'Color',col(ifile,:));
      hold on
      axis([min(t), max(t), ax(ii,1), ax(ii,2)]);
      ylabel(YL(ii))
      xlabel('time')
      title(TT(iterm))
      text(t(1),(1-0.1*ifile)*ax(ii,2),names(ifile),'Color',col(ifile,:))
      datetick('x','yyyy')
      if T(2) < 100; axis([min(t), min(t)+0.25*(max(t)-min(t)), ax(ii,1), ax(ii,2)]); end

    end   % loop over mass and wind terms
  end     % loop over components 1, 2, 3

end  % loop over datasets


