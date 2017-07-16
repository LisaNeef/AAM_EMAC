%---plot_tscales_tseries.m
%
%  make plots comparing the AAM excitation functions in time, for the different timescales 
%  the purpose is to examine differences / similarities between EMAC, the ERA reanalyses, and observations.
%  Make 5 plots: 1 for each timescale
%  In each plot, compare: EMAC, ERA, and GEO-OAM
%  Started 19 April 2010
%  Mods
%	27 may 2010 - change plot formatting and add start/stop years as inputs.
%  INPUT:
%    comp: component (1, 2, or 3)
%    iT: which timescale to look at: 
%	0: no filtering.
%	1: subseasonal
%	2: terannual
%	3: semiannual
%	4: annual
%	5: interannual
%    dataset: vector containing the datasets to look at:
%       1: EMAC
%       2: ERA
%	3: obs raw
%	4: obs-oam
%    y0, yf: start and stop years.
%function plot_tscales_tseries(comp,iT,dataset,y0,yf)

% temp inputs:
comp = 3;
iT = 3;
dataset = [1,3];
y0=1975;
yf=2000;


%---establish timescales

TT = zeros(5,2);
TT(1,:) = [30,90];      % subseasonal.
TT(2,:) = [90,150];     % terannual
TT(3,:) = [150,230];    % semiannual
TT(4,:) = [230,450];    % annual
%TT(5,:) = [730,1825];   % interannual
TT(5,:) = [730,2555];   % interannual
nT = size(TT,1);

if iT == 0, TT = [0,0]; iT = 1; end

%---some plot settings

Mcol = [0.7513    0.2551    0.5060];
Rcol = [0.6    0    0.8];
Gcol = zeros(1,3);
Ncol = 0.7*ones(1,3);
c = rand(4,3);
  c(1,:) = Mcol;
  c(2,:) = Rcol;
  c(3,:) = Gcol;
  c(4,:) = Ncol;
LH = zeros(length(dataset),1);
LW = 3;
TITLE1 = {'Subseasonal Variation', 'Terannual Variation', 'Semiannual Variation', 'Annual Variation', 'Interannual Variation'};
TITLE2 = {'\chi_1 (mas)', '\chi_2 (mas)', '\Delta LOD (ms)'};



%---loop over timescales and generate plots


  % load the filtered timeseries and combine equatorial terms
  % (note that LOD changes are converted to ms)
  [t,M,O,R,G,H,MJD] = aam_compare_tscales(TT(iT,1),TT(iT,2),'IERS','ERA');

  mjd0 = date2mjd(y0,1,1,0,0,0);
  if min(MJD) > mjd0, mjd0 = min(MJD), end
  mjdf = date2mjd(yf,1,1,0,0,0);
  if max(MJD) < mjdf, mjdf = max(MJD), end
  t0 = find(MJD == mjd0);
  tf = find(MJD == mjdf);

  %-- manipulate & combine timeseries
  N = G-O;
  nt = length(t);
  D = zeros(4,nT,nt);
  D(1,iT,:) = M(comp,:);
  D(2,iT,:) = R(comp,:);
  D(3,iT,:) = G(comp,:);
  D(4,iT,:) = N(comp,:);


  if comp == 3, D = D*1000; end	% conversion from s to ms
  names = {'EMAC','ERA','OBS','OBS-OAM'};
  
  figure(iT),clf
  hold on
  for is = 1:length(dataset)
    ss = dataset(is);
    LH(is) = plot(t(t0:tf),squeeze(D(ss,iT,t0:tf)),'Color',c(ss,:),'LineWidth',LW);
  end 
  datetick('x','YYYY')
  title(TITLE1(iT))
  legend(LH,names(dataset),'Location','best','Orientation','Horizontal')
  ylabel(TITLE2(comp))
