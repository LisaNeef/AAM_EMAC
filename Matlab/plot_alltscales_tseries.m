% ---plot_alltscales_tseries.m----------------
%
%  make a plot comparing the AAM excitation functions in time, for all timescales 
%  the purpose is to examine differences / similarities between EMAC, the ERA reanalyses, and observations.
%  Make 3 plots: 1 for each aam component.
%  In each plot, compare: EMAC, ERA, and GEO-OAM
%  Started 20 April 2010

%---load the filtered data

Tmin = 0;
Tmax = 20000;
eopset = 'IERS';
RAset = 'ERA';
[t,M,O,R,G,H] = aam_compare_tscales(Tmin,Tmax,eopset,RAset);


%-- manipulate & combine timeseries
N = G-O;
nt = length(t);
D = zeros(3,3,nt);
D(1,:,:) = M;
D(2,:,:) = R;
D(3,:,:) = N;

D(:,3,:) = D(:,3,:)*1000;	% conversion from s to ms

names = {'EMAC','ERA','OBS-OAM'};

%---some plot settings
Mcol = [0.7513    0.2551    0.5060];
Rcol = [0.0714    0.5216    0.0967];
Gcol = zeros(1,3);
Ncol = 0.7*Gcol;
c = rand(4,3);
  c(1,:) = Mcol;
  c(2,:) = Rcol;
  c(3,:) = Gcol;
LH = zeros(3,1);
LW = 2;
TITLE = {'\chi_1', '\chi_2', '\Delta LOD'};
DT = 0.3*(max(t)-min(t));	% the stretch of time we look at.
%---plots!

for ii = 1:3
  figure(ii),clf
  hold on
  for is = 1:3
    LH(is) = plot(t,squeeze(D(is,ii,:)),'Color',c(is,:),'LineWidth',LW);
  end 
  title(TITLE(ii))
  legend(LH,names,'Location','best','Orientation','Horizontal')
  datetick('x','YYYY')
  if ii < 3
     a = [max(t)-DT max(t) -150 150];
     axis(a)
     ylabel('mas')
  else
     a = [max(t)-DT max(t)  -2 2];
     ylabel('ms')
  end
end
 
