% very simple plot comparing the geodetic and modeled polar motion excitation functions
% 15 April 2010

[t,M,R,G,N] = aam_compare_tscales(0,4000);

mcol = 0.7*ones(3,1);
gcol = rand(3,1);


figure(1),clf
  mplot = plot(t,M(1,:),'Color',mcol);
  hold on
  gplot = plot(t,G(1,:),'Color',gcol,'LineWidth',2);
  datetick('x','YYYY')
  ylabel('(mas)')
  title('Polar Motion Excitation X1')
  lhandle = [mplot(1) gplot(1)];
  legend(lhandle,'EMAC','OBS')

% also a plot of the raw polar motion obs
dir = '/dsk/nathan/lisa/Data/IERS-ERP/';
ff = 'C04_1962_2000_notides.txt';
fname=[dir,ff];
eop = importdata(fname,' ',2);
mjd     = eop.data(:,1);
x       = eop.data(:,2);        % PM-x (mas)
[y, m, d] = mjd2date(mjd);
t2 = mjd*0;
nt = length(t2);
for ii=1:nt, t2(ii)=datenum([y(ii) m(ii) d(ii)]); end

figure(2),clf
  plot(t2,x,'LineWidth',2)
  datetick('x','YYYY')
  ylabel('(mas)')
  title('Polar Motion x')

