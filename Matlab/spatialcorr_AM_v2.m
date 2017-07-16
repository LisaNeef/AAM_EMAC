% spatialcorr_AM_v2.m--------------
%
%  Load components of AAM computed from an EMAC run (horizontal fields)
%  and compute the covariance and correlation with the total, 
%  where the total is filtered to certain timescales.
%  v2 started 18 Aug 2010
%
%-----------------------------------------------------------------
clear all;



addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/snctools/');
addpath('/home/ig/neef/MFiles/netcdf/mexcdf/mexnc/');
addpath('/home/ig/neef/MFiles/utilities/');


%--set inputs:
Tmin = 150;	
Tmax = 230;
comp = 1;	% which component to look at - choose 1 for equatorial complex, 3 for axial
fil_choice = 2;	% timeseries filter: 1=idealfilter, 2=butterworth

% load the data and read the grid.
load aam_mass_t7_T42L39_1980s_long.mat
X = Xw+Xm;
[ncomp,nt,nlat,nlon] = size(mm);

%---convert MJD array to a t axis that MATLAB likes
t = MJD*0;
nt = size(t,2);
[y, m, d] = mjd2date(MJD);
for ii=1:nt, t(ii)=datenum([y(ii) m(ii) d(ii)]); end

%-- generate filtered timeseries of the total AAM at the desired timescale
fmax = 1/Tmin;
fmin = 1/Tmax;
Fw = X*0;
Fm = X*0;

if fil_choice == 1
  % Option 1: MATLAB canned filter 
  Xtm1 = timeseries(Xm(1,:),t,'isDatenum',true);
  Xtm2 = timeseries(Xm(2,:),t,'isDatenum',true);
  Xtm3 = timeseries(Xm(3,:),t,'isDatenum',true);
  Xfm1 = idealfilter(Xtm1,[fmin,fmax],'pass');
  Xfm2 = idealfilter(Xtm2,[fmin,fmax],'pass');
  Xfm3 = idealfilter(Xtm3,[fmin,fmax],'pass');

  Xtw1 = timeseries(Xw(1,:),t,'isDatenum',true);
  Xtw2 = timeseries(Xw(2,:),t,'isDatenum',true);
  Xtw3 = timeseries(Xw(3,:),t,'isDatenum',true);
  Xfw1 = idealfilter(Xtw1,[fmin,fmax],'pass');
  Xfw2 = idealfilter(Xtw2,[fmin,fmax],'pass');
  Xfw3 = idealfilter(Xtw3,[fmin,fmax],'pass');

  %--recombine the filtered AAM:
  Fm(1,:) = squeeze(Xfm1.data)+i*squeeze(Xfm2.data);
  Fm(3,:) = squeeze(Xfm3.data);
  Fw(1,:) = squeeze(Xfw1.data)+i*squeeze(Xfw2.data);
  Fw(3,:) = squeeze(Xfw3.data);
end

if fil_choice == 2
  % Option 2: digital butterworth filter, degree 3
  [B,A] = butter(3,2*[fmin fmax]);
  Fw(1,:) = filter(B,A,Xw(1,:))+i*filter(B,A,Xw(2,:)); 
  Fw(3,:) = filter(B,A,Xw(3,:))+i*filter(B,A,Xw(3,:));
  Fm(1,:) = filter(B,A,Xm(1,:))+i*filter(B,A,Xm(2,:)); 
  Fm(3,:) = filter(B,A,Xm(3,:))+i*filter(B,A,Xm(3,:));
end



%--initialize large arrays holding correlations (R) and covariances (S)
Rm = zeros(nlat,nlon);
Cm = zeros(nlat,nlon);
Rw = zeros(nlat,nlon);
Cw = zeros(nlat,nlon);

%--compute covariances and correlations between components and the total

for ilat = 1:nlat
  for ilon = 1:nlon
    if comp == 1 
      yw = squeeze(ww(1,:,ilat,ilon))+i*squeeze(ww(1,:,ilat,ilon));
      ym = squeeze(mm(1,:,ilat,ilon))+i*squeeze(mm(1,:,ilat,ilon));
    else
      yw = squeeze(ww(3,:,ilat,ilon));
      ym = squeeze(mm(3,:,ilat,ilon));
    end
        [sigma_m,lags] = xcov(ym,Fm(comp,:),'coeff');		% covariance
        [rho_m,lags] = xcorr(ym,Fm(comp,:),'coeff');		% correlation
        [sigma_w,lags] = xcov(yw,Fw(comp,:),'coeff');		% covariance
        [rho_w,lags] = xcorr(yw,Fw(comp,:),'coeff');		% correlation
        %Rw(ilat,ilon) = rho_w(find(lags == 0));
        %Cw(ilat,ilon) = sigma_w(find(lags == 0));
        %Rm(ilat,ilon) = rho_m(find(lags == 0));
        %Cm(ilat,ilon) = sigma_m(find(lags == 0));
        Rw(ilat,ilon) = max(rho_w);
        Cw(ilat,ilon) = max(sigma_w);
        Rm(ilat,ilon) = max(rho_m);
        Cm(ilat,ilon) = max(sigma_m);
  end
end

%-- transform longitude from 0:360 to -180:180
dum = find(lon > 180)
lon2=lon;
lon2(dum) = -(360 - lon(dum));
[a,b] = sort(lon2);
lon = a;
Rm = Rm(:,b);
Cm = Cm(:,b);
Rw = Rw(:,b);
Cm = Cw(:,b);

%--plot maps of correlations and covariances

  figure(1),clf
  axesm miller; framem;gridm
  c = load('coast');
  caxis([0 1])
  plotm(c.lat,c.long,'Color',1*ones(3,1))
  contourfm(lat,lon,squeeze(abs(Rm)));
  title('mass term correlation')
  colorbar

  figure(2),clf
  axesm miller; framem;gridm
  plotm(c.lat,c.long,'Color',1*ones(3,1))
  caxis([0 1])
  contourfm(lat,lon,squeeze(abs(Rw)));
  title('wind term correlation')
  colorbar


  figure(3),clf
  axesm miller; framem;gridm
  c = load('coast');
  plotm(c.lat,c.long,'Color',1*ones(3,1))
  contourfm(lat,lon,squeeze(abs(Cm)));
  title('mass term covariance')
  colorbar

  figure(4),clf
  axesm miller; framem;gridm
  plotm(c.lat,c.long,'Color',1*ones(3,1))
  contourfm(lat,lon,squeeze(abs(Cw)));
  title('wind term covariance')
  colorbar

