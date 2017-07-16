% plot_spatialAM.m----------------
%
%  Contour plots of the covariance between spatial components of AM and the global value
%  22 April 2010




%---paths and file settings
addpath('/home/ig/neef/MFiles/AM/');
addpath('/home/ig/neef/MFiles/netcdf/');
addpath('/home/ig/neef/MFiles/utilities/');
addpath('/home/ig/neef/MFiles/utilities/m_map/');

%---compute AAM prefactors (see Barnes et al., 1983, p. 48)

d2r = pi/180.;         %  degrees to radian
LOD0 = 24*60*60;        % length of day (in s) as implied by earth rot. rate
r2m = (180/pi)*60*60*1000;    % conversion of radians to mas

%--load data: chi terms in space and time.
load spatialAM_1995_2000.mat
[nv,nlon,nlat,nt] = size(xxw);


% combine equatorial EFs to two complex functions
% and convert to equatorial terms to mas while converting axial term to LOD changes
% while also removing the mean in each case
dum = squeeze(xxm(3,1,1,:));
g = find(isfinite(dum) == 1);
meanm = mean(xxm(:,:,:,g),4);
meanw = mean(xxw(:,:,:,g),4);
XXw = zeros(2,nlon,nlat,nt);
XXm = zeros(2,nlon,nlat,nt);
for it = 1:nt
  XXm(2,:,:,it) = LOD0*(xxm(3,:,:,it) - meanm(3,:,:));
  XXw(2,:,:,it) = LOD0*(xxw(3,:,:,it) - meanw(3,:,:));
  XXm(1,:,:,it) = r2m*squeeze(xxm(1,:,:,it)-meanm(1,:,:)+i*xxm(2,:,:,it)-i*meanm(2,:,:));
  XXw(1,:,:,it) = r2m*squeeze(xxw(1,:,:,it)-meanw(1,:,:)+i*xxw(2,:,:,it)-i*meanw(2,:,:));
end


%---establish timescales

TT = zeros(5,2);
TT(1,:) = [30,90];      % subseasonal.
TT(2,:) = [90,150];     % terannual
TT(3,:) = [150,230];    % semiannual
TT(4,:) = [230,450];    % annual
TT(5,:) = [730,1825];   % interannual
nT = size(TT,1);


[t,M,O,R,G,H,MJD2] = aam_compare_tscales(TT(iT,1),TT(iT,2),'IERS','ERA');

%---slap everything on the same timeaxis.  Call these new functions Y
k0 = find(MJD2 == min(MJD)+0.5);
kf = length(MJD2);
nt2 = kf-k0+1;
ym = XXm(:,:,:,1:nt2);
yw = XXw(:,:,:,1:nt2);
Y = zeros(2,nt2);

for kk=k0:kf
  Y(1,kk-k0+1) = M(1,kk)+i*M(2,kk);
  Y(2,kk-k0+1) = M(3,kk);
end



%---cycle through points and compute covariances.

rho_w = zeros(2,nlon,nlat);
lag_w = zeros(2,nlon,nlat);
rho_m = zeros(2,nlon,nlat);
lag_m = zeros(2,nlon,nlat);

for ilat = 1:nlat
  for ilon = 1:nlon
    ymt = squeeze(ym(:,ilon,ilat,:));
    ywt = squeeze(yw(:,ilon,ilat,:));
    g = find(isfinite(ymt) == 1);
    for jj = 1:2
      [cm,lagsm] = xcorr(ymt(jj,:),Y(jj,:),'coeff');
      [cw,lagsw] = xcorr(ywt(jj,:),Y(jj,:),'coeff');
      [am,bm] = max(cm);
      [aw,bw] = max(cw);
      rho_w(jj,ilon,ilat) = aw;
      lag_w(jj,ilon,ilat) = lagsw(bw);
      rho_m(jj,ilon,ilat) = am;
      lag_m(jj,ilon,ilat) = lagsm(bm);
    end
  end
end


%----plots
TITLE1 = {'Subseasonal', 'Terannual', 'Semiannual', 'Annual', 'Interannual'};
TITLE2 = {'\chi_{eq}', '\Delta LOD'};

disp(TITLE1(iT))

dum = find(lon > 180);
lon1 = lon;
lon2 = lon;
lon2(dum) = lon(dum)-360;
[a,b] = sort(lon2);

for ii=1:2

  figure(ii),clf
  subplot(1,2,1)
  m_proj('miller');
  m_pcolor(lon2(b),lat,squeeze(real(rho_w(ii,b,:)))')
  shading interp
  m_coast('line','Color',zeros(3,1),'LineWidth',2);
  if ii == 1 
    title([TITLE1(iT) '\chi_{eq} Motion Term'])
  else
    title([TITLE1(iT) '\Delta LOD Motion Term'])
  end
  colorbar

  subplot(1,2,2)
  m_proj('miller');
  m_pcolor(lon2(b),lat,squeeze(real(rho_m(ii,b,:)))')
  shading interp
  m_coast('line','Color',zeros(3,1),'LineWidth',2);
  colorbar
  if ii == 1 
    title([TITLE1(iT) '\chi_{eq} Motion Term'])
  else
    title([TITLE1(iT) '\Delta LOD Motion Term'])
  end

end







