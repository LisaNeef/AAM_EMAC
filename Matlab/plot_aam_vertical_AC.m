% plot_aam_vertical_AC
%
% Make a pressure-time contour plot showing the vertical contributions to 
% the wind term of AAM
% 20 Dec 2010
%


%---------------------------------------------------

% temp inputs:
comp = 3;
lev = 1:27;

% load input data

datadir = '/dsk/nathan/lisa/ex/CCMval/mat/';
ff = 'aam_vol_vert_sliceCCMval_1980_long.mat';

load([datadir,ff])

X = squeeze(Xw(comp,:));
x = squeeze(ww(comp,:,:));
nlev = length(lev);

% rearrange the years
[y,m,d] = mjd2date(MJD);
ymax = max(y);
ymin = min(y);
ny = ymax-ymin+1;
xvert = zeros(nlev,365,ny);
Y = ymin:1:ymax;

for iy = 1:ny
  yt = find(y == Y(iy));
  xvert(:,:,iy) = x(yt(1:365),:)';
end
