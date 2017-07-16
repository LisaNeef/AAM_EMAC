%function plot_compare_aam_integrals(comp)
% plot routine to compute the AAM computed via two integration codes:
%	(1) aam_emac_massintegral.m
%	(2) aam_emac_volintegral.m
% ...and also compare to the same quantities from ERA
% 22 july 2010
% MODS
%  23 Mar 2011: several updates almost a year later..l.


% options:
useERA = 1;

% load aam timeseries for the same decade.
rundir = '/dsk/nathan/lisa/EMAC_ex/t7_T42L39/mat/'
fm = char([rundir,'aam_mass_t7_T42L39_1960s.mat']);
fv = char([rundir,'aam_vol_t7_T42L39_1960s.mat']);

load (fm)
  ntM = size(Xw,2);
  M = zeros(3,ntM);
  M(1,:) = squeeze(Xw(comp,:));
  M(2,:) = squeeze(Xm(comp,:));
  M(3,:) = squeeze(Xw(comp,:))+squeeze(Xm(comp,:));
  Mmjd = MJD;
load (fv)
  ntV = size(Xw,2);
  V = zeros(3,ntV);
  V(1,:) = squeeze(Xw(comp,:));
  V(2,:) = squeeze(Xm(comp,:));
  V(3,:) = squeeze(Xw(comp,:))+squeeze(Xm(comp,:));
  Vmjd = MJD;

MDT = M*0+NaN;
VDT = V*0+NaN;
rad2mas = (180/pi)*60*60*1000;
LOD0 = 24*60*60*1000;
gV = find(isfinite(V(iterm,:))==1);
gM = find(isfinite(V(iterm,:))==1);
if comp == 3
  for iterm = 1:3
    VDT(iterm,:) = LOD0*detrend(V(iterm,gV),'constant');
    MDT(iterm,:) = LOD0*detrend(M(iterm,gM),'constant');
  end
else
  for iterm = 1:3
    VDT(iterm,:) = rad2mas*detrend(V(iterm,gV),'constant');
    MDT(iterm,:) = rad2mas*detrend(M(iterm,gM),'constant');
  end
end


% make time axes for EMAC data
tM = Mmjd*0;
tV = Vmjd*0;
[yM, mM, dM] = mjd2date(Mmjd);
[yV, mV, dV] = mjd2date(Vmjd);
for ii=1:ntM
  if isfinite(Mmjd(ii)) == 1, tM(ii)=datenum([yM(ii) mM(ii) dM(ii)]); else tM(ii) = NaN; end
end
for ii=1:ntV
  if isfinite(Vmjd(ii)) == 1, tV(ii)=datenum([yV(ii) mV(ii) dV(ii)]); else tV(ii) = NaN; end
end


% load ERA data?
if useERA == 1
  [Xw_EI,Xm_EI,MJD_EI] = read_EFs('aam','ERAinterim',1);
  [Xw_E4,Xm_E4,MJD_E4] = read_EFs('aam','ERA40',1);
  s = find(MJD_E4 == min(MJD_EI));
  Xtemp_E4 = Xw_E4(:,1:s-1)+Xm_E4(:,1:s-1);
  Xtemp_EI = Xw_EI+Xm_EI;
  XR = [Xtemp_EI,Xtemp_E4];
  MJDR = [MJD_E4(1:s-1),MJD_EI];
  g_R = find(isfinite(XR(3,:)) == 1);
  XR(1:2,:) = XR(1:2,:)*rad2mas;                        % (conversion to mas)
  XR(3,:) = LOD0*(XR(3,:) - mean(XR(3,g_R)));           % (conversion to LOD)

  tR = MJDR*0;
  ntR = size(tR,2);
  [yR, mR, dR] = mjd2date(MJDR);
  for ii=1:ntR, tR(ii)=datenum([yR(ii) mR(ii) dR(ii)]); end
end

% let's pick some colors
Mcol = 0.7*ones(1,3);
Vcol = rand(3,1);
Rcol = 0*ones(3,1);

for iterm = 1:3
  figure(iterm),clf
    xM = squeeze(MDT(iterm,:));
    xV = squeeze(VDT(iterm,:));
    Mh = plot(tM,xM,'-','Color',Mcol);
    hold on
    Vh = plot(tV,xV,'-','Color',Vcol);
    datetick('x','YYYY')
    switch iterm
      case 1
        title(['\chi_',num2str(comp),'Motion Term'])
      case 2
        title(['\chi_',num2str(comp),'Mass Term'])
      case 3
        title(['\chi_',num2str(comp)])
    end
    if comp == 3, ylabel('ms'); else ylabel('mas'), end
    lh = [Mh(1) Vh(1)];
    legend(lh,'Mass Integral','Volume Integral')
end


% plot export settings


plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';


rend = 'opengl';
LW = 4;
ph = 6;        % paper height
pw = 20;        % paper width
fs = 23;        % fontsize


for ii = 1:3 
  if ii == 1, fig_name = [plot_dir,'comp_massint_volint_X',num2str(ii),'.png']; end
  if ii == 2, fig_name = [plot_dir,'comp_massint_volint_X',num2str(ii),'.png']; end
  if ii == 3, fig_name = [plot_dir,'comp_massint_volint_X',num2str(ii),'.png']; end
  exportfig(ii,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');
end







