% plot routine to compare the AAM computed for and emac run and the CCMval data
% since the CCMval data's AAM was computed with a volume integral, the emac one should also be.
% 12 Aug 2010
% see notes v3 p47


clear all;

% ---some options:
useERA = 1;	% set to 1 to compare both runs to ERA data

% load aam timeseries for the same decade.

load aam_mass_t7_T42L39_1960s.mat
  Aw = Xw;
  Am = Xm;
  Amjd = MJD;
load aam_mass_ccmval_1960s.mat
  Bw = Xw;
  Bm = Xm;
  Bmjd = MJD;

% convert X3's to LOD
LOD0 = 24*60*60*1000;
gm = find(isfinite(Aw(3,:)) == 1);
gv = find(isfinite(Bw(3,:)) == 1);
Aw(3,:) = LOD0*(Aw(3,:)-mean(Aw(3,gm)));
Am(3,:) = LOD0*(Am(3,:)-mean(Am(3,gm)));
Bw(3,:) = LOD0*(Bw(3,:)-mean(Bw(3,gv)));
Bm(3,:) = LOD0*(Bm(3,:)-mean(Bm(3,gv)));

% convert X1 and X2 to mas
rad2mas = (180/pi)*60*60*1000;
Aw(1:2,:) = Aw(1:2,:)*rad2mas;
Am(1:2,:) = Am(1:2,:)*rad2mas;
Bw(1:2,:) = Bw(1:2,:)*rad2mas;
Bm(1:2,:) = Bm(1:2,:)*rad2mas;

% make time axes for EMAC data
tA = Amjd*0;
tB = Bmjd*0;
ntA = size(tA,2);
ntB = size(tB,2);
[yA, mA, dA] = mjd2date(Amjd);
[yB, mB, dB] = mjd2date(Bmjd);
for ii=1:ntA
  if isfinite(Amjd(ii)) == 1, tA(ii)=datenum([yA(ii) mA(ii) dA(ii)]); else tA(ii) = NaN; end
end
for ii=1:ntB
  if isfinite(Bmjd(ii)) == 1, tB(ii)=datenum([yB(ii) mB(ii) dB(ii)]); else tB(ii) = NaN; end
end


% load ERA data?
if useERA == 1
  [Xw_EI,Xm_EI,MJD_EI] = read_EFs('aam','ERAinterim');
  [Xw_E4,Xm_E4,MJD_E4] = read_EFs('aam','ERA40');
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
Acol = rand(3,1);
Bcol = rand(3,1);
Rcol = 0*ones(3,1);

for ii = 1:3
  figure(ii),clf
  subplot(3,1,1)
    A = plot(tA,Aw(ii,:),'-','Color',Acol);
    hold on
    B = plot(tB,Bw(ii,:),'-','Color',Bcol);
    datetick('x','YYYY')
    title('Motion Term')
    if ii == 3, ylabel('ms'); else ylabel('mas'), end
    lh = [A(1) B(1)];
    legend(lh,'Run t7','CCMVal')
  subplot(3,1,2)
    A = plot(tA,Am(ii,:),'-','Color',Acol);
    hold on
    B = plot(tB,Bm(ii,:),'-','Color',Bcol);
    datetick('x','YYYY')
    title('Mass Term')
    if ii == 3, ylabel('ms'); else ylabel('mas'), end
    lh = [A(1) B(1)];
    legend(lh,'Run t7','CCMVal')
  subplot(3,1,3)
    R = plot(tR,XR(ii,:),'-','Color',Rcol);
    hold on
    A = plot(tA,Am(ii,:)+Aw(ii,:),'-','Color',Acol);
    B = plot(tB,Bm(ii,:)+Bw(ii,:),'-','Color',Bcol);
    datetick('x','YYYY')
    title('Total')
    if ii == 3, ylabel('ms'); else ylabel('mas'), end
    lh = [A(1) B(1) R(1)];
    axis([min(tA) max(tA) 1.5*min(XR(ii,:)) 1.5*max(XR(ii,:))])
    legend(lh,'Run t7','CCMVal','ERA')
end
