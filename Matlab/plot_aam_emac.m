% a quickie plot routine to check the AAM computed from emac runs  
% and compare to ERA.
% Need to have the output of aam_emac or aam_emac_massintegral loaded.
% 20 july 2010

% options:
useERA = 1;

% convert X3 to LOD
g = find(isfinite(Xw(3,:)) == 1);
LOD0 = 24*60*60*1000;
LODw = LOD0*(Xw(3,:)-mean(Xw(3,g)));
LODm = LOD0*(Xm(3,:)-mean(Xm(3,g)));

% convert X1 and X2 to mas
rad2mas = (180/pi)*60*60*1000;
Xeq_w = Xw(1:2,:)*rad2mas;
Xeq_m = Xm(1:2,:)*rad2mas;

% make time axis for EMAC data
t = MJD*0;
nt = size(t,2);
[y, m, d, h] = mjd2date(MJD);
for ii=1:nt
  if isfinite(MJD(ii)) == 1
    t(ii)=datenum([y(ii) m(ii) d(ii)]);
  else
    t(ii) = NaN;
  end
end

% load ERA data?
if useERA == 1
  [Xw_EI,Xm_EI,MJD_EI] = read_EFs('aam','ERAinterim');
  [Xw_E4,Xm_E4,MJD_E4] = read_EFs('aam','ERA40');
  s = find(MJD_E4 == min(MJD_EI));
  Xtemp_E4 = Xw_E4(:,1:s-1)+Xm_E4(:,1:s-1);
  Xtemp_EI = Xw_EI+Xm_EI;
  Xtemp_R = [Xtemp_EI,Xtemp_E4];
  MJDR = [MJD_E4(1:s-1),MJD_EI];
  g_R = find(isfinite(Xtemp_R(3,:)) == 1);
  Xeq_R = Xtemp_R(1:2,:)*rad2mas;                        % (conversion to mas)
  LOD_R = LOD0*(Xtemp_R(3,:) - mean(Xtemp_R(3,g_R)));      % (conversion to LOD)
end

% also load the obs.
% EOP-implied AAM terms: X1 and X2 in mas, dlodg in ms 
[X1g,X2g,dlodg,MJDG] = read_eops('IERS');
LOD_G = (dlodg-mean(dlodg))*1000;

% make time axes for comparison data
tR = MJDR*0;
tG = MJDG*0;
ntR = size(tR,2);
ntG = size(tG,1);
[yR, mR, dR] = mjd2date(MJDR);
[yG, mG, dG] = mjd2date(MJDG);
for ii=1:ntR, tR(ii)=datenum([yR(ii) mR(ii) dR(ii)]); end
for ii=1:ntG, tG(ii)=datenum([yG(ii) mG(ii) dG(ii)]); end


figure(1),clf
  subplot(2,1,1)
    plot(t,LODw,'.-','Color',rand(3,1))
    datetick('x','YYYY')
    title('\Delta LOD In EMAC Run t8-T42L39: Motion Term')
    ylabel('ms')
  subplot(2,1,2)
    plot(t,LODm,'.-','Color',rand(3,1))
    datetick('x','YYYY')
    ylabel('ms')
    title('\Delta LOD In EMAC Run t8-T42L39: Mass Term')

figure(2),clf
  subplot(2,1,1)
    plot(t,Xeq_w(1,:),'.-','Color',rand(3,1))
    datetick('x','YYYY')
    title('\chi_1 In EMAC Run t8-T42L39: Motion Term')
    ylabel('mas')
  subplot(2,1,2)
    plot(t,Xeq_m(1,:),'.-','Color',rand(3,1))
    datetick('x','YYYY')
    ylabel('mas')
    title('\chi_1 In EMAC Run t8-T42L39: Mass Term')

figure(3),clf
  subplot(2,1,1)
    M = plot(t,LODw+LODm,'.-','Color',rand(3,1));
    hold on
   % G = plot(tG,LOD_G,'Color',0.7*ones(3,1),'LineWidth',3);
    R = plot(tR,LOD_R,'Color',0*ones(3,1));
    datetick('x','YYYY')
    axis([min(t) max(t) 1.3*min(LODw+LODm) 1.3*max(LODw+LODm)])
    lhandle = [M(1) R(1)]; 
    legend(lhandle,'EMAC','ERA')
    title('Comparing LOD Variations')
  subplot(2,1,2)
    M = plot(t,Xeq_w(1,:)+Xeq_m(1,:),'.-','Color',rand(3,1));
    hold on
    R = plot(tR,Xeq_R(1,:),'Color',0*ones(3,1));
    G = plot(tG,X1g,'Color',0.7*ones(3,1),'LineWidth',3);
    datetick('x','YYYY')
    lhandle = [M(1) R(1) G(1)];
    legend(lhandle,'EMAC','ERA','OBS')
    title('Comparing \chi_1 Variations')
    axis([min(t) max(t) 1.3*min(Xeq_w(1,:)+Xeq_m(1,:)) 1.3*max(Xeq_w(1,:)+Xeq_m(1,:))])

