% make quick example plots of u and ps fields frome emac.


clear all;

load emac_uvpfields

pfield = squeeze(Ps(:,:,[500]));
ufield = squeeze(U300(:,:,[1]));


pmean = mean(Ps(:,:,1:365),3);
umean = mean(U300(:,:,1:365),3);

pvar= var(Ps(:,:,1:365),0,3);
uvar= var(U300(:,:,1:365),0,3);


% also straighten out longitude and sort.
dum = find(lon > 180);
lon1 = lon;
lon2 = lon;
lon2(dum) = lon(dum)-360;
[a,b] = sort(lon2);
lon = lon2(b);

XT = [-180 -120 -60 0 60 120 180];
YT = [-90 -60 -30 0 30 60 90];



figure(4),clf
    m_proj('miller');
    m_pcolor(lon,lat,pmean')
    shading interp
    m_coast('line','Color',zeros(3,1),'LineWidth',2);
    m_grid('xtick',XT,'tickdir','out','ytick',YT,'yaxislocation','right','box','fancy','fontsize',7);
    colorbar('Location','SouthOutside')
    title('EMAC Mean Surface Pressure, 1960')

figure(5),clf
    m_proj('miller');
    m_pcolor(lon,lat,umean')
    shading interp
    m_coast('line','Color',zeros(3,1),'LineWidth',2);
    m_grid('xtick',XT,'tickdir','out','ytick',YT,'yaxislocation','right','box','fancy','fontsize',7);
    colorbar('Location','SouthOutside')
    title('EMAC Mean Zonal Wind, 1960')
