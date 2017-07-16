%  quickie program to compare the diff between Xm when we use mslp and when we use aps
% 22 Feb 2010

addpath('/home/ig/neef/MFiles/netcdf/');


load 'testrun_1983_1985_mslp.mat'
Xm_mslp = Xm;
MJD_mslp = MJD;

load 'testrun_1983_1985_ps.mat'
Xm_ps = Xm;
MJD_ps = MJD;



figure(1),clf
  psplot = plot(MJD_ps,Xm_ps(1,:),'LineWidth',3,'Color',rand(3,1));
  hold on
  mslpplot = plot(MJD_mslp,Xm_mslp(1,:),'LineWidth',3,'Color',rand(3,1));
  lhandle = [psplot(1) mslpplot(1)];
  legend(lhandle,'using ps','using mslp')
  ylabel('X1')
  title('AEF-1 Mass Terms')

figure(2),clf
  psplot = plot(MJD_ps,Xm_ps(2,:),'LineWidth',3,'Color',rand(3,1));
  hold on
  mslpplot = plot(MJD_mslp,Xm_mslp(2,:),'LineWidth',3,'Color',rand(3,1));
  lhandle = [psplot(1) mslpplot(1)];
  legend(lhandle,'using ps','using mslp')
  ylabel('X2')
  title('AEF-2 Mass Terms')

%---let's also compare these two fields

dir = '/dsk/nathan/lisa/Data/CCMval';
ff_mslp = 'scout02____dm_mslp_199812.nc';
ff_ps = 'O3_199812_aps.nc';


[data_mslp name] = r_netcdf(dir,ff_mslp);
[data_ps name] = r_netcdf(dir,ff_ps);

mslp = data_mslp{4};
ps = data_ps{4};

figure(3),contourf(mslp(:,:,1)')
title('MSLP')
figure(4),contourf(ps(:,:,1)')
title('ps')




