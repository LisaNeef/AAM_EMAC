%  quickie program to compare the diff between X_model 
%  with and without the inverse barometer correction to surface pressure 
%  over the oceans.  (see notes vol 2 p 37)
%  1 March 2010

clear all;
addpath('/home/ig/neef/MFiles/netcdf/');


load 'aam_emac_1960_2000_noIB.mat'
Xm_reg = Xm;
MJD_reg = MJD;

load 'aam_emac_1960_2000_withIB.mat'
Xm_IB = Xm;
MJD_IB = MJD;

IBcol = rand(3,1);

figure(1),clf
  IBplot = plot(MJD_IB,Xm_IB(1,:),'LineWidth',2,'Color',IBcol);
  hold on
  REGplot = plot(MJD_reg,Xm_reg(2,:),'k','LineWidth',2);
  lhandle = [IBplot(1) REGplot(1)];
  legend(lhandle,'With IB','Without IB')
  ylabel('X1')
  title('AEF-1 Mass Terms')

figure(2),clf
  IBplot = plot(MJD_IB,Xm_IB(2,:),'LineWidth',2,'Color',IBcol);
  hold on
  REGplot = plot(MJD_reg,Xm_reg(2,:),'k','LineWidth',2);
  lhandle = [IBplot(1) REGplot(1)];
  legend(lhandle,'With IB','Without IB')
  ylabel('X2')
  title('AEF-2 Mass Terms')



