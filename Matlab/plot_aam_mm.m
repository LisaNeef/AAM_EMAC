% Comparing the monthly mean AAM terms in EMAC output.
% We want to see what the annual cycle looks like for each term and then compare that
% to how the components change regionally.
% The overlying goal is to see what regions most influence a certain timescale of variation.


clear all;

%---load data processed from an EMAC run
%   (data file generated by aam_emac_volint_regional_mm.m)
load aam_vol_comp_t7_allv2.mat


%---constants and conversions
d2r = pi/180.;         %  degrees to radian
rlat = lat*d2r;
rlon = lon*d2r;

Re = 6.371e6;           % radius of earth in meters
Q = 7.292115e-5;        % rotation rate of the earth in rad/s
g = 9.81;               % grav constant in m/s2
C = 7.1236e37;          % principal moment of inertia of the mantle in kgm2
CminusA = 2.34e35;       % C minus eq. MOI, same units
d2r = pi/180.;         %  degrees to radian

pm12 = -Re^4/(g*CminusA);
pw12 = -1.43*Re^3/(Q*g*CminusA);
pm3  = 0.7*Re^4/(g*C);
pw3  = Re^3/(Q*g*C);

rad2microas = (180/pi)*60*60*1e6;
rad2mas = (180/pi)*60*60*1e3;

LOD0 = double(24*60*60*1e3);            % standard LOD in milliseconds.

%---initialize empty arrays
Mass_term = zeros(3,12);
Wind_term = zeros(3,12);

%---cycle though months and do the remaining integrals (over lat and lon)

for im = 1:12

  xxm = squeeze(Xm_mean(:,im,:,:));
  int_xxm_dlat = -trapz(rlat,xxm, 2);
  int_xxm_dlat_dlon = trapz(rlon,int_xxm_dlat, 3);
  Mass_term(:,im) = int_xxm_dlat_dlon;

  xxw = squeeze(Xw_mean(:,im,:,:));
  int_xxw_dlat = -trapz(rlat,xxw, 2);
  int_xxw_dlat_dlon = trapz(rlon,int_xxw_dlat, 3);
  Wind_term(:,im) = int_xxw_dlat_dlon;

end


%---multiply by prefactors to get friendly units (mas for X1 and X2, ms for X3)

Xm = zeros(3,12);
Xw = zeros(3,12);

Xm(1:2,:) = rad2mas*pm12*Mass_term(1:2,:);
Xw(1:2,:) = rad2mas*pw12*Wind_term(1:2,:);

Xm3_anom = Mass_term(3,:) - mean(Mass_term(3,:));
Xw3_anom = Wind_term(3,:) - mean(Wind_term(3,:));

Xm(3,:) = LOD0*pm3*Xm3_anom;
Xw(3,:) = LOD0*pw3*Xw3_anom;


%---plot 12-month cycles for 6 aam terms!

LW = 3;
ifig = 1;
TT = {'\chi_1 Mean Annual Cycle','\chi_2 Mean Annual Cycle','\Delta LOD Mean Annual Cycle'};

for ii = 1:3

  figure(ii), clf
  m = plot(1:12,Xm(ii,:),'k','LineWidth',LW);
  hold on
  w = plot(1:12,Xw(ii,:),'k--','LineWidth',LW);
  t = plot(1:12,Xw(ii,:)+Xm(ii,:),'k-','LineWidth',LW,'Color',0.7*[1,1,1]);
  title(TT(ii))
  xlabel('Month')
  if ii < 3, ylabel('mas'), else ylabel('ms'), end
  legend([m(1) w(1) t(1)], 'Mass Term','Wind Term','Total',0)

end



