%----phasor_fit_annual.m---------------
%
% Compute fits to the annual and semiannual componetns of AAM in model and obs
% and plot as phasor diagrams
% (following the example of Gross et al 1993)
% started 15 April 2010 - see notes vol 2 p. 76

clear all;

% load data, or compute prograde and retrograde components using fit_annual.m
fA = 'fit_ann.mat';
fS = 'fit_sem.mat';
if exist(fA) == 0
  A = fit_annual('A');
  save 'fit_ann.mat'
else
  load 'fit_ann.mat'
end

if exist(fS) == 0
  S = fit_annual('S');
  save 'fit_sem.mat'
else
  load 'fit_sem.mat'
end


% convert to cartesian coordinates. Compare EMAC (M), ERA (R) and GEO-HAM-OAM (N)

Ax = zeros(2,4);	% x amplitudes (prograde, retrograde)
Ay = zeros(2,4); % y amplitudes (prograde, retrograde)
Sx = zeros(2,4);	% x amplitudes (prograde, retrograde)
Sy = zeros(2,4); % y amplitudes (prograde, retrograde)

for s = 1:4
  [Ax(1,:),Ay(1,:)] = pol2cart(A.Pp,A.Ap);
  [Ax(2,:),Ay(2,:)] = pol2cart(A.Pr,A.Ar);
  [Sx(1,:),Sy(1,:)] = pol2cart(S.Pp,S.Ap);
  [Sx(2,:),Sy(2,:)] = pol2cart(S.Pr,S.Ar);
end

% plot formatting stuff
TA = {'Annual Prograde', 'Annual Retrograde'};
TS = {'Semiannual Prograde', 'Semiannual Retrograde'};
LW = 4;
HW = .06;
HH = .15;
LH = zeros(4,2);
Mcol = [0.7513    0.2551    0.5060];
Rcol = [0.0714    0.5216    0.0967];
Gcol = zeros(1,3);
Ncol = 0.7*ones(1,3);
gr = 0.6*ones(3,1);

c = rand(4,3);
  c(1,:) = Mcol; 
  c(2,:) = Rcol;
  c(3,:) = Gcol;
  c(4,:) = Ncol;

figure(5),clf
for j = 1:2
  subplot(2,1,j)
  plot([0,Ax(j,1)],[0,Ay(j,1)],'Color',ones(1,3))
  hold on
  Mplot = plot_arrow(0,0,Ax(j,1),Ay(j,1),'Color',c(1,:),'LineWidth',LW,'facecolor',c(1,:),'headwidth',HW,'headheight',HH,'edgecolor',c(1,:));
  Rplot = plot_arrow(0,0,Ax(j,2),Ay(j,2),'Color',c(2,:),'LineWidth',LW,'facecolor',c(2,:),'headwidth',HW,'headheight',HH,'edgecolor',c(2,:));
 % Gplot = plot_arrow(0,0,Ax(j,3),Ay(j,3),'Color',c(3,:),'LineWidth',LW,'facecolor',c(3,:),'headwidth',HW,'headheight',HH,'edgecolor',c(3,:));
  Nplot = plot_arrow(0,0,Ax(j,4),Ay(j,4),'Color',c(4,:),'LineWidth',LW,'facecolor',c(4,:),'headwidth',HW,'headheight',HH,'edgecolor',c(4,:));
  LH = [Mplot(1),Rplot(1),Nplot(1)];
  legend(LH,A.names,'Location','Best')
  legend boxoff
  title(TA(j))
  xlabel('Real Part (mas)')
  ylabel('Imaginary Part (mas)')
  grid on
  %axis([-20 20 -20 20])
end

figure(6),clf
for j = 1:2
  subplot(2,1,j)
  plot([0,Sx(1,1)],[0,Sy(1,1)],'Color',ones(1,3))
  hold on
  Mplot = plot_arrow(0,0,Sx(j,1),Sy(j,1),'Color',c(1,:),'LineWidth',LW,'facecolor',c(1,:),'headwidth',HW,'headheight',HH,'edgecolor',c(1,:));
  Rplot = plot_arrow(0,0,Sx(j,2),Sy(j,2),'Color',c(2,:),'LineWidth',LW,'facecolor',c(2,:),'headwidth',HW,'headheight',HH,'edgecolor',c(2,:));
 % Gplot = plot_arrow(0,0,Sx(j,3),Sy(j,3),'Color',c(3,:),'LineWidth',LW,'facecolor',c(3,:),'headwidth',HW,'headheight',HH,'edgecolor',c(3,:));
  Nplot = plot_arrow(0,0,Sx(j,4),Sy(j,4),'Color',c(4,:),'LineWidth',LW,'facecolor',c(4,:),'headwidth',HW,'headheight',HH,'edgecolor',c(4,:));
  %LH = [Mplot(1),Rplot(1),Gplot(1),Nplot(1)];
  LH = [Mplot(1),Rplot(1),Nplot(1)];
  legend(LH,A.names,'Location','Best')
  legend boxoff
  title(TS(j))
  xlabel('Real Part (mas)')
  ylabel('Imaginary Part (mas)')
  grid on
  %axis(7*[-1 1 -1 1])
end



