% fit_annual.m-----------------
function F = fit_annual(tscale)
%
% Fit an annual cycle to AAM / PM / LOD, then compare obs (G) to reanalysis (R) EMAC (M) and 
% obs with OAM and HAM taken out (N).
% started 9 April 2010
% see notes vol 2 p 71
%
% Mods:
%  14 Apr 2010: make this into a function
%
% INPUTS:
%  tscale: choose 'Annual' or 'Semi'
%
% OUTPUTS:
%   structure F containing prograde and retrograde amplitudes and phases, for 
%   datasets 'EMAC','ERA','OBS','OBS-OAM'
%---------------------------------

% use aam_compare_tscales function to load the annual cycles of G, N, R, and M

if tscale == 'A' 
  [t,M,O,R,G,H] = aam_compare_tscales(360,370,'JPL','ERA');
  %[t,M,O,R,G,H] = aam_compare_tscales(230,450,'IERS','ERA');
end
if tscale == 'S' 
  [t,M,O,R,G,H] = aam_compare_tscales(180,185,'JPL','ERA');
  %[t,M,O,R,G,H] = aam_compare_tscales(150,230,'IERS','ERA');
end
nt = length(t);
N = G-O;	% net = geo - ocean

Y = zeros(4,3,length(t));
for ii = 1:3, Y(:,ii,:) = [M(ii,:); R(ii,:); G(ii,:); N(ii,:)]; end
names = {'EMAC','ERA','OBS','OBS-OAM'};

% choose a starting points for the fits
% TODO - set starting values for LOD!
Xm0 = 0;	% mean
m0 = 0;		% trend in mas/day
A0 = 20;	% prograde and retrograde amplitudes in mas
ph_0 = 0;	% prograde and retrograde phases in radians
if tscale =='S'
  A0 = 5;	% prograde and retrograde amplitudes in mas
end

p0 = [Xm0, m0, A0, ph_0, A0, ph_0];

% define some output arrays
Yfit = Y*0;		% fited annual curves - 3 components, 4 datasets
Pfit = zeros(4,6);	% fitted means, trends, amplitudes and phases

%---- compute least squares fit for each component and each dataset

% define the optimization problem
problem.x0 = p0;
problem.xdata = t;
problem.solver = 'lsqcurvefit';
problem.options = optimset;
%problem.options.MaxFunEvals = 1e4;
%problem.options.JacobPattern = sparse(ones(size(t,2),6));
%problem.options.ActiveConstrTol = sqrt(eps);
problem.lb = [0,-5, 0, -pi, 0,-pi ];
problem.ub = [50,5, 50,pi , 50,pi ];
if tscale == 'A', problem.objective = @annual; end
if tscale == 'S', problem.objective = @semiannual; end

% cycle through 4 datasets and optimize
for iset = 1:4
  problem.ydata = squeeze(Y(iset,1:2,:));
  pfit = lsqcurvefit(problem);
  if tscale == 'A', Yfit(iset,1:2,:) = annual(pfit,t); end
  if tscale == 'S', Yfit(iset,1:2,:) = semiannual(pfit,t); end
  Pfit(iset,:) = pfit;
end

% fill into a structure, also converting phases to angles

rad2ang = 180/pi;

F.names = names;
F.emac.true= squeeze(Y(1,1:2,:));
F.emac.fit= squeeze(Yfit(1,1:2,:));
F.era.true = squeeze(Y(2,1:2,:));
F.era.fit = squeeze(Yfit(2,1:2,:));
F.geo.true = squeeze(Y(3,1:2,:));
F.geo.fit = squeeze(Yfit(3,1:2,:));
F.net.true = squeeze(Y(4,1:2,:));
F.net.fit = squeeze(Yfit(4,1:2,:));
F.Ap = Pfit(:,3);
F.Ar = Pfit(:,5);
F.Pp = Pfit(:,4)*rad2ang;
F.Pr = Pfit(:,6)*rad2ang;
F.t  = t;
F.m = Pfit(:,2);
F.Xm = Pfit(:,1);


