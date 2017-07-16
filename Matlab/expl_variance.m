%--expl_variance.m----------------
%
%  Compute the percentage of the observed variance in EOP that is explained by ERA reanalysis and EMAC
%  where variace explained is VE = (1-var(residual)/var(obs))*100
%  Options:
%	1. Raw EOP observations
%	2. EOPs minus OAM
%	3. EOPs minus OAM and HAM
%
%  See notes vol. 2 p. 79	April 21, 2010

%---file settings and paths

clear all;



%---establish timescales

TT = zeros(5,2);
TT(1,:) = [30,90]; 	% subseasonal.
TT(2,:) = [90,150];	% terannual
TT(3,:) = [150,230];	% semiannual
TT(4,:) = [230,450];	% annual
TT(5,:) = [730,1825];	% interannual
nT = size(TT,1);

% array for explained variances:
EVM = zeros(nT,2);
EVR = zeros(nT,2);
pM = zeros(nT,2);
pR = zeros(nT,2);

%---loop over timescales and compute variance explained in each case

for iT = 1:nT

  % load the filtered timeseries and combine equatorial terms
  % (note that LOD changes are converted to ms)
  [t,M0,O0,R0,G0,H0,MJD] = aam_compare_tscales(TT(iT,1),TT(iT,2),'IERS','ERA');
  M = zeros(2,length(t));	R = M;	G = M;	
  M(1,:) = M0(1,:)+i*M0(2,:);	M(2,:) = 1000*M0(3,:);
  R(1,:) = R0(1,:)+i*R0(2,:);	R(2,:) = 1000*R0(3,:);
  G(1,:) = G0(1,:)+i*G0(2,:);	G(2,:) = 1000*G0(3,:);
  O(1,:) = O0(1,:)+i*O0(2,:);	O(2,:) = 1000*O0(3,:);
  N = G-O;

  % compute residual  
  DM = N-M;
  DR = N-R;

  % compute explained variance for equatorial and axial terms
  EVM(iT,:) = 100*([1,1] - var(DM')./var(N')); 	% expl variance for EMAC
  EVR(iT,:) = 100*([1,1] - var(DR')./var(N'));	% expl variance for ERA

  % also compute the correlation for each
  for ii = 1:2
    [cm,lag] = xcorr(M(ii,:),N(ii,:),'coeff');
    [cr,lag] = xcorr(R(ii,:),N(ii,:),'coeff');
    [am,bm] = max((cm));
    [ar,br] = max((cr));
    pM(iT,ii) = am;
    pR(iT,ii) = ar;
  end
end


