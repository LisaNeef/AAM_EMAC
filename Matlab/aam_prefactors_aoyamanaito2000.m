% compute the prefactors for the AAM excitation functions as given in Aoyama & Naito, 2000
% 23 Aug 2010
% see notes vol. 3 p. 59

clear all;

% these are the constants given in the Aoyama and Naito (2000) paper:

k2p 	= -0.310;	% 2nd order load Love nr for whole earth 
k2pm 	= -0.245;	% 2nd order load Love nr for mantle
k2	= -0.348;	% 2nd order rotational Love nr
kr	= 0.998;	% correction for axial component due to rotational deformation.

% additional values not given in this paper, taken as best I could from Dobslaw et al 2010
k 	= 0.938         
C  = 8.0365e37;         % 3,3 component of Earth inertia tensor
A  = 8.0101e37;         % 1,1 component of Earth inertia tensor
Cm = 7.1237e37;         % 3.3 component of mantle inertia tensor
Am = 7.0999e37;         % 1,1


% numerators and denominators of the prefactors (eqs A4a-A4b in Aoyama et al)

num12I = A*k*(1+k2p);
den12I = k-k2;
num12h = k;
den12h = k-k2;

pref12_I = num12I/den12I;
pref12_h = num12h/den12h;

pref3_I = kr*(1+k2pm);
pref3_h = kr;


disp(['Prefactor equatorial I:   ' num2str(pref12_I)])
disp(['Prefactor axial I:        ' num2str(pref3_I)])
disp(['Prefactor equatorial h:   ' num2str(pref12_h)])
disp(['Prefactor axial h:        ' num2str(pref3_h)])
