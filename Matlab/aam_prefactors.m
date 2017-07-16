% compute the prefactors for the AAM excitation functions as given in Dobslaw et al 2010
% 23 July 2010
% see notes vol. 3 ps. 40 & 54


% these are the constants in Henryk's paper:

k2 = 0.295;		% deg. 2 rotational love number
ks = 0.938		% secular (fluid limit) Love number
kl = -0.301;		% load Love number
C  = 8.0365e37;		% 3,3 component of Earth inertia tensor
A  = 8.0101e37;		% 1,1 component of Earth inertia tensor
Cm = 7.1237e37;		% 3.3 component of mantle inertia tensor
Am = 7.0999e37;		% 1,1

% numerators and denominators of the prefactors (eq 5 in Dobslaw et al)

num = 1+kl;
den1 = 1-k2/ks;
den3 = 1+(4*k2)/(3*ks)*(C-A)/C;

pref12_I = num/den1;
pref3_I  = num/den3;
pref12_h = 1/den1;
pref3_h  = 1/den3;

disp(['Prefactor equatorial I:   ' num2str(pref12_I)])
disp(['Prefactor axial I:        ' num2str(pref3_I)])
disp(['Prefactor equatorial h:   ' num2str(pref12_h)])
disp(['Prefactor axial h:        ' num2str(pref3_h)])
