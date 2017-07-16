function fac = eam_prefactors(comp,var)
% 
% Compute the prefactors needed to get the AAM integral for some variable
% var (= 'U','V', or 'PS') and a given AAM component, into the appropriate 
% AEF units, that is, nothing for X3 and radians for X1 and X2.
%
% NOTE: 
%   Presently I'm taking the prefactors given in Gross (2009).  They aren't
%   necessarily the correct ones, but result from his assumption about love 
%   numbers.
%
% INPUTS:
%   comp: angular momentum component: 'X1', 'X2', or 'X3'
%   var: variable for which weighting is needed: 'U','V', or 'PS'
%
%  Lisa Neef / 9 Dec 2011.
%------------------------------------------------------------------------


%% load constants from the Gross (2009) paper
aam_constants_gross
Re = Re_m;

%% compute the appropriate weighting functions

switch comp
    case {'X1','X2'}
        switch var
            case {'U','V'}
                fac = (-1.591*Re^3)/(Q*g*CminusA);
            case 'PS'
                fac = (-1.098*Re^4)/(g*CminusA);
        end
    case 'X3'
        switch var
            case {'U','V'}
                fac = (0.997*Re^3)/(Q*g*Cm);
            case 'PS'
                fac = (0.748*Re^4)/(g*Cm);
        end
    case 'none'
        fac = 1;
end



