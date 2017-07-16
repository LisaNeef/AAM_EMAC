function Xout = annual(p,t)
% an annual-varying function plus a mean term and trend
%to fit to the annual cycle of polar motion / lod
%
% Output: 
%	Xout:	X1 and X2 equatorial excitation functions  
% Inputs
%	t: 	time vector
%	p: 	function parameters (to be solved for with least squares)
%		p = [Xm, m, Ap, ph_p, Ar, ph_r]
%		Xm = mean state
%		m = trend
%		Ap = prograde component
%		ph_p = prograde phase
%		Ar = retrograde component
%		ph_r = retrograde phase
%----------------------------------------------------------

%--decompose the input vector
Xm = p(1);
m = p(2);
Ap = p(3);
ph_p = p(4);
Ar = p(5);
ph_r = p(6);

w = 2*pi/365;	% annual frequency in rad/day

%---compute the fit function

A = Xm+t*0;			% mean
B = m*t;			% trend
C = Ap*exp(i*ph_p)*exp(i*w*t);  % prograde
D = Ar*exp(i*ph_r)*exp(-i*w*t);  % retrograde

Xreal = real(A+B+C+D);
Ximag = imag(A+B+C+D);

Xout = [Xreal;Ximag];

end
