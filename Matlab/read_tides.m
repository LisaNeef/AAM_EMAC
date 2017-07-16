%----read_tides.m-------------------
%
%  This function returns X1, X2, and LOD change due to Ocean tides of different
%  frequencies from the IERS (SB for Tides) website
%  We read in amplitudes and phases of the prograde and retrograde components 
%  at each amplitude.
%  http://bowie.gsfc.nasa.gov/ggfc/tides/eop.html
%
%  Lisa Neef |  8 Mar 2010
%
%  Notes:
%  - so far just hard-coding in the values from dataset eop_tpxo6
function Xout = read_tides

 
%---structure: tide amplitudes and phases hardcoded

tpxo.name = 'TPXO EOP Data';
tpxo.tides = ['Q1', 'O1', 'P1', 'K1', 'N2', 'N2', 'M2', 'S2', 'K2'];
tpxo.periods = [26.868,25.819,24.066,23.934,12.658,12.421,12.000,11.967];
tpxo.X1.amp = 
