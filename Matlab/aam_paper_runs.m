%  function [runs,names,ef,no_ib] = aam_paper_runs;
% this function loads the names of the available datasets for the AAM paper, 
% and some of their metadata.
%
% MODS:
%  2 Mar 2011: change  the vector ef to be a flag that holds the various run types:
%	1: emac runs
%	2: ERA reanalysis sets
%	3: geo data
function [runs,names,ef,no_ib] = aam_paper_runs;

r1  = 'ref2_T31L39';             n1  = 'EMAC31';
r2  = 't7_T42L39';               n2  = 'EMAC42';
r3  = 'ref_T63L39';              n3  = 'EMAC63';
r4  = 'CCMval';                  n4  = 'CCMVal';
r5  = 'ERA40';                   n5  = 'ERA-40';
r6  = 'ERAinterim';              n6  = 'ERA-In';
r7  = 'OAM40';                   n7  = 'OAM-40';
r8  = 'HAM40';                   n8  = 'HAM-40';
r9  = 'CAMs4';                   n9  = 'CAM-S4';
r10 = 'GEO';                     n10 = 'ERP';
r11 = 'NET1';                    n11 = 'OBS';		% primary observed AEF estimate
r12 = 'CCMval';                  n12 = 'CCMval-noIB';
r13 = 'OAMin';                   n13 = 'OAM-in';
r14 = 'HAMin';                   n14 = 'HAM-in';
r15 = 'CAMt4';                   n15 = 'CAM-T4';
r16 = 'CAMs1';                   n16 = 'CAM-S1';
r17 = 'CAMt1';                   n17 = 'CAM-T1';
r18 = 'OAMec';			 n18 = 'OAM-ECCO';
r19 = 'NET2';                    n19 = 'OBS. AEF';		% another observed AEF estimate
r20 = 'UNC';			 n20 = 'Obs. error';

runs= {r1;r2;r3;r4; r5; r6; r7; r8; r9; r10; r11; r12; r13; r14; r15; r16; r17; r18; r19; r20};
names= {n1;n2;n3;n4; n5; n6; n7; n8; n9; n10; n11; n12; n13; n14; n15; n16; n17; n18; n19; n20};

nruns = length(runs);

ef = zeros(1,nruns);
  emac = [1,2,3,4,12];
  ef(emac) = 1;

  era = [5,6,7,8,13,14,18];
  ef(era) = 2;

  geo = [10,11,19];
  ef(geo) = 3;

  unc = 20;
  ef(unc) = 4;

  no_ib = zeros(1,nruns);
  no_ib(12) = 1;


