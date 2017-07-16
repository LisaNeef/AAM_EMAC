function col = aam_paper_colors;
%col = aam_paper_colors;
% these are the colors to use in the AAM Diagnosis Paper
% started 26 Nov 2010
% Mods:
%  - 6 Dec 2010: add other colors for additional run options.

col = zeros(19,3);
%col = colormap(jet(17));

col_era = 0*ones(1,3);
col_era_interim = 0.7*ones(1,3);
col_geo = [0.7060    0.0318    0.2769];
col_net = col_geo;

col(1:4,:) 	= colormap(cool(4));
%col(1:4,:) 	= colormap(jet(4));
col(5:8,:) 	= colormap(bone(4));		% ERA-40 and ERA-interim, incl. OAM and HAM
col(10,:) 	= col_geo;
col(11,:) 	= col_net;
col(12,:) 	= 0.7*col(4,:);
col(13:14,:) 	= 0.7*(col(7:8,:));		% OAM and HAM estimates from ERA Interim
cam = [9,15:17];
col(cam,:) 	= colormap(summer(4));		% CAM estimates
col(18,:)	= col(7,:);			% other OAM estimate
col(19,:) 	= col_geo;			% the other AEF residual estimate

%col(1:4,:) = colormap(summer(4));		% CCMVal runs
%col(5,:) = col_era;
%col(6,:) = col_era_interim;
%col(7:8,:) = colormap(autumn(2));
%col(12:14,:) = colormap(jet(3));


%--this is the order of runs and colors:

%r1 = 'ref2_T31L39';             n1 = 'EMAC31';
%r2 = 't7_T42L39';               n2 = 'EMAC42';
%r3 = 'ref_T63L39';              n3 = 'EMAC63';
%r4 = 'CCMval';                  n4 = 'CCMVal';
%r5 = 'ERA40';                   n5 = 'ERA-40';
%r6 = 'ERAinterim';              n6 = 'ERA-In';
%r7 = 'OAM40';                   n7 = 'OAM-40';
%r8 = 'HAM40';                   n8 = 'HAM-40';
%r9 = 'CAMs4';                   n9 = 'CAM';
%r10= 'GEO';                     n10= 'GEO';
%r11= 'NET';                     n11= 'Residual';
%r12= 'CCMval';                  n12= 'CCMval-noIB';
%r13= 'OAMin';                   n13= 'OAM-in';
%r14= 'HAMin';                   n14= 'HAM-in';
%r15= 'CAMt4';                   n15= 'CAM';
%r16= 'CAMs1';                   n16= 'CAM';
%r17= 'CAMt1';                   n17= 'CAM';
%r18 = 'OAMec';                   n18 = 'OAM-ECCO';
%r19 = 'NET2';                    n19 = 'OBS. AEF';              % another observed AEF estimate

