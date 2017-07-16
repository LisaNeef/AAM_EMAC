%plot_compare_LOD_CAM_simple.m
%
% started 22 Dec 2010


%---general settings
  Ty = [7,20]; T = Ty*365;
  date_start = [1960,1,1];
  date_stop = [1999,12,31];

  fil_order = 2;          % greater than 2 seems to not work for ERA data
  filtype = 1;            % butterworth filter



%---load the relevant timeseries
[X,XF,MJD] = retrieve_AAM_filtered('CAMs4',T,filtype,fil_order,date_start,date_stop, 0,0);
Xcam = squeeze(X(3,:));

[X,XF,MJD] = retrieve_AAM_filtered('GEO',T,filtype,fil_order,date_start,date_stop, 0,0);
Xgeo = squeeze(X(3,:));



%---plot settings





%---make plots 



%---export plots



