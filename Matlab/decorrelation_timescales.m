% This code uses the CCMVal data to roughly estimate the decorrelation time
% for each timescale on which I focus in the AAM paper.
% started 31 Aug 2011.


%---settings that stay the same
date_start = [1962,1,1];
date_stop = [1999,1,1];
fo = 2;		% greater than 2 seems to not work for ERA data
ft = 1;			% butterworth filter




%---cycle through timescales, load the filtered data, compute
%decorreltation time,

for tscale = 1:5
    
    switch tscale
        case 1
           T = [30,90];             % subseasonal
           date_start = [1997,1,1]; 
           date_stop = [1999,12,31];
           axx = [-30 30 -1 1];
         
        case 2
           Tm = [23,34];            % QBO
           Ty = Tm*(1/12); T = Ty*365; 
            axx = [-T(1) T(1) -1 1];
       case 3
           Ty = [2,7]; T = Ty*365;  % IAV
            axx = [-T(1) T(1) -1 1];
        case 4                      % long-term
           Ty = [7,20]; T = Ty*365;
            axx = [-T(2) T(2) -1 1];
        case 5
            T = [230,450];          % annual
            axx = [-T(1) T(1) -1 1];
    end

    [X,XF,MJD] =  retrieve_AAM_filtered('CCMval',T,ft,fo,date_start,date_stop, 1,0,'p');
    
    x1 = squeeze(XF(3,1,:));
    x3 = squeeze(XF(3,3,:));
    
    [p1,lags] = xcorr(x1,x1,'coeff');
    [p3,lags] = xcorr(x3,x3,'coeff');
    
    figure(tscale)
    subplot(1,2,1)
      plot(lags,p1,lags,lags*0,'LineWidth',2),title('\chi_1')
      grid minor
      axis(axx)
    subplot(1,2,2)
      plot(lags,p3,lags,lags*0,'LineWidth',2),title('\chi_3')
      grid minor
      axis(axx)
end