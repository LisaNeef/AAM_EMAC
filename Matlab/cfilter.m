function [out] = cfilter(d,fs,fcol,fcoh,odn,unit)

% INPUT
% d         data to be filtered
% fs        sample frequency
% fcol      cut off wave length, lower limit 
%               (inverse of cut off frequency low pass)
% fcoh      cut off wavelength , upper limit 
%               (inverse of cut off frequency high pass)
% odn       order n filer
% unit      unit of the time axis

%        lower f | higher f 
% ---------------f--------------
%                |
% <--------------| low pass
%                |
%      high pass |------------>

%  band pass 
%   high pass |-------| low pass


% fcol = [7]
% fcoh = [14]
% fs = 12
% odn = 3
% unit='year'


t = [0:length(d)-1]/fs;

if isempty(fcol) 
    ftype = 'high';
    fil = inv(fcoh);
    tit = ['pass ','0 - ',num2str(fcoh)];
elseif isempty(fcoh)
    ftype = 'low';
    fil = inv(fcol);
    tit = ['pass ',num2str(fcol),' - Inf'];
elseif fcoh > fcol
    ftype = 'bandpass';
    fil = [inv(fcoh) inv(fcol)];
    tit = ['pass ',num2str(fcol),' - ',num2str(fcoh)];
else
    error('error in cut off frequencies')
end

%%
[b a] = butter(odn,2*fil./fs,ftype);
x = filtfilt(b,a,d);
out(:,1) = x;

%h = dfilt.df2(b,a);
%hc = cascade(h,h); % hc magnitude = h magnitude squared
%hfvt = fvtool(h,hc,'FrequencyScale','log');
%legend(hfvt,'24 dB/Octave','48 dB/Octave')

%figure
%plot(t,d)
%hold on
%plot(t,x,'r','linewidth',2)
%legend('Input','Filtered');xlabel(unit);grid on
%title(['Butterworth filter design ',ftype,' pass; ',tit,' ',unit])



%%
[b,a] = cheby1(odn,1,2*fil./fs,ftype);            % Chebyshev Type I
x = filtfilt(b,a,d);
out(:,2) = x;

%figure
%plot(t,d)
%hold on
%plot(t,x,'r','linewidth',2)
%legend('Input','Filtered');xlabel(unit);grid on
%title(['Chebyshev Type I filter design ',ftype,' pass; ',tit,' years'])

%h = dfilt.df2(b,a);
%hc = cascade(h,h); % hc magnitude = h magnitude squared
%hfvt = fvtool(h,hc,'FrequencyScale','log');
%legend(hfvt,'24 dB/Octave','48 dB/Octave')


%%
[b,a] = cheby2(odn,10,2*fil./fs,ftype);          % Chebyshev Type II
x = filtfilt(b,a,d);
out(:,3) = x;

%figure
%plot(t,d)
%hold on
%plot(t,x,'r','linewidth',2)
%legend('Input','Filtered');xlabel(unit);grid on
%title(['Chebyshev Type II filter design ',ftype,' pass; ',tit,' years'])

%h = dfilt.df2(b,a);
%hc = cascade(h,h); % hc magnitude = h magnitude squared
%hfvt = fvtool(h,hc,'FrequencyScale','log');
%legend(hfvt,'24 dB/Octave','48 dB/Octave')




