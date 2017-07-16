% this code compares the power spectrum of ERA-40 atmospheric excitation
% functions, before and after a certain date.
% The goal is to see whether the spectrum has changed as the ERA-40
% analysis becomes filled with better observations.
%
% Lisa Neef
% started 23 June 2011.




%---temporary inputs:
comp = 2;               % vector component of AAM to plot.
cutoff_year = 1980;     % year where we separate
term = 1;               % which terms to plot?  1 = wind, 2 = mass, 3 = total

%---load ERA data

[Xw,Xm,mjd] = read_EFs('aam','ERA40',1);

switch term
    case 1
        X = Xw(comp,:);
    case 2
        X = Xm(comp,:);
    case 3
        X = Xw(comp,:)+Xm(comp,:);
end


%---separate timeseries into pre and post 1980.

mjd_cutoff = date2mjd(cutoff_year,1,1,0,0,0);
old = find(mjd <= mjd_cutoff);
new = find(mjd > mjd_cutoff);

X_old = X(old);
X_new = X(new);


%---compute frequency spectra.

fs = 365;               % sampling frequency
n_old = length(old);    % Length of signal 1
n_new = length(new);    % Length of signal 2

%t_old = (0:n_old-1)/fs;         % Period vector in years.
%t_new = (0:n_new-1)/fs;         % Period vector in years.

NFFT_old = 2^nextpow2(n_old); % Next power of 2 from length of y
NFFT_new = 2^nextpow2(n_new); % Next power of 2 from length of y

P_old = fft(X_old,NFFT_old)/n_old;
P_new = fft(X_new,NFFT_new)/n_new;

f_old = fs/2*linspace(0,1,NFFT_old/2+1);
f_new = fs/2*linspace(0,1,NFFT_new/2+1);

t_old = f_old.^-1;              % period vector (years)
t_new = f_new.^-1;              % period vector

S1 = 2*abs(P_old(1:NFFT_old/2+1));
S2 = 2*abs(P_new(1:NFFT_new/2+1));


%---plot settings

lh = zeros(1,2);
col_old = 0*ones(1,3);
col_new = rand(1,3);

% plot limits in years.
tmin = 14/365;      
tmax = 15/365;

xmin = 10^-2;
xmax = 10;
ymax = 10^-8;
ymin = 10-13;

%---plot!

figure(1),clf
    lh(1) = loglog(t_old,S1,'Color',col_old);
    hold on
    lh(2) = loglog(t_new,S2,'Color',col_new);

    legend(lh,'ERA-40 Old','ERA-40 New')
    xlabel('Period (y)')
%    axis([xmin,xmax,ymin,ymax])
grid on

%---export.