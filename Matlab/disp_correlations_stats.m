function disp_correlations_stats(comp,tscale)
% This function displays, for each dataset available, the correlation to observations, 
% the effective degrees of freedom, the decorrelation timescale, and the minimum correlation
% that's necessary for 95% significance.
% 14 Mar 2011
%
% MODS:
%
% INPUTS:
%	comp: the vector component of excitation
%	tscale: the timescale of interest.
%       tscale = 1: 30-90 day (subseasonal) variations.
%       tscale = 2: 24-30 month variations (QBO).
%       tscale = 3: 2-7 year variation.
%       tscale = 4: 7-20 year variations.


exp_type = 1;
[p,lag,Nsamples,tau,runs,C,L,R2] = correlations_aam_filtered_v2(comp,tscale,exp_type);

if tscale == 1, tscalename = 'Subseasonal'; end
if tscale == 2, tscalename = 'Quasi-Biennial'; end
if tscale == 3, tscalename = 'Interannual'; end
if tscale == 4, tscalename = 'Long-Term'; end

nR = length(p);
Neff = Nsamples./tau;
psig = erfinv(0.95)*sqrt(2./Neff);

disp(['Timescale:     ',tscalename])

disp('Data Set  Decorrelation Time    N_eff      95% Sig. Limit')
for ii = 1:nR
  disp([runs(ii), num2str(tau(ii)),num2str(Neff(ii)),num2str(psig(ii))]);
end



