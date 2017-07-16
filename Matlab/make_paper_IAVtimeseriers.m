% Plot the timeseries of AEFs on the interannual timescale -- final plot for paper.
% stated 12 May 2011

comp = 3;
tscale = 3;
plot_type = 2;
add_SOI = 0;
terms = 3;
error_case = 1; % set to 1 to compare two residual estiamtes to get the error.

%for comp = 1:2
%  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms,error_case)
%end


comp = 3;
  add_SOI = 1;
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms,error_case)

