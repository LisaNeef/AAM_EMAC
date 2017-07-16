% Plot the timeseries of AEFs on the interannual timescale -- final plot for paper.
% stated 12 May 2011


tscale = 1;
plot_type = 2;
terms = 3;
add_SOI = 0;

for comp = 1:3
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms)
end
