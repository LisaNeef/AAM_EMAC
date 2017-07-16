clear all;

plot_type = 2;
terms = 3;
error_case = 1;
tscale = 3;


% X2 plot for the paper
plot_label = '(a) \chi_2 Interannual Variation';
add_SOI = 0;
comp = 2;
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label)


% X3 plot for the paper
plot_label = '(b) \chi_3 Interannual Variation';
comp = 3;
add_SOI = 1;
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label)

