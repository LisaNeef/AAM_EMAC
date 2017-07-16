clear all;

plot_type = 2;
add_SOI = 0;
terms = 3;
error_case = 1;
tscale = 1;


% X2 plot for the paper
plot_label = '(a) \chi_2 Subseasonal Variation';
comp = 2;
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label)


% X3 plot for the paper
plot_label = '(b) \chi_3 Subseasonal Variation';
comp = 3;
  plot_aam_filtered_compare_v2(comp,tscale,plot_type,add_SOI,terms, error_case,plot_label)

