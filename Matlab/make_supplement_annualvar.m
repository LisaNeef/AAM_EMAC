% Make the plots for the supplementary material -- annual cycle variance in
% wind components.
%
% 11 July 2011.
%-----------------------------------------------

clear all;

% (1) latitude slices.

comp = 3;
term = 1;
lonplot = 0;

var_to_plot = 'u';
text = '(a) u (300 hPa) Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp = 1;
var_to_plot = 'X';
text = '(b) \chi_1 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp =  2;
text = '(c) \chi_2 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp = 3;
text = '(d) \chi_3 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);



% (2) longitude slices.

comp = 3;
term = 1;
lonplot = 1;

var_to_plot = 'u';
text = '(a) u (300 hPa) Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp = 1;
var_to_plot = 'X';
text = '(b) \chi_1 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp =  2;
text = '(c) \chi_2 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);

comp = 3;
text = '(d) \chi_3 Wind Term Variance';
plot_aam_acvar(comp,term,lonplot,var_to_plot,text);
