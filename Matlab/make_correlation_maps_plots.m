clear all;
runid = 'CCMval';

plot_dir = '/home/ig/neef/Documents/Plots/aam_run_comparison/';

tscale = 3;
comp = 1;
term = 2;
extras = 1;
abs_val =  1;
fil_order = 2;
plot_label = '';

figure(term),clf
plot_correlation_maps_bootstrap(runid,tscale,comp,term,abs_val,fil_order,plot_label,0,extras)
