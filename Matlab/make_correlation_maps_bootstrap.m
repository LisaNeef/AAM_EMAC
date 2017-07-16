clear all;
runid = 'CCMval';
var_name = 'EF';
total_name = 'EF';
fil_order = 2;

nboot = 500;
alpha = 0.05;


tscale = 3;
comp = 3;
  correlation_map_bootstrap( runid, tscale, comp, var_name, total_name, fil_order, nboot, alpha)
