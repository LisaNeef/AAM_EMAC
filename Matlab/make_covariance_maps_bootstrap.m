clear all;
runid = 'CCMval';
var_name = 'EF';
total_name = 'EF';
fil_order = 2;

nboot = 500;
alpha = 0.05;


tscale = 2;
comp = 3;
for tscale = [2,1,3]
    for comp = [3,1,2]
          covariance_map_bootstrap( runid, tscale, comp, var_name, total_name, fil_order, nboot, alpha)
    end
end
