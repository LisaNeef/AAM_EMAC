fil_order = 2;
runid = 'CCMval';

% for local-global correlations in aam

var_name = 'EF';
total_name = 'EF';
for tscale = 3
  for comp = 1
    correlation_map( runid, tscale, comp, var_name, total_name, fil_order );
  end
end


% for correlations between local aam contributions and observations.
var_name = 'EF';
total_name = 'ERP';
for tscale = 2
  for comp = 1
    correlation_map( runid, tscale, comp, var_name, total_name, fil_order );
  end
end


