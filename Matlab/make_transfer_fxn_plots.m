% make plots for the paper supplementary material of the various transfer functions.


v = 'u';
for ii = 1:3
  plot_transfer_functions(ii,v);
end

v = 'v';
for ii = 1:2
  plot_transfer_functions(ii,v);
end

v = 'm';
for ii = 1:3
  plot_transfer_functions(ii,v);
end

