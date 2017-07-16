% make plots of the correlations between models and the geodetic residual

exp_type = 1;

for itscale = 1:3
  plot_correlations_bar_v2(itscale)
end
