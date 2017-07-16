% make slice plots of the variance contributions per latitude and longitude band.

for comp = 1:3 
  for ilon = 0:1
    plot_aam_acvar(comp,ilon,'X')
  end
end

for ilon = 0:1
  plot_aam_acvar(comp,ilon,'u')
end
for ilon = 0:1
  plot_aam_acvar(comp,ilon,'p')
end
