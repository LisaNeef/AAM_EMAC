% Make sample plots comparing AAM in EMAC, and that implied by observations, for the NATHAN wesbsite.
% ** Use a fixed fontsize of 27 points when exporting
% 9 June 2010

comp = 3;
dataset = [1,3,4];

% (a) comparison of subseasonal variations
plot_tscales_tseries(comp,1,dataset,1987,1990);
ylabel({'\Delta LOD (ms)'},'FontSize',28,'FontName','times');
title('Subseasonal Variation (30-90 days)','FontSize',28,'FontName','times')
annotation('textbox',[0.02 0.8 0.07 0.18],...
    'String',{'(a)'},'FontSize',28, 'FitBoxToText','on','EdgeColor',[1 1 1],'LineWidth',0);
legend('boxoff')

% (b) comparison of interannual variations
plot_tscales_tseries(comp,5,dataset,1962,1999);
title('Interannual Variation (2-7 years)')
ylabel({'\Delta LOD (ms)'},'FontSize',28,'FontName','times');
annotation('textbox',[0.02 0.8 0.07 0.18],...
    'String',{'(b)'},'FontSize',28, 'FitBoxToText','on','EdgeColor',[1 1 1],'LineWidth',0);
legend('boxoff')

