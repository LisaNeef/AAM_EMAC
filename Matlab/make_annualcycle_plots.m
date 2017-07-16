% commands to make the annual cycle plots for the AAM diagnosis paper.

%---plots with comparison to obs.
comp_geo = 1;
for ii = 1:3
    switch ii
        case 1
            plot_label = '(a) \chi_1';
        case 2
            plot_label = '(b) \chi_2';
        case 3
            plot_label = '(c) \chi_3';
    end
    plot_aam_ann_cycle_compare(ii,3,comp_geo,plot_label);
end

comp_geo = 0;
for ii = 1:3
    switch ii
        case 1
            plot_label_m = '(a) \chi_1 Mass Term';
            plot_label_w = '(b) \chi_1 Wind Term';
        case 2
            plot_label_m = '(c) \chi_2 Mass Term';
            plot_label_w = '(d) \chi_2 Wind Term';
        case 3
            plot_label_m = '(e) \chi_3 Mass Term';
            plot_label_w = '(f) \chi_3 Wind Term';
    end
    plot_aam_ann_cycle_compare(ii,1,comp_geo,plot_label_w);
    plot_aam_ann_cycle_compare(ii,2,comp_geo,plot_label_m);
end
