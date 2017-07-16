% make the final plots of spectra for the aam paper
%


% (1) comparison of ERA-40 AAM, OAM, and HAM to GEO

plot_type = 4;
spect_type = 1;
BW = 0;
terms = [1,2];
annotate = 1;

comp = 1;
plot_letter = '(a) \chi_1';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 2;
plot_letter = '(b) \chi_2';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 3;
plot_letter = '(c) \chi_3';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

% (2) comparison of ERA-40 and EMAC runs

plot_type = 1;
annotate = 0;

comp = 1;
terms = 2;
plot_letter = '(a) \chi_1 Mass Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 2;
terms = 2;
plot_letter = '(c) \chi_2 Mass Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 1;
terms = 1;
plot_letter = '(b) \chi_1 Wind Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 2;
terms = 1;
plot_letter = '(d) \chi_2 Wind Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)


comp = 3;
terms = 2;
plot_letter = '(e) \chi_3 Mass Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)

comp = 3;
terms = 1;
plot_letter = '(f) \chi_3 Wind Term';
plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate)


