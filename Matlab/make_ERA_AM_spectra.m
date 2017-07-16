% plot spectra of the ERA-40 AM contributors (AAM, OAM, HAM), compared to
% the total, scaling PSD by frequency, so that different features are more
% visible.
% this image is intended for talks and stuff.


%% figure settings

clear all;
annotate = 0;
terms = 3;
plot_letter = ' ';
plot_type = 4;
spect_type = 1;
BW = 0;
extra_labels = 1;
comp = 3;
export = 0;
label_x_axis = 1;
label_y_axis = 1;


%% plot code

figure(1),clf

plot_aam_spectra(comp,terms,plot_type,spect_type,BW,plot_letter,annotate,export,extra_labels,label_x_axis,label_y_axis)


%% export

plot_dir = '/home/ig/neef/Documents/Plots/ERA/';
fig_name = [plot_dir,'spectrum_X',num2str(comp),'.png'];

LW = 3;
ph = 6;        % paper height
pw = 17;        % paper width
fs = 20;        % fontsize

exportfig(1,fig_name,'width',pw,'height',ph,'fontmode','fixed', 'fontsize',fs,'color','cmyk','LineMode','fixed','LineWidth',LW,'format','png');





  
