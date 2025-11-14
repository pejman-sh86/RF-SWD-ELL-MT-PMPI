function [] = plot_rjhist_auv_buck_shear_track(istart);

files=dir('*_sample.mat');
for ifile = istart:length(files);

   disp('plotting');disp(files(ifile).name);
   sample = files(ifile).name;
   [p_mean,p_mead] = plot_rjhist_auv_buck_shear(sample);
   close all;
   compute_buckingham_disp(sample);
   close all;

end;

return;
