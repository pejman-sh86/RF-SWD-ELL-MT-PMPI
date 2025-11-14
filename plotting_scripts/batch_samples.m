function []=batch_samples();

files=dir('*2400_sample.mat');
for ifile=1:length(files);
  filename = files(ifile).name;
  filename
  plot_rjhist_auv(filename);
  close all;
end;


return;
