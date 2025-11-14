function [] = convert_sample(sample)

A = dlmread(sample);
matfile = strrep(sample,'.txt','.mat');
save(matfile,'A','-v7.3');

return;
