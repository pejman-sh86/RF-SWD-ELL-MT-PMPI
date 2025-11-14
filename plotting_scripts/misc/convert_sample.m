function [] = convert_sample(sample,iburn)

A = load(sample);
% A = importdata(sample);
matfile = strrep(sample,'.txt','.mat');

A(1:iburn,:) = [];

%save(sample,'A','-ascii');
save(matfile,'A', '-v7.3');

return;
