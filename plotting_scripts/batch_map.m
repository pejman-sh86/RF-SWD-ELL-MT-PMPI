function []=batch_samples();

files=dir('*2400_sample.mat');
for ifile=1:length(files);
  filename = files(ifile).name;
  mapfile = strrep(filename,'sample.mat','map.dat');
  load(filename);
  [a,b]=max(A(:,1));
  disp([mapfile,num2str(a)]);
  map=A(b,4:end-6);
  save(mapfile,'map','-ascii');
  clear A a b map;
end;


return;
