function [] = convert_sample(istart);

files=dir('*_particles.txt');
for ifile = istart:length(files);

   disp('converting');
   disp(files(ifile).name);
   sample = files(ifile).name;
   A = load(sample);
   %A(1:5000,:)=[];
   matfile = strrep(sample,'.txt','.mat');
   save(matfile,'A');

end;

return;
