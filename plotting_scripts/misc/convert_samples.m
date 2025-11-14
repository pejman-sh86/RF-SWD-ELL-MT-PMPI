function [] = convert_sample(istart);

files=dir('*_sample.txt');
for ifile = istart:length(files);

   disp('converting');
   disp(files(ifile).name);
   sample = files(ifile).name;
   A = load(sample);
   nsmp=size(A,1);
   nburnin = round(nsmp/3);
   A(1:nburnin,:) = [];
   nsmp=size(A,1);
   if(nsmp>100000);
     idx = randsample(nsmp,100000);
     idx = sort(idx);
   else
       idx = [1:nsmp];
   end;
   A = A(idx,:);
%   if(nsmp>80000);
%     A=A(1:2:end,:);
%     A(1:end-40000,:)=[];
%   elseif(nsmp>40000);
%     A(1:end-40000,:)=[];
%   end;
   matfile = strrep(sample,'.txt','.mat');
   save(matfile,'A');
   save(sample,'A','-ascii');

end;

return;
