function []=rf_print_map2(filebase)

  infile = strcat(filebase,'_voro_sample.mat');
  infile2 = strcat(filebase,'_postlog.dat');

  load(infile);
  AS = A;
  clear A;
  logL2 = load(infile2,'-ascii');
  AS(:,1) = logL2(:,1);

  for i=min(AS(:,4)):max(AS(:,4))
    idx=find(AS(:,4)==i);
    [logLmx,ilogLmx]=max(AS(idx,1));
    disp([i,logLmx,ilogLmx,length(idx)]);
    clear idx;
  end

return