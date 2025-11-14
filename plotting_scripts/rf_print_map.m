function []=rf_print_map(filebase,kmx);

  %imt = 1;
  infile = strcat(filebase,'_voro_sample.mat');
  outfile = strcat(filebase,'_map_voro.dat');
  load(infile);
  AS = A;
  clear A;

  for i=min(AS(:,4)):max(AS(:,4));
    idx=find(AS(:,4)==i);
    [logLmx,ilogLmx]=max(AS(idx,1));
    disp([i,logLmx,ilogLmx,length(idx)]);
    clear idx;
  end;

  idx=find(AS(:,4)==kmx);
  [xmap,imap]=max(AS(idx,1));
  map = AS(idx(imap),4:end-6);
   %if imt==1
     %  map(62) = AS(idx(imap),end-11);
   %end
  save(outfile,'map','-ascii');

disp('Updated MAP file');

return;
