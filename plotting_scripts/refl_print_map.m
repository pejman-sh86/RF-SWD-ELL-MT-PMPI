function []=refl_print_map(filebase,kmx);

  infile = strcat(filebase,'_sample.mat');
  outfile = strcat(filebase,'_map.dat');
  load(infile);
  AS = A;
  clear A;

  for i=min(AS(:,4)):max(AS(:,4));
    idx=find(AS(:,4)==i);
    [logLmx,ilogLmx]=max(AS(idx,1));
    disp([i,logLmx,idx(ilogLmx),length(idx)]);
    clear idx;
  end;

  idx=find(AS(:,4)==kmx);
  [xmap,imap]=max(AS(idx,1));
  map = AS(idx(imap),4:end-6);
  idxa = [5:4:(map(1)+1)*4+5];idxa(end)=idxa(end)-1;
%  map(2:28)
%  map(find(map(idxa)<=-3))=-2.;
%  map(2:28)
%  if();map();end;
  save(outfile,'map','-ascii');

disp('Updated MAP file');

return;
