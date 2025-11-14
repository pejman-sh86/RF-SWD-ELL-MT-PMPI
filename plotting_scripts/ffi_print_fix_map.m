function []=ffi_print_fix_map(filebase,kmx);

  parfile = strcat(filebase,'_parameter.dat');
  infile = strcat(filebase,'_sample.mat');
  datfile = strcat(filebase,'.hdf5');
  outfile = strcat(filebase,'_map.dat');
  [IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
   NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF]=ffi_read_fix_parfile(parfile);

  NFPMX = NVMX*NPV;

  load(infile);
  NSTN = h5readatt(datfile,'/Observed_data','N_sta');
  AS = A;
  clear A;

  [logLmx,idx]=max(AS(:,1))
  map = AS(idx,5:end-6);
  kmap = AS(idx,4);
  for ivo = 1:kmap;
    mapvoro(ivo,1:NPV) = map((ivo-1)*NPV+1:ivo*NPV);
  end;
  maphyp   = map(NVMX*NPV+1:NVMX*NPV+NMISC);
  mapar    = map(NVMX*NPV+NMISC+1:NVMX*NPV+NMISC+NSTN);
  mapsd    = map(NVMX*NPV+NMISC+NSTN+1:NVMX*NPV+NMISC+NSTN+NSTN);
  save(outfile,'kmap','-ascii');
  save(outfile,'mapvoro','-ascii','-append');
  save(outfile,'maphyp','-ascii','-append');
  save(outfile,'mapar','-ascii','-append');
  save(outfile,'mapsd','-ascii','-append');
  
disp('Updated MAP file');

return;
