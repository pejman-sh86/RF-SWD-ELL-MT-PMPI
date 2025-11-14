function []=ffi_print_map(filebase,kmx);

  parfile = strcat(filebase,'_parameter.dat');
  infile = strcat(filebase,'_sample.mat');
  datfile = strcat(filebase,'.hdf5');
  datfilecGPS = strcat(filebase,'_cGPS.hdf5');
  outfile = strcat(filebase,'_map.dat');
  [IMAP,ICOV,I_WP,I_cGPS,I_GPS,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
   NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile);

  NFPMX = NVMX*NPV;

  load(infile);
  NSTN = h5readatt(datfile,'/Observed_data','N_sta');
  NcGPS = h5readatt(datfilecGPS,'/Observed_data','N_sta');
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
  %idx(imap)=10;
  map = AS(idx(imap),5:end-6);
  kmap = AS(idx(imap),4);
  for ivo = 1:kmap;
    mapvoro(ivo,1:NPV) = map((ivo-1)*NPV+1:ivo*NPV);
  end;
  maphyp   = map(NVMX*NPV+1:NVMX*NPV+NMISC);
  mapsd    = map(NVMX*NPV+NMISC+1:NVMX*NPV+NMISC+NSTN);
  mapsdGPS = map(NVMX*NPV+NMISC+NSTN+1:NVMX*NPV+NMISC+NSTN+3);
  mapsdcGPS= map(NVMX*NPV+NMISC+NSTN+3+1:NVMX*NPV+NMISC+NSTN+3+NcGPS);
  mapalp   = map(NVMX*NPV+NMISC+NSTN+NcGPS+3+1);
  save(outfile,'kmap','-ascii');
  save(outfile,'mapvoro','-ascii','-append');
  save(outfile,'maphyp','-ascii','-append');
  save(outfile,'mapsd','-ascii','-append');
  save(outfile,'mapsdcGPS','-ascii','-append');
  save(outfile,'mapsdGPS','-ascii','-append');
  save(outfile,'mapalp','-ascii','-append');

disp('Updated MAP file');

save tmp
return;
