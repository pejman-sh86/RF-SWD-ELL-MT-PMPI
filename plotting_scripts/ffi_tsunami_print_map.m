function []=ffi_print_map(filebase,kmx);

  parfile = strcat(filebase,'_parameter.dat');
  infile = strcat(filebase,'_sample.mat');
  datfile = strcat(filebase,'.hdf5');
  outfile = strcat(filebase,'_map.dat');
  [IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
   NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
   MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

  NFPMX = NVMX*NPV;

  load(infile);
  NSTN = h5readatt(datfile,'/Observed_data','N_sta');
  AS = A;
  clear A;
  j = 1;
  for i=min(AS(:,4)):max(AS(:,4));
    idx=find(AS(:,4)==i);
    k(j) = i;
    if(idx);
      [logLmx(j),ilogLmx(j)]=max(AS(idx,1));
      map = AS(idx(ilogLmx(j)),:);
      idead(j) = length(find(map(5:NFPMX+4)<-99.));
    else
      logLmx(j)  = NaN;
      ilogLmx(j) = NaN;
      idead(j)   = NaN;
    end;
    ialive(j) = i*(NPV-2)-idead(j);
    disp([i,ialive(j),logLmx(j),ilogLmx(j),length(idx)]);
    clear idx;
    j = j + 1;
  end;

  figure;hold on;
  plot(k,logLmx)
  figure;hold on;
  plot(ialive,logLmx,'xk')

  idx=find(AS(:,4)==kmx);
  [xmap,imap]=max(AS(idx,1));
  map = AS(idx(imap),5:end-6);
  kmap = AS(idx(imap),4);
  idead = 0;
  for ivo = 1:kmap;
    mapvoro(ivo,1:NPV) = map((ivo-1)*NPV+1:ivo*NPV);
    idxdead=find(mapvoro(ivo,1:NPV)<-99.);
    idead = idead + length(idxdead);
    clear idxdead;
  end;
  maphyp = map(NVMX*NPV+1:NVMX*NPV+NMISC);
  mapar = map(NVMX*NPV+NMISC+1:NVMX*NPV+NMISC+NSTN);
  mapsd = map(NVMX*NPV+NMISC+NSTN+1:NVMX*NPV+NMISC+NSTN+NSTN);
  save(outfile,'kmap','-ascii');
  save(outfile,'mapvoro','-ascii','-append');
  save(outfile,'maphyp','-ascii','-append');
  save(outfile,'mapar','-ascii','-append');
  save(outfile,'mapsd','-ascii','-append');

disp('No. nodes');disp(kmap);
disp('No. displacement parameters');disp(kmap*(NPV-2));
disp('No. alive displacement parameters');disp(kmap*(NPV-2)-idead);
disp('No. dead displacement parameters');disp(idead);
disp('Updated MAP file');

return;
