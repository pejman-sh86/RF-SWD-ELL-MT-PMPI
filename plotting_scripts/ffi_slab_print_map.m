function []=ffi_mtw_print_map(filebase,kmx);

  parfile = strcat(filebase,'_parameter.dat');
  infile = strcat(filebase,'_sample.mat');
  datfile = strcat(filebase,'.hdf5');
  datfilecGPS = strcat(filebase,'_cGPS.hdf5');
  outfile = strcat(filebase,'_map.dat');
  [IMAP,ICOV,I_WP,I_cGPS,I_GPS,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
   NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile);

  load(infile);
  NSTN = h5readatt(datfile,'/Observed_data','N_sta');
  NTW  = h5readatt(datfile,'/Sensitivity_kernel','N_time_win');
  if(I_cGPS == 1);
    NcGPS = h5readatt(datfilecGPS,'/Observed_data','N_sta');
  end;
  AS = A;
  clear A;
  NFPMX = NVMX*NPV*NTW;

  [xmap,imap]=max(AS(:,1));
  ist = 4+NTW;
  map = AS(imap,ist:end);
  kmap = AS(imap,4:ist-1);
  ktmp = kmap(1);
  save(outfile,'ktmp','-ascii');
  ipar = 1;
  for itw = 1:NTW;
    for ivo = 1:kmap(itw);
      %disp([ivo,kmap(itw),ipar,ipar+NPV-1])
      x(itw).mapvoro(ivo,1:NPV) = map(ipar:ipar+NPV-1);
      ipar = ipar + NPV;
    end;
    mapvoro = x(itw).mapvoro;
    if(itw > 1);
      ktmp = kmap(itw);
      save(outfile,'ktmp','-ascii','-append');
    end
    save(outfile,'mapvoro','-ascii','-append');
  end;

  maphyp     = map(NFPMX+1:NFPMX+NMISC);
  if(ICOV == -1 | ICOV == 3);
      mapsd      = zeros(1,NSTN);
      mapsdGPS   = zeros(1,3);
      if(I_cGPS == 1);
        mapsdcGPS  = zeros(1,NcGPS);
      else
        mapsdcGPS  = zeros(1,1);
      end;
      mapsdcGPSs = zeros(1,3);
      mapalp     = map(NFPMX+NMISC+1);
  else;
      mapsd      = map(NFPMX+NMISC+1:NFPMX+NMISC+NSTN);
      mapsdGPS   = map(NFPMX+NMISC+NSTN+1:NFPMX+NMISC+NSTN+3);
      if(I_cGPS == 1);
        mapsdcGPS  = map(NFPMX+NMISC+NSTN+3+1:NFPMX+NMISC+NSTN+3+NcGPS);
      else
        mapsdcGPS  = zeros(1,1);
      end;
      mapsdcGPS  = map(NFPMX+NMISC+NSTN+3+1:NFPMX+NMISC+NSTN+3+NcGPS);
      mapsdcGPSs = map(NFPMX+NMISC+NSTN+3+NcGPS+1:NFPMX+NMISC+NSTN+3+NcGPS+3);
      mapalp     = map(NFPMX+NMISC+NSTN+NcGPS+3+3+1);
  end
  save(outfile,'maphyp','-ascii','-append');
  save(outfile,'mapsd','-ascii','-append');
  save(outfile,'mapsdcGPS','-ascii','-append');
  save(outfile,'mapsdcGPSs','-ascii','-append');
  save(outfile,'mapsdGPS','-ascii','-append');
  save(outfile,'mapalp','-ascii','-append');

disp('Updated MAP file');

save tmp
return;
