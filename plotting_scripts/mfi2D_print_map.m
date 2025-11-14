function []=mfi2d_print_map(kmx,ISED,IWATER);

%ISED = 1;
%IWATER = 0;

if(ISED == 1);
  NVMX  = 20;
  NPV   = 5;
  NFPMX = NVMX*NPV;
  load rjmh_pecan_sample.mat
  AS = A;
  clear A;
end;
if(IWATER == 1);
  NVMXW = 20;
  NPVW  = 3;
  NFPMXW= NVMXW*NPVW;
  load rjmh_pecan_samplew.mat
  AW = A;
  clear A;
end;

if(ISED == 1)
  for i=min(AS(:,4)):max(AS(:,4));
    idx=find(AS(:,4)==i);
    [logLmx,ilogLmx]=max(AS(idx,1));
    if(IWATER == 1);
      disp([i,AW(idx(ilogLmx),4),logLmx,idx(ilogLmx),length(idx)]);
    else;
      disp([i,logLmx,ilogLmx,length(idx)]);
    end;
    clear idx;
  end;
else
  clear idx;
  for i=min(AW(:,4)):max(AW(:,4));
    idx=find(AW(:,4)==i);
    [logLmx,ilogLmx]=max(AW(idx,1));
    disp([i,AW(idx(ilogLmx),4),logLmx,idx(ilogLmx),length(idx)]);
    clear idx;
  end;
end;

%% Write sediment map:
if(ISED == 1);
  idx=find(AS(:,4)==kmx);
  [xmap,imap]=max(AS(idx,1))
  map = AS(idx(imap),5:NFPMX+4);
  kmap = AS(idx(imap),4);
  for ivo = 1:kmap;
    mapvoro(ivo,1:5) = map((ivo-1)*NPV+1:ivo*NPV);
  end;
  save('rjmh_pecan_map.dat','kmap','-ascii')
  save('rjmh_pecan_map.dat','mapvoro','-ascii','-append')
end;

%% Write water map:
if(IWATER == 1);
  if(ISED==0);
    idx=find(AW(:,4)==kmx);
    [xmap,imap]=max(AW(idx,1))
  end;
  mapw = AW(idx(imap),5:NFPMXW+4);
  kmapw = AW(idx(imap),4);
  for ivo = 1:kmapw;
    mapvorow(ivo,1:3) = mapw((ivo-1)*NPVW+1:ivo*NPVW);
  end;
  save('rjmh_pecan_mapw.dat','kmapw','-ascii')
  save('rjmh_pecan_mapw.dat','mapvorow','-ascii','-append')
end;
disp('Updated MAP file');

return;
