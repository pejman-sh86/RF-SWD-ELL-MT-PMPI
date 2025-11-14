function []=cmt_print_map(filebase,kmx);

parfile = strcat(filebase,'_parameter.dat');
infile = strcat(filebase,'_sample.mat');
outfile = strcat(filebase,'_map.dat');
cmtfile = strcat(filebase,'_cmt.eqs');
cmtlinfile = strcat(filebase,'_cmtlin.eqs');

[IMAP,ICOV,NVMX,ILOC,IAR,IEXCHANGE,IDBLCPL,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay,minMw,minstr,mindip,minrak,...
 maxMw,maxstr,maxdip,maxrak]=cmt_read_parfile(parfile);

[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

isyn = 0;
NCMT = 5;
if(IDBLCPL == 1);
  NCMT2 = 4;
else
  NCMT2 = 5;
end;
NLOC = 4;
NPV = NCMT2+NLOC;

%%
%% All thses colutions are Mrr, Mtt, Mff, Mrt, Mrf, Mtf
%% and scale is 10^22
%%
if(filebase(1:5) == 'maule')
  %% Maule linear solution:
  voro_lin = [1.07801073, -2.20406943e-03, -1.07580666,...
             -2.11090666e-03,-2.97344950, -1.50030506e-01,...
             -35.509, -72.3329, 15.5, 57.0];
elseif(filebase(1:5) == 'tohok');
  %% Tohoku linear solution:
  voro_lin = [1.73932071, -1.51540195e-01, -1.58778051,...
              2.21198538,  4.92560869,     -0.60684836,...
              37.52, 142.71, 11.5, 68.0];
elseif(filebase(1:4) == 'hg12');
  %% Tohoku linear solution:
  voro_lin = [4.016515e-02, -2.520619e-02, -1.495896e-02, ...
              4.012309e-02, -2.104478e-02, 2.230925e-02, ...
              52.54999999999999, -131.43999999999997, 17.5, 28.0];
elseif(filebase(1:4) == 'bo13');
  %% Tohoku linear solution:
  voro_lin = [4.016515e-02, -2.520619e-02, -1.495896e-02, ...
              4.012309e-02, -2.104478e-02, 2.230925e-02, ...
              52.54999999999999, -131.43999999999997, 17.5, 28.0];
elseif(filebase(1:4) == 'am16');
  %% Tohoku linear solution:
  voro_lin = [4.016515e-02, -2.520619e-02, -1.495896e-02, ...
              4.012309e-02, -2.104478e-02, 2.230925e-02, ...
              52.54999999999999, -131.43999999999997, 17.5, 28.0];
elseif(isyn >= 1);
  %% Sim:
end;

NFPMX = NVMX*NPV;
NFPMX
NSTN

load(infile);
AS = A;
clear A;

for i=min(AS(:,4)):max(AS(:,4));
  idx=find(AS(:,4)==i);
  [logLmx,ilogLmx]=max(AS(idx,1));
  disp([i,logLmx,ilogLmx,length(idx)]);
  clear idx;
end;

[xmap,imap]=max(AS(:,1))
map = AS(imap,5:end-6);
kmap = AS(imap,4);
for ivo = 1:kmap;
  mapvoro(ivo,1:NPV) = map((ivo-1)*NPV+1:ivo*NPV);
end;
mapar = map(NVMX*NPV+1:NVMX*NPV+NSTN);
mapsd = map(NVMX*NPV+NSTN+1:NVMX*NPV+NSTN+NSTN);
save(outfile,'kmap','-ascii','-double');
save(outfile,'mapvoro','-ascii','-append','-double');
save(outfile,'mapar','-ascii','-append','-double');
save(outfile,'mapsd','-ascii','-append','-double');

return;

cmt = [mapvoro(:,7),mapvoro(:,6),mapvoro(:,8),...
       mapvoro(:,1:2),-mapvoro(:,1)-mapvoro(:,2),...
       mapvoro(:,3:5)];
cmt_lin = [voro_lin(8),voro_lin(7),voro_lin(9),voro_lin(1:6)];

% 30 X Y 201002270635A

tmp = [cmt(1,:),22,0,0];
save(cmtfile,'tmp','-ascii');
%if(NVMX > 1);
  for ivo = 1:kmap;
    tmp = [cmt(ivo,:),22,0,0,mapvoro(ivo,NPV)];
    save(cmtfile,'tmp','-ascii','-append');
  end;
%end;

tmp = [cmt_lin,22,0,voro_lin(10)];
save(cmtlinfile,'tmp','-ascii');
save(cmtlinfile,'tmp','-ascii','-append');

disp('Updated MAP file');

return;
