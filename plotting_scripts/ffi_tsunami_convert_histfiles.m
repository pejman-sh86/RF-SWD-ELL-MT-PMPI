function []=ffi_tsunami_convert_histfiles(filebase);

parfile  = strcat(filebase,'_parameter.dat');
datfile  = strcat(filebase,'.hdf5');

NSTN  = h5readatt(datfile,'/Observed_data','N_sta');
deltt = h5readatt(datfile,'/Observed_data','Sample_rate');
NRAN  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx   = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx   = h5readatt(datfile,'/Sensitivity_kernel','max_y');
NTW   = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
%NTW   = 1;
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

%infilesl1G= strcat(filebase,'_ensavesl1G.dat');
%infilesl1G2= strcat(filebase,'_ensavesl1G2.dat');
infileVr  = strcat(filebase,'_ensavedelay.dat');
histfileVr  = strcat(filebase,'_histdelay.dat');

fid = fopen('plotpar.txt');
%% Discard first 8 entries, then read NHIST
for i=1:9;
  s = fgets(fid);
end;
NHIST = sscanf(s, '%d'); % convert to number

sl1_ens = zeros(NDEP,NRAN,NTW);
for itw=1:NTW;
  infilesl1 = strcat(filebase,'_ensavesl1_',num2str(itw,'%02i'),'.dat');
  sl1_ens(:,:,itw) = load(infilesl1);
end;
sl1_tot_ens = zeros(NDEP,NRAN);
infilesl1_tot = strcat(filebase,'_ensavesl1_tot.dat');
sl1_tot_ens = load(infilesl1_tot);
%sl1G_ens = load(infilesl1G);
%sl1G2_ens = load(infilesl1G2);
if(NPV > 2+NTW);Vr_ens = load(infileVr);end;

disp('loading Slip 1');
sl1_hst = zeros(NDEP*NHIST,NRAN,NTW);
for itw=1:NTW;
  histfilesl1 = strcat(filebase,'_histsl1_',num2str(itw,'%02i'),'.dat');
  sl1_hst(:,:,itw) = load(histfilesl1);
end;
sl1_tot_hst = zeros(NDEP*NHIST,NRAN);
histfilesl1_tot = strcat(filebase,'_histsl1_tot.dat');
sl1_tot_hst = load(histfilesl1_tot);
if(NPV > 2+NTW);
  Vr_hst = zeros(NDEP*NHIST,NRAN);
  disp('loading Rup. slowness');
  Vr_hst = load(histfileVr);
end;
bins  = load('bins.dat');

save tmp

sl1_hst2 = zeros(NDEP,NRAN,NHIST,NTW);
sl1_tot_hst2 = zeros(NDEP,NRAN,NHIST);
Vr_hst2 = zeros(NDEP,NRAN,NHIST);
for i=1:NHIST;
  for itw=1:NTW;
    sl1_hst2(1:NDEP,1:NRAN,i,itw) = sl1_hst((i-1)*NDEP+1:i*NDEP,:,itw);
  end;
  sl1_tot_hst2(1:NDEP,1:NRAN,i) = sl1_tot_hst((i-1)*NDEP+1:i*NDEP,:);
  if(NPV > 2+NTW);
    Vr_hst2(1:NDEP,1:NRAN,i) = Vr_hst((i-1)*NDEP+1:i*NDEP,:);
  end;
end;

clear sl1_hst sl1_tot_hst Vr_hst i;
save hists.mat;

return;
