function []=ffi_convert_histfiles(filebase);
%filebase = 'bohol_data';

parfile  = strcat(filebase,'_parameter.dat');
datfile  = strcat(filebase,'.hdf5');

NSTN  = h5readatt(datfile,'/Observed_data','N_sta');
deltt = h5readatt(datfile,'/Observed_data','Sample_rate');
NRAN  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx   = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx   = h5readatt(datfile,'/Sensitivity_kernel','max_y');
[IMAP,ICOV,I_WP,I_GPS,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile)

infilesl1 = strcat(filebase,'_ensavesl1.dat');
infilesl2 = strcat(filebase,'_ensavesl2.dat');
infileslt = strcat(filebase,'_ensavet.dat');
infilealf = strcat(filebase,'_ensavealf.dat');
infileVr  = strcat(filebase,'_ensavedelay.dat');
histfilesl1 = strcat(filebase,'_histsl1.dat');
histfilesl2 = strcat(filebase,'_histsl2.dat');
histfilet = strcat(filebase,'_histst.dat');
histfilealf = strcat(filebase,'_histsalf.dat');
histfileVr  = strcat(filebase,'_histdelay.dat');

fid = fopen('plotpar.txt');
%% Discard first 8 entries, then read NHIST
for i=1:9;
  s = fgets(fid);
end;
NHIST = sscanf(s, '%d'); % convert to number

sl1_ens = load(infilesl1);
sl2_ens = load(infilesl2);
sltot_ens = load(infileslt);
alf_ens = load(infilealf);
if(NPV == 5);Vr_ens = load(infileVr);end;

disp('loading Slip 1');
sl1_hst = load(histfilesl1);
disp('loading Slip 2');
sl2_hst = load(histfilesl2);
disp('loading Total Slip');
sltot_hst = load(histfilet);
disp('loading rake angles');
alf_hst = load(histfilealf);
if(NPV == 5);
  disp('loading Rup. slowness');
  Vr_hst = load(histfileVr);
end;
bins  = load('bins.dat');

sl1_hst2 = zeros(NDEP,NRAN,NHIST);
NDEP
NRAN
NHIST
size(sl1_hst)

for i=1:NHIST;
  i
  sl1_hst2(1:NDEP,1:NRAN,i) = sl1_hst((i-1)*NDEP+1:i*NDEP,:);
  sl2_hst2(1:NDEP,1:NRAN,i) = sl2_hst((i-1)*NDEP+1:i*NDEP,:);
  sltot_hst2(1:NDEP,1:NRAN,i) = sltot_hst((i-1)*NDEP+1:i*NDEP,:);
  alf_hst2(1:NDEP,1:NRAN,i) = alf_hst((i-1)*NDEP+1:i*NDEP,:);
  if(NPV == 5);
    Vr_hst2(1:NDEP,1:NRAN,i) = Vr_hst((i-1)*NDEP+1:i*NDEP,:);
  end;
end;

clear sl1_hst sl2_hst sltot_hst alf_hst Vr_hst i;
save hists.mat;

return;
