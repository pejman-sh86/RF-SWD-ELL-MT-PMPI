function []=ffi_tsunami_slip_slice(filename);

filebase = strrep(filename,'_sample.mat','');
datfile  = strcat(filebase,'.hdf5');
parfile  = strcat(filebase,'_parameter.dat');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,...
 NBUF,MAXDISP,MINDISP]=ffi_tsunami_wl_read_parfile(parfile);


load ensemble_sl.mat;
NENS = size(sl_ens,3);
%NENS = 1000;
NTW  = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
NRAN = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NTW  = cast(NTW,'like',1);
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
Rmx     = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx     = h5readatt(datfile,'/Sensitivity_kernel','max_y');

deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);

sl_ensgf = zeros(size(sl_meangf,1),size(sl_meangf,2));
sl_ensgf_sum = zeros(size(sl_meangf,1),size(sl_meangf,2));

for iens=1:NENS;
    if(rem(iens,1000)==0);iens,end;
    [dep,ran,sl_ensgf(:,:)]=ffi_tsunami_gf_smoother(deltr,sl_ens(:,:,iens));
    sl_ensgf_sum = sl_ensgf_sum+sl_ensgf;
end;
save tmp.mat dep ran sl_ensgf sl_ensgf_sum;
NRAN2 = length(ran);
NDEP2 = length(dep);
slbin=[MINDISP:(MAXDISP-MINDISP)/399:MAXDISP];
NBIN = length(slbin);
n1sl = zeros(size(sl_meangf,1),size(sl_meangf,2),NBIN);
for ir=1:NRAN2
  if(rem(ir,100)==0);ir,end;
  for iz=1:NDEP2
    [n1sl(iz,ir,:),slout]=hist(squeeze(sl_ensgf(iz,ir,:)),slbin);
    n1sl(iz,ir,:)=n1sl(iz,ir,:)/trapz(slout,n1sl(iz,ir,:));
  end;
end;
%save('tmp.mat','dep','ran','sl_ensgf','n1sl','slout','slbin','-v7.3');

skip = 64;
nx= floor(NRAN2/skip);
ny= 1;
ML = .045;
MR = .03;
MB = .07;
MT = .02;
SP = .02;
PAD = 0;
FNT = 12;
figslice_sl1=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8])
jj = 0;
for iran=2:skip:NRAN2;
  jj = jj + 1;
  subaxis(ny,nx,jj,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);  
  imagesc(slout,dep,squeeze(n1sl(:,iran,:)));hold on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  set(gca,'TickDir','out')
  if(iran == 2);
    ylabel('Northing (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  xlabel('Displacement (m)');
  box on;
  set(gca,'XLim',[-1.8 8],'YLim',[0 Zmx],'YDir','reverse')
end;

return;
