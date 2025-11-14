function [] = ffi_plot_single_voro(filename);

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF]=ffi_read_parfile(parfile);

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc');
rake   = h5readatt(datfile,'/Sensitivity_kernel','rake_comp');
dhyp   = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin  = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax  = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
mu     = h5read(datfile,'/Rigidity/mu');

%% Prior bounds:
minlim(1:2) = [ 0.,  0.];
maxlim(1:2) = [Rmx, Zmx];

load(filename);
m = A(:,5:84);
k = A(:,4);
NSMP = length(k);

%% Sub fault size in range and depth
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
%% Newton meter scaling for moment
mu = repmat(mu,1,NRAN);
Nmscale = double(mu)*1.e6*deltr*deltz;

deltrn = double(deltr)/double(Rmx);
deltzn = double(deltz)/double(Zmx);

x_evn = [deltrn/2.:deltrn:1.-deltrn/2.];
z_evn = [deltzn/2.:deltzn:1.-deltrn/2.];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';

j = 1;
NSUB = 1;
NENS = NSMP/NSUB
sl1_ens = zeros(NDEP,NRAN,NENS);
sl2_ens = zeros(NDEP,NRAN,NENS);
for ismp = NSUB:NSUB:1%NSMP;
  if(rem(ismp,1000) == 0);disp(ismp);end;
  voro = zeros(k(ismp),NPV);
  for ivo = 1:k(ismp);
    voro(ivo,:) = m(ismp,(ivo-1)*NPV+1:ivo*NPV);
  end;
  voro
  x  = voro(:,1);
  z  = voro(:,2);
  xn = voro(:,1)/Rmx;
  zn = voro(:,2)/Zmx;
  %% Plot MAP
  for iz = 1:NDEP;
    dz = z_evn(iz)-zn;
    for ir = 1:NRAN;
       dx = x_evn(ir)-xn;
       d  = sqrt(dx.*dx + dz.*dz);
       [dmin,iv] = min(d);
       sl1_evmap(iz,ir) = voro(iv,3);
       sl2_evmap(iz,ir) = voro(iv,4);
       %Sr_evmap(iz,ir)  = voro(iv,5);
    end;
  end;
  sl1_ens(:,:,j) = sl1_evmap;
  sl2_ens(:,:,j) = sl2_evmap;
  j = j + 1;

%  imagesc(x_ev,z_ev,sl1_evmap);shading flat;
end;
figure();imagesc(x_ev,z_ev,mean(sl1_ens,3));shading flat;
figure();imagesc(x_ev,z_ev,mean(sl2_ens,3));shading flat;



return;
