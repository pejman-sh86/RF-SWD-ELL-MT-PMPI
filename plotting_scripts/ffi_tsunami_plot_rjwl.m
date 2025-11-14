function []=ffi_tsunami_plot_rjwl(filename);

isave   = 1;
ifull   = 1;
ieps    = 1;
idatfit = 1;
isyn    = 0;
igauss  = 0;
inonstat= 1;

gridfact = 10.;

if(isyn == 1);
  %%
  %%  For simulations:
  %%
  %% double prior trans-D test
  gauge = [801,802,803,806,804,807,21418,21401,21413,5741,5742,5861,5862,1002,1006];
  time1 = [200, 200, 200, 200, 200, 200,1380,3000,3900, 600, 900, 150, 150, 150, 150];
  time2 = [3600,3600,3600,4600,3600,3600,3600,6800,6700,3600,3600,1200,2000,2000,2000];
  ylimmax = [ 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
  ylimmin = [-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];

  %% trans-D test simulation
  %gauge = [ 205, 801, 802, 803, 806, 804, 807,21418,21401,21413,5741,5742,5861,5862,1002,1006];
  %time1 = [ 100, 300, 300, 300, 300, 300, 300, 1380, 3000, 3900, 600, 900, 300, 300, 250, 250];
  %time2 = [3800,3600,2400,3600,4600,3600,3600, 3600, 6800, 6700,3600,3600,1200,2000,2000,2000];
  %ylimmax = [ 8, 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
  %ylimmin = [-2,-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];

else;
  %%
  %%  For Tohoku updated with new station (205)
  %%
  gauge = [ 205, 801, 802, 803, 806, 804, 807,21418,21401,21413,5741,5742,5861,5862,1002,1006];
  %% For original manuscript submission
  %time1 = [ 100, 300, 300, 300, 300, 300, 300, 1380, 3000, 3900, 600, 900, 300, 300, 250, 250];
  %time2 = [3800,3600,2400,3600,4600,3600,3600, 3600, 6800, 6700,3600,3600,1200,2000,2000,2000];
  %% For COV3 run (removed later parts... avoid non-linearity?)
  time1 = [ 100, 300, 300, 300, 300, 300, 300, 1380, 3000, 3900, 600, 900, 300, 300, 250, 250];
  time2 = [3800,3000,2400,2400,3600,2800,3200, 3200, 5600, 6400,3600,3600,1200,2000,2000,2000];  
  ylimmax = [ 8, 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
  ylimmin = [-2,-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];
end;
dtime = ceil((time2-time1)/400)*100;
damp = round((ylimmax-ylimmin)/4);

%% For prediction of unused stations:
%gauge = [202,203,205,602,613];
%time1 = [ 900, 900,  10,1500,1800];
%time2 = [7200,7200,7200,7200,7200];
%ylimmax = [ 3, 4.5, 7, 2.0, 1.5];
%ylimmin = [-2,-2.0,-4,-1.0,-1.0];

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');
imagefile_sl1= strcat(filebase,'_sl1_image.mat');
imagefile_sl2= strcat(filebase,'_sl2_image.mat');
imagefile_sl3= strcat(filebase,'_sl3_image.mat');
imagefile_sr = strcat(filebase,'_sr_image.mat');

imagefile_sl1_txt= strrep(imagefile_sl1,'.mat','.txt');
imagefile_sl2_txt= strrep(imagefile_sl2,'.mat','.txt');
imagefile_sl3_txt= strrep(imagefile_sl3,'.mat','.txt');
imagefile_sr_txt = strrep(imagefile_sr,'.mat','.txt');

[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,...
 NBUF,MAXDISP,MINDISP]=ffi_tsunami_wl_read_parfile(parfile);

MINDISP = -3;
MAXDISP =  8.5;
NSTN    = h5readatt(datfile,'/Observed_data','N_sta');
deltt   = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc = h5read(datfile,'/Observed_data/Hypo_Loc_Cart');
mat_excl= h5read(datfile,'/Sensitivity_kernel/mat_excl');
mat_excl = fliplr(mat_excl);
dhyp    = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx     = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx     = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin   = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax   = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
NTW     = h5readatt(datfile,'/Sensitivity_kernel','num_tw');

NTW = cast(NTW,'like',1);
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
NSTN = cast(NSTN,'like',1);
NTW

Vrmean = Vrmin+(Vrmax-Vrmin)/2.;

NRAN2 = NRAN * gridfact;
NDEP2 = NDEP * gridfact;

%% Prior bounds:
minlim(1:2) = [ 0.,  0.];
maxlim(1:2) = [Rmx, Zmx];

%gaugelatlon(1,1:2) = [];

if(filebase(1:5) == 'maule')
  %% Maule
  idxdat = [11:6:NSTN];
  xplt = [8 15 25 35 40];
  dvslip = 5;
elseif(filebase(1:5) == 'tohok');
  %% Tohoku
  %idxdat = [floor(NSTN/5):floor(NSTN/5):NSTN]
  idxdat = [1:NSTN];
  %idxdat = [   1, 4, 7,    12,14,  16];
  %ymin   = [-3.1,-6,-3.5,-1.5,-1,-0.4];
  %ymax   = [ 2.6, 6, 6.0, 4.2, 2, 1.0];
  if(NRAN <= 8)
    xplt = [2 4 6];
  else
    xplt = [10 11 12];
  end
  dvslip = 10;
elseif(isyn >= 1);
  %% Sim:
  idxdat = [11:6:NSTN];
  xplt = [8 15 25 35 40];
  if(isyn == 2);xplt = [13 14];end;
  dvslip = 5;
end;

minlim(NPV) = Vrmin;
maxlim(NPV) = Vrmax;
maxpert = maxlim-minlim;

%% Sub fault size in range and depth
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);

deltrn = double(deltr)/double(Rmx);
deltzn = double(deltz)/double(Zmx);

% Input files
mapfile        = strcat(filebase,'_map.dat');
covfile        = strcat(filebase,'_covmat.mat');
sdfile         = strcat(filebase,'_post_replicastddev.dat');
reparfile      = strcat(filebase,'_post_replicaar.dat');
repfile        = strcat(filebase,'_post_replica.dat');
vrmxfile       = strcat(filebase,'_t_Vrmx.txt');
ffdelfile      = strcat(filebase,'_tdelay.txt');

% Output files
plotfilek      = strcat(filebase,'_khist.');
plotfilelogL   = strcat(filebase,'_logL.');
plotfilemap1   = strcat(filebase,'_map1.');
plotfilemap2   = strcat(filebase,'_map2.');
plotfilehyp    = strcat(filebase,'_hyp.');
plotfilerupvel = strcat(filebase,'_rupvel.');
plotfilemaptot = strcat(filebase,'_maptot.');
plotfileens_sl1   = strcat(filebase,'_ens_sl1.');
plotfileens_vr   = strcat(filebase,'_ens_vr.');
plotfileens_sl1gf= strcat(filebase,'_ens_sl1GF.');
plotfileens_vrgf   = strcat(filebase,'_ens_vrGF.');
plotfileenstot = strcat(filebase,'_enstot.');
plotfiletru    = strcat(filebase,'_tru.');
plotfiletrugf  = strcat(filebase,'_trugf.');
plotfileslice_sl1  = strcat(filebase,'_slice_sl1.');
plotfileslice_sr  = strcat(filebase,'_slice_sr.');
plotfile95ci      = strcat(filebase,'_95ci.');
plotfilenodes  = strcat(filebase,'_nodes.');
plotfileMw     = strcat(filebase,'_Mw.');
plotfiledata   = strcat(filebase,'_data.');
plotfiledatasel= strcat(filebase,'_datasel.');
plotfileaxxsel= strcat(filebase,'_autocovsel.');
plotfileressel= strcat(filebase,'_ressel.');
plotfileres    = strcat(filebase,'_residuals.');
plotfilerestot = strcat(filebase,'_totalres.');
plotfilereshist    = strcat(filebase,'_reshist.');
plotfilerestothist = strcat(filebase,'_totalreshist.');
plotfileacovr  = strcat(filebase,'_autocovraw.');
plotfileacovtot= strcat(filebase,'_autocovtot.');
plotfilearpred = strcat(filebase,'_arpred.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

%%
%% Load sample file
%%
convert_sample(imagefile_sl1_txt);
load(filename);AS = A;
load(imagefile_sl1);Aimg_sl1 = A;
if(NTW > 1);
  convert_sample(imagefile_sl2_txt);
  load(imagefile_sl2);Aimg_sl2 = A;
end;
if(NTW > 2);
  convert_sample(imagefile_sl3_txt);
  load(imagefile_sl3);Aimg_sl3 = A;
end;
if(IVRUP == 1);
  convert_sample(imagefile_sr_txt);
  load(imagefile_sr);Aimg_sr = 10.^A;
end;
clear A;
NSMP  = size(Aimg_sl1,1);

%%
%% Load true model
%%
sltru2 = zeros(NDEP,NRAN,NTW);
if(isyn >= 1);
  sltru = load('true_sl1.txt');
  if(IVRUP==1);srtru = load('true_slowr.txt');end;
  for itw = 1:NTW;
    sltru2(:,:,itw) = sltru((itw-1)*NDEP+1:itw*NDEP,:);
  end;
  sltru_tot = sum(sltru2,3);
end;
if(idatfit == 1);
  dat = h5read(datfile,'/Observed_data/displacements');
  replin = h5read(datfile,'/Sensitivity_kernel/synthetic_displacements');
  %dat = dat';
  %replin = replin';
  disp('size of dat');
  disp(size(dat))
  disp('size of replin')
  disp(size(replin))
  
  if(isyn == 1);
    dat = dat';
    replin = replin';
  end;
  NTSMP = cast(h5read(datfile,'/Observed_data/Ntraces'),'like',1);
  NDAT = length(dat);
  rep=dlmread(repfile);
  NDSMP=size(rep,1);
  if(IAR == 1);repar=dlmread(reparfile);end;
  %% For rupture contours:
  if(IVRUP == 1);
    ffsrmn=load(vrmxfile);
    ffdel=load(ffdelfile);
  end;
end;
%%
%% Compute residual errors:
%%
if(idatfit == 1);
  res = dat(1,:)-rep(end,:);
end;
if(ICOV == 2 | ICOV == 4);
    load(covfile);
end;

for istn = 1:NSTN;
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  res(istart:iend) = res(istart:iend)-mean(res(istart:iend));
  resraw(istart:iend) = res(istart:iend)/std(res(istart:iend));
  if(ICOV == 2);
    L3 = chol(F(istn).Cd);
    resstd(istart:iend) = inv(L3')*res(istart:iend)';
    %if(ICOV == 3);
    %  resstd(istart:iend) = resstd(istart:iend)/sqrt(xi(istn));
    %end;
  else
    if(inonstat == 0);
      resstd(istart:iend) = res(istart:iend)./sd(istn,1:NTSMP(istn));
    else;
      resstd(istart:iend) = res(istart:iend)/std(res(istart:iend));
    end;
  end;
end; 
if(IAR == 1);
  resstd = dat(1,:)-rep(end,:)-repar(3,:);
  for istn = 1:NSTN;
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    resstd(istart:iend) = resstd(istart:iend)-mean(resstd(istart:iend));
    resstd(istart:iend) = resstd(istart:iend)/std(resstd(istart:iend));
  end; 
end;
for istn = 1:NSTN;
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  resstd(istart:iend) = resstd(istart:iend)-mean(resstd(istart:iend));
  resstd(istart:iend) = resstd(istart:iend)/std(resstd(istart:iend));
  std(resstd(istart:iend));
end; 

x_evn = [deltrn/2.:deltrn:1.-deltrn/2.];
z_evn = [deltzn/2.:deltzn:1.-deltzn/2.];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';

minlimmisc = [ hyp_loc(1)-dhyp, hyp_loc(2)-dhyp, Vrmin];
maxlimmisc = [ hyp_loc(1)+dhyp, hyp_loc(2)+dhyp, Vrmax];

minlimar = -0.5;
maxlimar =  1.0;

dx_ev = x_ev(2)-x_ev(1);

thinstep = 1;
logLmin = min(AS(:,1))-(max(AS(:,1))-min(AS(:,1)))/10;
logLmax = max(AS(:,1))+(max(AS(:,1))-min(AS(:,1)))/10;

k_sl(:,1) = AS(:,2);
if(NTW>1);k_sl(:,2) = AS(:,3);end;
if(NTW>2);k_sl(:,3) = AS(:,4);end;
if(IVRUP == 1);k_sr = AS(:,5);end;
logL = AS(:,1);
[maxL,imap] = max(logL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute DIC
%%
%ii = 1;
%for i=min(AS(:,4)):max(AS(:,4));
%  idx=find(AS(:,4)==i);
%  if(idx);
%  [logLmx(ii),ilogLmx(ii)]=max(AS(idx,1));
%  numk(ii) = length(idx);
%  logLmxidx(ii) = idx(ilogLmx(ii));
%  %disp([i,logLmx(ii),ilogLmx(ii),numk(ii)]);
%  ii = ii + 1;
%  end;
%  clear idx;
%end;
%D = -2.*logL;
%[tmp,idxkhat] = max(numk);
%Dmhat = D(logLmxidx(idxkhat));
%Dbar = mean(D);
%Pd = Dbar - Dmhat;
%DIC = Dmhat + 2 * Pd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figk_sl1=figure;
for itw=1:NTW;
  subplot(2,3,itw);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1);
  stairs(k_sl(:,itw),'k');
  ylabel('No. wavelet coefficients');
  xlabel('rjMCMC step');
  set(gca,'XLim',[0 size(AS,1)]);

  subplot(2,3,itw+3);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  [n,lim]=hist(k_sl(:,itw),[1:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim-0.5;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;
  xlabel('No. wavelet coefficients');
  ylabel('Probability');
  set(gca,'XLim',[0.5 NVMX+0.5]);
  set(gca,'YLim',[0.   1.0]);
end;
%%
%% K PLOT
%%
if(IVRUP==1);
  figk_sr=figure;
  subplot(2,1,1);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs(k_sr,'k')
  ylabel('No. nodes V_rup');
  xlabel('rjMCMC step');
  set(gca,'XLim',[0 size(k_sr,1)])

  subplot(2,1,2);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  [n,lim]=hist(k_sr,[1:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim-0.5;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;
  xlabel('No. nodes V_rup');
  ylabel('Probability');
  set(gca,'XLim',[0.5 NVMX+0.5]);
  set(gca,'YLim',[0.   1.0]);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%
figlogL=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
plot(AS(:,1),'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(AS(:,1))])
if(logLmin == logLmax);
   set(gca,'YLim',[0 2])
else;
   set(gca,'YLim',[logLmin logLmax])
end;

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(AS(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
if(logLmin == logLmax);
   set(gca,'YLim',[0 2])
else;
   set(gca,'YLim',[logLmin logLmax])
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Hypocentre plot
%%
hyp = AS(:,6:5+NMISC);
maphyp = hyp(imap,:);

fighyp = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
loc = [0.1,  0.55,  0.1, 0.55;
       0.555, 0.555, 0.14, 0.14 ];
spw = 0.44;
sph = 0.405;
spw2 = 0.15;
sph2 = 0.15;
subplot('Position',[loc(1,1) loc(2,1) spw sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim1]=hist(hyp(:,1),[minlimmisc(1):(maxlimmisc(1)-minlimmisc(1))/200:maxlimmisc(1)]);n = [0, n, 0];lim1 = [lim1(1) lim1 lim1(end)];
n = n/sum(n);
lim1 = lim1-(lim1(3)-lim1(2))/2.;
[xx,yy]=stairs(lim1,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim1,n,'k');
clear n;
set(gca,'XLim',[minlimmisc(1) maxlimmisc(1)]);
set(gca,'YLim',[0 0.08],'XTickLabel',[],'TickDir','out');
subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim2]=hist(hyp(:,2),[minlimmisc(2):(maxlimmisc(2)-minlimmisc(2))/200:maxlimmisc(2)]);n = [0, n, 0];lim2 = [lim2(1) lim2 lim2(end)];
n = n/sum(n);
lim2 = lim2-(lim2(3)-lim2(2))/2.;
[xx,yy]=stairs(n,lim2,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim2,'k');
clear n;
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[0 0.08],'YLim',[minlimmisc(2) maxlimmisc(2)]);
set(gca,'YTickLabel',[]);

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[minlimmisc(1) maxlimmisc(1)],'YLim',[minlimmisc(2) maxlimmisc(2)]);
set(gca,'TickDir','out');
clear n lim1 lim2;

if(IVRUP == 0 & NMISC == 3);
  figrupvel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.5 3])
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  [n,lim1]=hist(1./hyp(:,3),[minlimmisc(NMISC):(maxlimmisc(NMISC)-minlimmisc(NMISC))/300:maxlimmisc(NMISC)]);
  n = [0, n, 0];lim1 = [lim1(1) lim1 lim1(end)];
  n = n/sum(n);
  lim1 = lim1-(lim1(3)-lim1(2))/2.;
  [xx,yy]=stairs(lim1,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim1,n,'k');
  xlabel('Rupture velocity (km/s)');
  ylabel('Probability (s/km)');
  clear n lim1;
  set(gca,'XLim',[minlimmisc(NMISC) maxlimmisc(NMISC)]);
  set(gca,'YLim',[0 0.15],'TickDir','out');
  [minlimmisc(NMISC),maxlimmisc(NMISC)]
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% AR PLOT
%%
if(IAR == 1);
  figar = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  alpha = AS(:,9:8+NSTN);
  for istn=1:NSTN;
    subplot('Position',[loc(1,istn) loc(2,istn) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    [n,lim]=hist(alpha(:,istn),30);n = [0, n, 0];lim = [lim(1) lim lim(end)];
    n = n/sum(n);
    lim = lim-(lim(3)-lim(2))/2.;
    [xx,yy]=stairs(lim,n,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n,'k');
    clear n lim;

    set(gca,'YLim',[0 .3],'XLim',[minlimar maxlimar]);
    xlabel('AR(1)');
    set(gca,'XTick',[-1:.2:1]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;

  figar2 = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  for istn=1:NSTN;
    subplot('Position',[loc(1,istn) loc(2,istn) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    idxar = ones(size(alpha(:,istn)));
    idx = find(alpha(:,istn) < -0.5);
    idxar(idx) = 0;
    lim = [-1:1:2];
    [n]=hist(idxar,lim);%n = [0, n, 0];lim = [lim(1) lim lim(end)];
    n = n/sum(n);
    lim = lim - (lim(2)-lim(1))/2;
    [xx,yy]=stairs(lim,n,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n,'k');
    clear n lim;
    if(istn==1);ylabel('Prob. AR(1) off/on');end;
    xlabel('');
    if(istn>1);set(gca,'YTickLabel',[]);end;
    xticklabel = ({'off','on'});
    set(gca,'XTick',[0,1],'XTickLabel',xticklabel);
    set(gca,'XLim',[-1 2])
    set(gca,'YLim',[0 1.1])
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Load previously computed displacement image
%%
sl_mean = zeros(NDEP,NRAN,NTW);
sl_ens = zeros(NDEP,NRAN,NSMP,NTW);
for ismp=1:NSMP
  for iz=1:NDEP
    sl_ens(iz,:,ismp,1) = Aimg_sl1(ismp,(iz-1)*NRAN+1:iz*NRAN);
  end;
  sl_ens(:,:,ismp,1) = sl_ens(:,:,ismp,1).*mat_excl;
end;
sl_mean(:,:,1) = mean(sl_ens(:,:,:,1),3);
sl_med(:,:,1) = median(sl_ens(:,:,:,1),3);
slbin=[MINDISP:(MAXDISP-MINDISP)/799:MAXDISP];
for ir=1:NRAN
  for iz=1:NDEP
    [n1sl(iz,ir,:,1),slout]=hist(squeeze(sl_ens(iz,ir,:,1)),slbin);
    n1sl(iz,ir,:,1)=squeeze(n1sl(iz,ir,:,1))/trapz(slout,squeeze(n1sl(iz,ir,:,1)));    
    sl_prcnt = prctile(squeeze(sl_ens(iz,ir,:,1)), [2.5,97.5]);
    sl_hpd(iz,ir,1) = sl_prcnt(2)-sl_prcnt(1);
  end;
end;
if(NTW > 1);
  for ismp=1:NSMP
    for iz=1:NDEP
      sl_ens(iz,:,ismp,2) = Aimg_sl2(ismp,(iz-1)*NRAN+1:iz*NRAN);
    end;
    sl_ens(:,:,ismp,2) = sl_ens(:,:,ismp,2).*mat_excl;
  end;
  sl_mean(:,:,2) = mean(sl_ens(:,:,:,2),3);
  sl_med(:,:,2) = median(sl_ens(:,:,:,2),3);
  slbin=[MINDISP:(MAXDISP-MINDISP)/799:MAXDISP];
  for ir=1:NRAN
    for iz=1:NDEP
      [n1sl(iz,ir,:,2),slout]=hist(squeeze(sl_ens(iz,ir,:,2)),slbin);
      n1sl(iz,ir,:,2)=squeeze(n1sl(iz,ir,:,2))/trapz(slout,squeeze(n1sl(iz,ir,:,2)));
    end;
  end;
end;
if(NTW > 2);
  for ismp=1:NSMP
    for iz=1:NDEP
      sl_ens(iz,:,ismp,3) = Aimg_sl3(ismp,(iz-1)*NRAN+1:iz*NRAN);
    end;
    sl_ens(:,:,ismp,3) = sl_ens(:,:,ismp,3).*mat_excl;
  end;
  sl_mean(:,:,3) = mean(sl_ens(:,:,:,3),3);
  sl_med(:,:,3) = median(sl_ens(:,:,:,3),3);
  slbin=[MINDISP:(MAXDISP-MINDISP)/799:MAXDISP];
  for ir=1:NRAN
    for iz=1:NDEP
      [n1sl(iz,ir,:,3),slout]=hist(squeeze(sl_ens(iz,ir,:,3)),slbin);
      n1sl(iz,ir,:,3)=squeeze(n1sl(iz,ir,:,3))/trapz(slout,squeeze(n1sl(iz,ir,:,3)));
    end;
  end;
end;
if(1==2);
NY1=1;NY2=8;NX1=1;NX2=8;
ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
if(NDEP > 8);
  NY1=9;NY2=16;NX1=1;NX2=8;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
end;
if(NDEP > 16);
  NY1=1;NY2=8;NX1=9;NX2=16;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
  NY1=9;NY2=16;NX1=9;NX2=16;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);

  NY1=17;NY2=24;NX1=1;NX2=8;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
  NY1=17;NY2=24;NX1=9;NX2=16;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);

  NY1=25;NY2=32;NX1=1;NX2=8;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
  NY1=25;NY2=32;NX1=9;NX2=16;
  ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);
end;
end;

if(IVRUP == 1);
  sr_mean = zeros(NDEP,NRAN,NTW);
  sr_ens = zeros(NDEP,NRAN,NSMP);
  for ismp=1:NSMP
    for iz=1:NDEP
      sr_ens(iz,:,ismp) = Aimg_sr(ismp,(iz-1)*NRAN+1:iz*NRAN);
    end;
    sr_ens(:,:,ismp) = sr_ens(:,:,ismp).*mat_excl;
  end;
  sr_mean(:,:,1) = mean(sr_ens,3);
  sr_med(:,:,1) = median(sr_ens,3);
  srbin=[.0:(3.5)/799:3.5];
  for ir=1:NRAN
    for iz=1:NDEP
      [n1sr(iz,ir,:),slout]=hist(squeeze(sr_ens(iz,ir,:)),slbin);
      n1sr(iz,ir,:)=squeeze(n1sr(iz,ir,:))/trapz(slout,squeeze(n1sr(iz,ir,:)));
    end;
  end;
  %save ensemble_sr.mat sr_ens;
end;

nx = NTW;
if(isyn==1);nx=nx*2;end;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
%%
%% Ensemble mean figure slip components
%%
fig_postmeansl1=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

for itw=1:NTW;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(x_ev,z_ev,sl_mean(:,:,itw));shading flat;
%  for iplt=1:length(xplt);
%    plot([x_ev(xplt(iplt)) x_ev(xplt(iplt))],[z_ev(1) z_ev(end)],'-k')
%    plot([x_ev(xplt(iplt)) x_ev(xplt(iplt))],[z_ev(1) z_ev(end)],'--w')
%  end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
  set(gca,'CLim',[MINDISP MAXDISP],'FontSize',14)
  colormap(darkb2r(MINDISP,MAXDISP));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  cbfreeze(cb);
  set(get(cb,'ylabel'),'String', 'Posterior mean (m)','FontSize',14);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;

if(IVRUP == 1);
  fig_postmeanVr=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  axis equal;
  imagesc(x_ev,z_ev,sr_mean);shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  !v=[floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
  !contour(x_ev,z_ev,ffsrmn-ffdel,v,'-k','LineWidth',2);
  !contour(x_ev,z_ev,ffsrmn-ffdel,v,'--w','LineWidth',2);
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  set(gca,'TickDir','out')
  polarmap;caxis([-1,1]*max(abs(caxis))+Vrmean);
  set(gca,'CLim',[0.5 3.5],'FontSize',14);
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'Mean rup. vel. (km/s)','FontSize',14);
  xlabel('Trench perpendicular (km)');ylabel('Trench parallel (km)');
  box on;
end;
%%
%% Ensemble mean figure with GF smoothing
%%
fig_postmeansl1_gf=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

for itw=1:NTW;
  if(igauss == 0);
    [dep,ran,sl_meangf(:,:,itw)]=ffi_tsunami_gf_cos_smoother(deltr,sl_mean(:,:,itw));
  else;
    [dep,ran,sl_meangf(:,:,itw)]=ffi_tsunami_gf_gauss_smoother(deltr,sl_mean(:,:,itw));
  end;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(ran,dep,sl_meangf(:,:,itw));shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
  set(gca,'CLim',[MINDISP,MAXDISP],'FontSize',14)
  colormap(darkb2r(MINDISP,MAXDISP));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  cbfreeze(cb);
  set(get(cb,'ylabel'),'String', 'Posterior mean (m)','FontSize',14);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;

if(IVRUP == 1);
  if(igauss == 0);
    [dep,ran,sr_meangf]=ffi_tsunami_gf_cos_smoother(deltr,sr_mean(:,:));
  else
    [dep,ran,sr_meangf]=ffi_tsunami_gf_gauss_smoother(deltr,sr_mean(:,:));
  end;
  fig_postmeanVr_gf=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(ran,dep,sr_meangf);shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  set(gca,'TickDir','out')
  polarmap;caxis([-1,1]*max(abs(caxis))+Vrmean);
  set(gca,'CLim',[0.5 3.5],'FontSize',14);
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'Posterior mean (m)','FontSize',14);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;
save ensemble_sl.mat sl_ens sl_hpd sl_meangf dep ran;
%tmp = imresize(sl_meangf,1/3);
tmp=sl_meangf;
N2 = size(tmp,1);
M2 = size(tmp,2);
jj=1;
for ir=M2:-1:1;
    tmp_hi((jj-1)*N2+1:jj*N2) = tmp(:,ir);
    jj = jj + 1;
end
save('displacement_hi.txt','tmp_hi','-ascii');

tmp=squeeze(sl_med);
N2 = size(tmp,1);
M2 = size(tmp,2);
jj=1;
for ir=M2:-1:1;
    tmp_lo((jj-1)*N2+1:jj*N2) = tmp(:,ir);
    jj = jj + 1;
end
save('displacement_lo.txt','tmp_lo','-ascii');

tmp=squeeze(sl_hpd);
N2 = size(tmp,1);
M2 = size(tmp,2);
jj=1;
for ir=M2:-1:1;
    tmp_hpd((jj-1)*N2+1:jj*N2) = tmp(:,ir);
    jj = jj + 1;
end
save('displacement_95ci.txt','tmp_hpd','-ascii');

%%
%% Uncertainty map (95% CIs)
%%
fig95ci=figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
%imagesc([1:16]*15,[1:32]*15,sl_hpd);
imagesc(x_ev,z_ev,sl_hpd);
colormap(jet);
cb=colorbar('northoutside');
set(get(cb,'ylabel'),'String', 'Displacement 95% CI (m)','FontSize',14);
axis('equal');
set(gca,'XTick',[0:50:500],'YTick',[0:50:500],'TickDir','out');
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
if(isyn == 1);set(gca,'CLim',[0 1.2]);end;
xlabel('Trench perpendicular (km)');
ylabel('Trench parallel (km)');


%%
%% True figure
%%
if(isyn >= 1);
fig_trusl=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
for itw=1:NTW;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(x_ev,z_ev,sltru2(:,:,itw));shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  set(gca,'TickDir','out')
  set(gca,'CLim',[MINDISP,MAXDISP],'FontSize',14)
  colormap(darkb2r(MINDISP,MAXDISP));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'True displacement (m)','FontSize',14);
  cbfreeze(cb);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;

fig_truslgf=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
for itw=1:NTW;
  if(igauss == 0);
    [dep,ran,sltru2gf(:,:,itw)]=ffi_tsunami_gf_cos_smoother(deltr,sltru2(:,:,itw));
  else;
    [dep,ran,sltru2gf(:,:,itw)]=ffi_tsunami_gf_gauss_smoother(deltr,sltru2(:,:,itw));
  end;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(ran,dep,sltru2gf(:,:,itw));shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  set(gca,'TickDir','out')
  set(gca,'CLim',[MINDISP,MAXDISP],'FontSize',14)
  colormap(darkb2r(MINDISP,MAXDISP));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'True displacement (m)','FontSize',14);
  cbfreeze(cb);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;

if(IVRUP == 1);
  fig_truslow=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  axis equal;
  imagesc(x_ev,z_ev,srtru);shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  set(gca,'TickDir','out')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  polarmap;caxis([-1,1]*max(abs(caxis))+Vrmean);
  set(gca,'CLim',[0.5 3.5],'FontSize',14);
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'True rup. vel. (km/s)','FontSize',14);
  xlabel('Trench perpendicular (km)');ylabel('Trench parallel (km)');
  box on;
end;
end;
%% Total displacement
fig_postmeansltot=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
axis equal;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(ran,dep,sum(sl_meangf,3));shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
%plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'TickDir','out')
%colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
set(gca,'CLim',[MINDISP,MAXDISP],'FontSize',14)
colormap(darkb2r(MINDISP,MAXDISP));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
cbfreeze(cb);
set(get(cb,'ylabel'),'String', 'Posterior mean (m)','FontSize',14);
ylabel('Trench parallel (km)');
xlabel('Trench perpendicular (km)');
box on;
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
freezeColors;
  
if(isyn == 1);
  subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(ran,dep,sum(sltru2gf,3));shading flat;
%  for iplt=1:length(xplt);
%    plot([x_ev(xplt(iplt)) x_ev(xplt(iplt))],[z_ev(1) z_ev(end)],'-k')
%    plot([x_ev(xplt(iplt)) x_ev(xplt(iplt))],[z_ev(1) z_ev(end)],'--w')
%  end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
  set(gca,'CLim',[MINDISP,MAXDISP],'FontSize',14)
  colormap(darkb2r(MINDISP,MAXDISP));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  cbfreeze(cb);
  set(get(cb,'ylabel'),'String', 'True displacement (m)','FontSize',14);
  ylabel('Trench parallel (km)');
  xlabel('Trench perpendicular (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  xlim=get(gca,'XLim');set(gca,'XTickLabel',xlim(2)-get(gca,'XTick'));
  freezeColors;
end;
%%
%% Uncertainty
%%
nx = NRAN/2;
ny = NTW;
xim = 0.015;
yim = 0.015;
xymarg = [0.08 0.015 0.015 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

figslice_sl1=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
jj = 0;
for itw=1:NTW;
for iran=2:2:NRAN;
  jj = jj + 1;
  subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
  imagesc(slout,z_ev,squeeze(n1sl(:,iran,:,itw)));hold on;
  if(isyn == 1);
    stairs([sltru2(1,iran,itw); squeeze(sltru2(:,iran,itw)); sltru2(end,iran,itw)],...
           [z_ev(1)-deltz;z_ev;z_ev(end)+deltz]+deltz/2,'--w');
  end;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  set(gca,'TickDir','out')
  if(iran == 1);
    ylabel('Trench parallel (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  if(itw == NTW);
    xlabel('Displacement (m)');
  end;
  box on;
  set(gca,'XLim',[-1.8 8],'YLim',[0 Zmx],'YDir','reverse')
end;
end;
if(IVRUP == 1);
  figslice_sr=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  jj = 0;
  for iran=2:2:NRAN;
    jj = jj + 1;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    imagesc(slout,z_ev,squeeze(n1sr(:,iran,:)));hold on;
    if(isyn == 1);
      stairs([srtru(1,iran); squeeze(srtru(:,iran)); srtru(end,iran)],...
             [z_ev(1)-deltz;z_ev;z_ev(end)+deltz]+deltz/2,'--w');
    end;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    set(gca,'TickDir','out')
    if(iran == 1);
      ylabel('Trench parallel (km)');
    else;
      set(gca,'YTickLabel',[]);
    end;
    if(itw == NTW);
      xlabel('Rup. vel. (km/s)');
    end;
    box on;
    set(gca,'XLim',[0.5 3.5],'YLim',[0 Zmx],'YDir','reverse')
  end;
end;

if(idatfit == 1);
  %%
  %% Selected Data plots:
  %%
  figdatsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 9])
  nx = 4;
  ny = 4;
  xim = 0.1/nx;
  yim = 0.1/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
  idx2 = [];
  deltt = double(deltt);
  
  
  for jdat=1:length(idxdat);
    idat = idxdat(jdat);
    h(jdat)=subaxis(ny,nx,jdat,'Spacing',0.01,'Padding',0,'Paddingbottom',0.032,...
        'Paddingleft',0.025,'ML', 0.04,'MR',0.01,'MB',.04,'MT',.012);hold on;box on;
    %subplot('Position',[loc_s(1,jdat) loc_s(2,jdat) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1,'TickDir','out')
    iend = sum(NTSMP(1:idat));
    istart = iend-NTSMP(idat)+1;
    idx2 = [idx2,[istart:iend]];
    sigma = std(dat(1,istart:iend)-rep(end,istart:iend));
    %disp([istart,iend,deltt]);
    plot(deltt*[0:NTSMP(idat)-1]+time1(idat),replin(1,istart:iend),'-k','LineWidth',1.5);
    plot(deltt*[0:NTSMP(idat)-1]+time1(idat),replin(1,istart:iend),'--w','LineWidth',1.5);
    p1=plot(deltt*[0:NTSMP(idat)-1]+time1(idat),dat(1,istart:iend),'b','LineWidth',1.5);
    for isub=1:10:NDSMP;
      p3=plot(deltt*[0:NTSMP(idat)-1]+time1(idat),rep(isub,istart:iend),':r');
    end;
    if(jdat > nx*(ny-1));
      xlabel('Time (s)');
    end;
    set(gca,'XTick',[time1(idat):dtime(idat):time2(idat)]);
    %set(gca,'YTick',[-8:dtime(idat):8]);
    text(time1(idat)+(time2(idat)-time1(idat))/12,ylimmax(idat)-(ylimmax(idat)-ylimmin(idat))/12,num2str(gauge(idat)),'FontSize',16);
    if(jdat == 1 | jdat == nx+1 | jdat == 2*nx+1 | jdat == 3*nx+1);
      ylabel('Displacement (m)');
    end;
    set(gca,'YLim',[ylimmin(jdat),ylimmax(jdat)],'XLim',[time1(idat) time2(idat)]);
    box on;
  end;
  if(ICOV == 2 | ICOV == 4);
  %%
  %%  TOTAL Autocovariance
  %%
  figaxxsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
  idx2 = [];
  for jdat=1:length(idxdat);
    idat = idxdat(jdat);
    subplot('Position',[loc_s(1,jdat) loc_s(2,jdat) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:idat));
    istart = iend-NTSMP(idat)+1;

    axx = xcorr(resraw(istart:iend),'coeff');
    plot(deltt*[-NTSMP(idat)+1:NTSMP(idat)-1],axx,'k','LineWidth',1.5,'color', [0.6 0.6 0.6]);
    clear axx;

    axxtot = xcorr(resstd(istart:iend),'coeff');
    plot(deltt*[-NTSMP(idat)+1:NTSMP(idat)-1],axxtot,'k','LineWidth',1.5);
    plot(deltt*[-max(NTSMP(:)),max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axxtot;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    set(gca,'XTick',[-1600:400:1600]);
    if(jdat > 12);
      xlabel('lag (s)');
    else
      set(gca,'XTickLabel',[]);
    end;
    if(jdat == 1 | jdat == 5 | jdat == 9 | jdat == 13);
      ylabel('');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  end;
  
  if(ICOV == 2| IAR == 1);
    %%
    %%  TOTAL RESIDUAL HISTOGRAMS
    %%
    x = -6.2:.4:6.2;
    xx = -6.25:.01:6.25;
    figressel = figure();hold on;box on;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
    idx2 = [];
    for jdat=1:length(idxdat);
      idat = idxdat(jdat);
      subplot('Position',[loc_s(1,jdat) loc_s(2,jdat) spw_s sph_s]);hold on;box on;
      set(gca,'FontSize',14,'layer','top','LineWidth',1)
      iend = sum(NTSMP(1:idat));
      istart = iend-NTSMP(idat)+1;

      [n1,xout] = hist(resstd(istart:iend),x);
      n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
      area = sum(n1) * (xout(3)-xout(2));
      n1 = n1/area;
      [xx,yy]=stairs(xout,n1);
      patch(xx,yy,[0.7,0.7,0.7]);
      stairs(xout,n1,'k','LineWidth',1.5);
      clear n1 xout;

      [n1,xout] = hist(resraw(istart:iend),x);
      n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
      area = sum(n1) * (xout(3)-xout(2));
      n1 = n1/area;
      stairs(xout,n1,'r','LineWidth',1);
      clear n1 xout;

      nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
      plot(xx,nd,'-b','LineWidth',1)

      set(gca,'XLim',[-4.5 4.5]);
      %if(jdat == 3 | jdat == 4);
      if(jdat > 12);
        xlabel('Residuals (std. dev.)');
      else
        set(gca,'XTickLabel',[]);
      end;
      ylabel('');set(gca,'YTickLabel',[]);
      box on;
    end;
  end;
end;
if(isave == 1)
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext2),'-dpng');
  saveas(figlogL,strcat(plotfilelogL,plotext1),'fig');
  print(figk_sl1,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
  saveas(figk_sl1,strcat(plotfilek,plotext1),'fig');
  %print(figmap1,'-r250',strcat(plotfilemap1,plotext2),'-dpng');
  %saveas(figmap1,strcat(plotfilemap1,plotext1),'fig');
  %print(figmap2,'-r250',strcat(plotfilemap2,plotext2),'-dpng');
  %saveas(figmap2,strcat(plotfilemap2,plotext1),'fig');
  print(fig_postmeansl1,'-painters','-r250',strcat(plotfileens_sl1,plotext2),'-dpng');
  saveas(fig_postmeansl1,strcat(plotfileens_sl1,plotext1),'fig');
  print(fig_postmeansl1_gf,'-painters','-r250',strcat(plotfileens_sl1gf,plotext2),'-dpng');
  saveas(fig_postmeansl1_gf,strcat(plotfileens_sl1gf,plotext1),'fig');
  print(fig_postmeansltot,'-painters','-r250',strcat(plotfileenstot,plotext2),'-dpng');
  saveas(fig_postmeansltot,strcat(plotfileenstot,plotext1),'fig');
  if(IVRUP == 1);
    print(fig_postmeanVr,'-painters','-r250',strcat(plotfileens_vr,plotext2),'-dpng');
    saveas(fig_postmeanVr,strcat(plotfileens_vr,plotext1),'fig');
    print(fig_postmeanVr_gf,'-painters','-r250',strcat(plotfileens_vrgf,plotext2),'-dpng');
    saveas(fig_postmeanVr_gf,strcat(plotfileens_vrgf,plotext1),'fig');    
  end;
  print(fig95ci,'-painters','-r250',strcat(plotfile95ci,plotext2),'-dpng');
  saveas(fig95ci,strcat(plotfile95ci,plotext1),'fig');
  %print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext2),'-dpng');
  %saveas(fighpd,strcat(plotfilehpd,plotext1),'fig');
  if(isyn >= 1);
    print(fig_trusl,'-painters','-r250',strcat(plotfiletru,plotext2),'-dpng');
    saveas(fig_trusl,strcat(plotfiletru,plotext1),'fig');
    print(fig_truslgf,'-painters','-r250',strcat(plotfiletrugf,plotext2),'-dpng');
    saveas(fig_truslgf,strcat(plotfiletrugf,plotext1),'fig');
  end;
  print(figslice_sl1,'-painters','-r250',strcat(plotfileslice_sl1,plotext2),'-dpng');
  saveas(figslice_sl1,strcat(plotfileslice_sl1,plotext1),'fig');
  if(IVRUP == 1);
    print(figslice_sr,'-painters','-r250',strcat(plotfileslice_sr,plotext2),'-dpng');
    saveas(figslice_sr,strcat(plotfileslice_sr,plotext1),'fig');
  end;
  print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext2),'-dpng');
  saveas(figdatsel,strcat(plotfiledatasel,plotext1),'fig');
  if(ICOV == 2 | ICOV == 4);
      print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext2),'-dpng');
      saveas(figaxxsel,strcat(plotfileaxxsel,plotext1),'fig');
      print(figressel,'-painters','-r250',strcat(plotfileres,plotext2),'-dpng');
      saveas(figressel,strcat(plotfileres,plotext1),'fig');
  end;
  if(ieps == 1);
    print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-depsc');
    print(figk_sl1,'-painters','-r250',strcat(plotfilek,plotext3),'-depsc');
    print(fig_postmeansl1,'-painters','-r250',strcat(plotfileens_sl1,plotext3),'-depsc');
    print(fig_postmeansl1_gf,'-painters','-r250',strcat(plotfileens_sl1gf,plotext3),'-depsc');
    print(fig_postmeansltot,'-painters','-r250',strcat(plotfileenstot,plotext3),'-depsc');
    print(fig95ci,'-painters','-r250',strcat(plotfile95ci,plotext3),'-depsc');
    if(IVRUP == 1);
      print(fig_postmeanVr,'-painters','-r250',strcat(plotfileens_vr,plotext3),'-depsc');
      print(fig_postmeanVr_gf,'-painters','-r250',strcat(plotfileens_vrgf,plotext3),'-depsc');
    end;
    if(isyn >= 1);
      print(fig_trusl,'-painters','-r250',strcat(plotfiletru,plotext3),'-depsc');
      print(fig_truslgf,'-painters','-r250',strcat(plotfiletrugf,plotext3),'-depsc');
    end;
    print(figslice_sl1,'-painters','-r250',strcat(plotfileslice_sl1,plotext3),'-depsc');
    if(IVRUP == 1);
      print(figslice_sr,'-painters','-r250',strcat(plotfileslice_sr,plotext3),'-depsc');
    end;
    print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext3),'-depsc'); 
    if(ICOV == 2 | ICOV == 4);
        print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext3),'-depsc');
        print(figressel,'-painters','-r250',strcat(plotfileres,plotext3),'-depsc');
    end;
  end;
end

return;
