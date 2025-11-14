function []=ffi_plot_rjvoro(filename);

isave   = 1;
ieps    = 0;
idatfit = 1;
isyn    = 0;
inonstat= 1;
ialf    = 0;

alfmin = 59;
NSDGPS = 4; %% 3 sd for x, y, z, and 1 for alp (Okada "constant")

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,I_WP,I_cGPS,I_GPS,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile)

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
sfgrid = h5read(datfile,'/Sensitivity_kernel/subfaultgrid');

NSF = NRAN*NDEP;

%% Prior bounds:
minlim(1:2) = [min(sfgrid(2,:)),min(sfgrid(1,:))];
maxlim(1:2) = [max(sfgrid(2,:)),max(sfgrid(1,:))];

%% Mesh to grid data on:
[xq,yq] = meshgrid(minlim(1):.05:maxlim(1), minlim(2):.05:maxlim(2));
xxq = xq(1,:);
yyq = yq(:,1);

if(filebase(1:5) == 'maule')
  %% Maule
  idxdat = [11:6:NSTN];
  %xplt = [8 15 25 35 40];
  xplt = [8 15 25 40];
  dvslip = 7.5;
  dvslip2 = 10.;
elseif(filebase(1:5) == 'tohok');
  %% Tohoku
  idxdat = [10:10:NSTN];
  xplt = [3 7 15 22];
  dvslip = 15;
  dvslip2 = 15.;
%% Sim:
elseif(isyn >= 1);
  idxdat = [11:6:NSTN];
  xplt = [8 15 22];
  if(isyn == 2);xplt = [13 14];end;
  dvslip = 5;
  dvslip2 = 5.;
end;

fid = fopen('plotpar.txt');
%% Discard first 8 entries, then read NHIST
for i=1:12;
  s = fgets(fid);
  if(i == 9);NHIST = sscanf(s, '%d');end; % convert to number
  if(i ==11);minlim(3:5) = sscanf(s, '%f');end; % convert to number
  if(i ==12);maxlim(3:5) = sscanf(s, '%f');end; % convert to number
end;
minlim(5) = Vrmin;
maxlim(5) = Vrmax;
maxpert = maxlim-minlim;

%% Sub fault size in range and depth
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
%% Newton meter scaling for moment
Nmscale = double(mu)*1.e6*deltr*deltz;
%deltrn = double(deltr)/double(Rmx);
%deltzn = double(deltz)/double(Zmx);

% Output files
mapfile        = strcat(filebase,'_map.dat');
covfile        = strcat(filebase,'_covmat.mat');
sdfile         = strcat(filebase,'_nonstat_sd.dat');
reparfile      = strcat(filebase,'_replicaar.dat');
repfile        = strcat(filebase,'_replica.dat');
repfileGPS     = strcat(filebase,'_replicaGPS.dat');
vrmxfile       = strcat(filebase,'_t_Vrmx.txt');
ffdelfile      = strcat(filebase,'_tdelay.txt');
plotfilek      = strcat(filebase,'_khist.');
plotfilelogL   = strcat(filebase,'_logL.');
plotfilemap    = strcat(filebase,'_map.');
plotfilehyp    = strcat(filebase,'_hyp.');
plotfilerupvel = strcat(filebase,'_rupvel.');
plotfilemaptot = strcat(filebase,'_maptot.');
plotfileens    = strcat(filebase,'_ens.');
plotfileenstot = strcat(filebase,'_enstot.');
plotfilerake   = strcat(filebase,'_rake.');
plotfilealf   = strcat(filebase,'_alpha.');
plotfileraketru= strcat(filebase,'_raketrue.');
plotfilehpd    = strcat(filebase,'_hpd.');
plotfilehpdtot= strcat(filebase,'_hpdtot.');
plotfiletru    = strcat(filebase,'_tru.');
plotfiletrutot = strcat(filebase,'_trutot.');
plotfileslice  = strcat(filebase,'_slice.');
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
NFPMX = NVMX*NPV;
load(filename);AS = A;
clear A;
ms = AS(:,5:NFPMX+4);
NSMP  = size(ms,1);

%%
%% Load true model
%%
if(isyn >= 1);
  sl1tru = load('true_sl1.txt');
  sl2tru = load('true_sl2.txt');
  Vrtru = load('true_slowr.txt');
  sltottru = sqrt(sl1tru.^2+sl2tru.^2);
end;
if(idatfit == 1);
  dat = h5read(datfile,'/Observed_data/displacements');
  dat = dat';
  NTSMP = h5read(datfile,'/Observed_data/Ntraces');
  NDAT = length(dat);
  if(I_WP == 1);rep=dlmread(repfile);end
  if(I_GPS == 1)repGPS=dlmread(repfileGPS);end;
  if(IAR == 1);repar=dlmread(reparfile);end;
  %% For rupture contours:
  ffsrmn=load(vrmxfile);
  ffdel=load(ffdelfile);
end;
if(ICOV >= 3);
  load(covfile);
  if(inonstat == 1);sd = dlmread(sdfile);end;
end;
if(ICOV == 2 | ICOV == 4);
  xi = AS(:,NVMX*NPV+1+4+NMISC+NSTN:NVMX*NPV+NSTN+NSTN+4+NMISC);
end;
%%
%% Compute residual errors:
%%
if(I_WP == 1);res = dat(1,:)-rep(3,:);end;
if(I_GPS == 1);resGPS = repGPS(3,:)-repGPS(2,:);end;

%%
%%  MAP model from file:
%%
map=dlmread(mapfile);
kmap = map(1,1);
mapvoro = map(2:end-3,:);
maphyp = map(end-2,1:2);

if(I_WP == 1);
for istn = 1:NSTN;
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  res(istart:iend) = res(istart:iend)-mean(res(istart:iend));
  resraw(istart:iend) = res(istart:iend)/std(res(istart:iend));
  if(ICOV >= 3);
    L3 = chol(F(istn).Cd);
    resstd(istart:iend) = inv(L3')*res(istart:iend)';
    if(ICOV == 4);
      resstd(istart:iend) = resstd(istart:iend)/sqrt(xi(istn));
    end;
  else
    if(inonstat == 0);
      resstd(istart:iend) = res(istart:iend)./sd(istn,1:NTSMP(i));
    else;
      resstd(istart:iend) = res(istart:iend)/std(res(istart:iend));
    end;
  end;
end; 
if(IAR == 1);
  resstd = dat(1,:)-rep(3,:)-repar(3,:);
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
end; 

x_ev = sfgrid(2,:);
z_ev = sfgrid(1,:);

minlimmisc = [ hyp_loc(1)-dhyp, hyp_loc(2)-dhyp, Vrmin];
maxlimmisc = [ hyp_loc(1)+dhyp, hyp_loc(2)+dhyp, Vrmax];

minlimar = -0.5;
maxlimar =  1.0;

thinstep = 1;
logLmin = min(AS(:,1))-(max(AS(:,1))-min(AS(:,1)))/10;
logLmax = max(AS(:,1))+(max(AS(:,1))-min(AS(:,1)))/10;

%% Discard burnin:
k    = AS(:,4);
logL = AS(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute DIC
%%
ii = 1;
for i=min(AS(:,4)):max(AS(:,4));
  idx=find(AS(:,4)==i);
  [logLmx(ii),ilogLmx(ii)]=max(AS(idx,1));
  numk(ii) = length(idx);
  logLmxidx(ii) = idx(ilogLmx(ii));
  %disp([i,logLmx(ii),ilogLmx(ii),numk(ii)]);
  ii = ii + 1;
  clear idx;
end;
D = -2.*logL;
[tmp,idxkhat] = max(numk);
Dmhat = D(logLmxidx(idxkhat));
Dbar = mean(D);
Pd = Dbar - Dmhat;
DIC = Dmhat + 2 * Pd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute EBIC
%%
EBIC = Dbar + (Pd*log(NDAT));
ElogL = mean(logL);
EBIC2 = -2*ElogL + (Pd*log(NDAT));
IC = -2*ElogL + 2*Pd;

%disp([DIC,EBIC,EBIC2,IC]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figk=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
idx = find(AS(:,end)==1);
stairs([1:length(idx)],AS(idx,4),'k')
ylabel('No. nodes in partition');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(idx)])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(AS(:,4),[1:1:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. nodes');
ylabel('Probability');
set(gca,'XLim',[0.5 NVMX+0.5]);
set(gca,'YLim',[0.   1.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%
figlogL=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)

%for i=1:max(AS(:,end));
i=1;
  idx = find(AS(:,end)==i);
  plot(AS(idx(1:thinstep:end),1),'k');
  clear idx;
%   plot([1:thinstep:length(AS(:,1))],AS(1:thinstep:end,1),'k');
%end;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(AS(:,1))])
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(AS(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
set(gca,'YLim',[logLmin logLmax])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Acceptance rate plot
%%
nx = 3;
ny = 2;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
 
figaccept=figure;hold on;box on;
title('Acceptance Rate');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
idx = find(AS(:,end)==1);
j = 1;
for i=1:NPV;
  subplot('Position',[loc(1,j) loc(2,j) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(AS(idx,4+NVMX*NPV+NMISC+NSTN+NSTN+NSDGPS+i),'k');
  plot([1:length(idx)],0.2*ones(size(idx)),'--k');
  plot([1:length(idx)],0.3*ones(size(idx)),'--k');

  if(j==1 | j==4);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  if(j>3);xlabel('rjMCMC step');else;set(gca,'XTickLabel',[]);end;
  set(gca,'XLim',[0 length(idx)],'YLim',[0 0.5])
  j=j+1;if(j==3);j=j+1;end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Hypocentre plot
%%
NPV*NVMX+5
hyp = AS(:,NVMX*NPV+1+4:NVMX*NPV+4+NMISC);
save tmp hyp;

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
[n,c]=hist3(hyp(:,1:2),{lim1,lim2});
n = n';
imagesc(c{1},c{2},n);
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[minlimmisc(1) maxlimmisc(1)],'YLim',[minlimmisc(2) maxlimmisc(2)]);
set(gca,'TickDir','out');
clear n lim1 lim2;

if(NMISC == 3);
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

  alpha = AS(:,NVMX*NPV+1+4+NMISC:NVMX*NPV+NSTN+4+NMISC);
  for i=1:NSTN;
    subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    [n,lim]=hist(alpha(:,i),30);n = [0, n, 0];lim = [lim(1) lim lim(end)];
    n = n/sum(n);
    lim = lim-(lim(3)-lim(2))/2.;
    [xx,yy]=stairs(lim,n,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n,'k');
    clear n lim;

    set(gca,'YLim',[0 .3],'XLim',[minlimar maxlimar]);
    xlabel('AR(1)');
    set(gca,'XTick',[-1:.2:1]);
    if(i == 1);
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

  for i=1:NSTN;
    subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    idxar = ones(size(alpha(:,i)));
    idx = find(alpha(:,i) < -0.5);
    idxar(idx) = 0;
    lim = [-1:1:2];
    [n]=hist(idxar,lim);%n = [0, n, 0];lim = [lim(1) lim lim(end)];
    n = n/sum(n);
    lim = lim - (lim(2)-lim(1))/2;
    [xx,yy]=stairs(lim,n,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n,'k');
    clear n lim;
    if(i==1);ylabel('Prob. AR(1) off/on');end;
    xlabel('');
    if(i>1);set(gca,'YTickLabel',[]);end;
    xticklabel = ({'off','on'});
    set(gca,'XTick',[0,1],'XTickLabel',xticklabel);
    set(gca,'XLim',[-1 2])
    set(gca,'YLim',[0 1.1])
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Plot and write MAP model to file:
%%
x  = mapvoro(:,1);
z  = mapvoro(:,2);
xn = mapvoro(:,1);
zn = mapvoro(:,2);
sl1  = mapvoro(:,3);
sl2  = mapvoro(:,4);
Vr  = 1./mapvoro(:,5);
%% Plot MAP
NSF
size(z_ev)
for iz = 1:NSF;
  dz = z_ev(iz)-zn;
  dx = x_ev(iz)-xn;
  d  = sqrt(dx.*dx + dz.*dz);
  [dmin,iv] = min(d);
  sl1_evmap(iz) = mapvoro(iv,3);
  sl2_evmap(iz) = mapvoro(iv,4);
  Vr_evmap(iz) = 1./mapvoro(iv,5);
end;
%%
%% Load previously computed and converted ensembles and histograms
%%
sl1_ens = zeros(NDEP,NRAN);
sl2_ens = zeros(NDEP,NRAN);
Vr_ens = zeros(NDEP,NRAN);
if(ialf == 1);alf_ens = zeros(NDEP,NRAN);end;

sl1_hst2 = zeros(NDEP,NRAN,NHIST);
sl2_hst2 = zeros(NDEP,NRAN,NHIST);
Vr_hst2 = zeros(NDEP,NRAN,NHIST);
if(ialf == 1);alf_hst2 = zeros(NDEP,NRAN,NHIST);end;

ffi_slab_convert_histfiles(filebase);
load hists.mat
%%
%% Compute 95% HPD map.
%%
sl1_hpd = zeros(size(sl1_ens));
sl2_hpd = zeros(size(sl2_ens));
sltot_hpd = zeros(size(sltot_ens));
alf_hpd = zeros(size(alf_ens));
Vr_hpd = zeros(size(Vr_ens));
nfsl1 = zeros(NDEP,NRAN,2);
nfsl2 = zeros(NDEP,NRAN,2);
nfVr = zeros(NDEP,NRAN,2);
nfsltot = zeros(NDEP,NRAN,2);
nfalf = zeros(NDEP,NRAN,2);
disp('Starting HPDs');

%% Compute CIs:
for ir = 1:NRAN;
  if(rem(ir,10) == 0);disp(ir);end;
  for iz = 1:NDEP;
    %% 95% HPD
    [nfsl1(iz,ir,:)] = hpd2(sl1_hst2(iz,ir,:),bins(1,:),95);
    sl1_hpd(iz,ir) = abs(nfsl1(iz,ir,2)-nfsl1(iz,ir,1));
    [nfsl2(iz,ir,:)] = hpd2(sl2_hst2(iz,ir,:),bins(2,:),95);
    sl2_hpd(iz,ir) = abs(nfsl2(iz,ir,2)-nfsl2(iz,ir,1));
    [nfsltot(iz,ir,:)] = hpd2(sltot_hst2(iz,ir,:),bins(4,:),95);
    sltot_hpd(iz,ir) = abs(nfsltot(iz,ir,2)-nfsltot(iz,ir,1));
    if(NMISC == 2);
      [nfVr(iz,ir,:)] = hpd2(Vr_hst2(iz,ir,:),bins(3,:),95);
      Vr_hpd(iz,ir) = abs(nfVr(iz,ir,2)-nfVr(iz,ir,1));
    end;
    if(ialf == 1);
      [nfalf(iz,ir,:)] = hpd2(alf_hst2(iz,ir,:),bins(5,:),95);
      alf_hpd(iz,ir) = abs(nfalf(iz,ir,2)-nfalf(iz,ir,1));
    end;
  end;
end;
disp('Done HPDs');
minslip = min(min(sltot_ens));
maxslip = max(max(sltot_ens));
minrake = min(min(alf_ens));
maxrake = max(max(alf_ens));
minnf1 = min(min(min(nfsl1)));
minnf2 = min(min(min(nfsl2)));
minnf3 = min(min(min(nfVr)));
minnft = min(min(min(nfsltot)));
maxnf1 = max(max(max(nfsl1)));
maxnf2 = max(max(max(nfsl2)));
maxnf3 = max(max(max(nfVr)));
maxnft = max(max(max(nfsltot)));
xplt
maxnftsl = max(max(max(nfsltot(:,xplt,:))));

minnfalf = min(min(min(nfalf)));
maxnfalf = max(max(max(nfalf)));

nx = 1;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

%%
%%
%%  MAP PLOTS
%%
%%

vq = griddata(z_ev,x_ev,sl1_evmap,yq,xq);
figure
mesh(xq,yq,vq);

figmap=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(yyq,xxq,vq');shading flat;
plot(z_ev,x_ev,'.w','Markersize',7);
plot(z_ev,x_ev,'.k','Markersize',4);

for i=1:kmap;
  plot(z(i),x(i),'dw','Markersize',7);
  plot(z(i),x(i),'dk','Markersize',5);
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
%plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 maxnf1],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'Map slip 1 (m)','FontSize',14);
%xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');

save tmp;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
vq = griddata(x_ev,z_ev,sl2_evmap,xq,yq);
imagesc(xxq,yyq,vq);shading flat;
for i=1:kmap;
  plot(x(i),z(i),'.w','Markersize',7)
  plot(x(i),z(i),'.k','Markersize',5)
end;
%[vx,vz]=voronoi(xn,zn);plot(vx*Rmx,vz*Zmx,'w','LineWidth',2);
%for i=1:NDEP;plot([x_ev(1) x_ev(end)],[z_ev(i) z_ev(i)],'-w');end;
%for i=1:NRAN;plot([x_ev(i) x_ev(i)],[z_ev(1) z_ev(end)],'-w');end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
colorbar;set(gca,'CLim',[0 maxnf2],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'MAP slip 2 (m)','FontSize',14);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse','TickDir','out')
if(NPV == 5);
  set(gca,'XTickLabel',[])
else;
  xlabel('Along-strike distance (km)');
end;
ylabel('Along-dip distance (km)');

if(NPV == 5);
subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
vq = griddata(x_ev,z_ev,Vr_evmap,xq,yq);
imagesc(xxq,yyq,vq);shading flat;
for i=1:kmap;
  plot(x(i),z(i),'.w','Markersize',7)
  plot(x(i),z(i),'.k','Markersize',5)
end;
[vx,vz]=voronoi(xn,zn);plot(vx*Rmx,vz*Zmx,'w','LineWidth',2);
for i=1:NDEP;plot([x_ev(1) x_ev(end)],[z_ev(i) z_ev(i)],'-w');end;
for i=1:NRAN;plot([x_ev(i) x_ev(i)],[z_ev(1) z_ev(end)],'-w');end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
v = [floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
%imagesc(ffsrmn-ffdel);hold on;
contour(x_ev,z_ev,ffsrmn-ffdel,v,'-w');
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'MAP rup. vel. (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
end;

%%
%%
%%  MAP PLOT TOTAL SLIP
%%
%%
figmaptot=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
if(isyn == 2);
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
else;
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
end;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sqrt(sl1_evmap.^2+sl2_evmap.^2));shading flat;
for i=1:kmap;
  plot(x(i),z(i),'.w','Markersize',7)
  plot(x(i),z(i),'.k','Markersize',5)
end;
[vx,vz]=voronoi(xn,zn);plot(vx*Rmx,vz*Zmx,'w','LineWidth',2);
for i=1:NDEP;plot([x_ev(1) x_ev(end)],[z_ev(i) z_ev(i)],'-w');end;
for i=1:NRAN;plot([x_ev(i) x_ev(i)],[z_ev(1) z_ev(end)],'-w');end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'TickDir','out')
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
if(isyn == 2);
  set(gca,'XLim',[225 475],'YLim',[0 150],'YDir','reverse')
  colorbar;set(gca,'CLim',[0 maxnft],'FontSize',14)
else
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  colorbar;set(gca,'CLim',[0 maxnft],'FontSize',14)
end;
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'Map slip 1 (m)','FontSize',14);

%clear AS;
t1 = tic;
voro = zeros(NVMX,5,NSMP);
nodes = zeros(NVMX*NSMP,2);
for ivo = 1:NVMX;
  idx(ivo) = (ivo-1)*NPV+1;
end;
idxk = [1:NSMP];
%idxk = find(k == 5);
for ismp = 1:length(idxk);
  for ivo = 1:NVMX;
    voro(ivo,1:NPV,ismp) = ms(idxk(ismp),idx(ivo):idx(ivo)+NPV-1);
  end;
  nodes((ismp-1)*NVMX+1:ismp*NVMX,:) = voro(:,1:2,ismp);
end;
clear idxk;
nodes(find(nodes(:,1) == 0),:)=[];

%%
%% Nodal density
%%
fignodes=figure();hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
nx = 1;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
nodes(:,1) = Rmx - nodes(:,1);
[hdens]=cloudPlot(nodes(:,1),nodes(:,2),[0 600 0 300],true,[601 301]);
%colormap( 1-hot );
colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YLim',[0 600],'YLim',[0 300]);
set(gca,'CLim',[0 1.5],'FontSize',14)
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;

%%
%% Compute moment:
%%
disp('Moment:');
%M0=sqrt(sum(sum(sl1_ens.*Nmscale))^2+sum(sum(sl2_ens.*Nmscale))^2)
M0=sum(sum(sqrt((sl1_ens.*Nmscale).^2+(sl2_ens.*Nmscale).^2)))
Mw = 2/3*log10(M0)-6

if(isyn > 0);
  disp('true Moment:');
  M0_t=sum(sum(sqrt((sl1tru.*Nmscale).^2+(sl2tru.*Nmscale).^2)))
  Mw_t = 2/3*log10(M0_t)-6.
end;

M0=load(M0file);M0(end)=[];
Mw=2/3*log10(exp(M0))-6;
figMw = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4])
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n1,xout]=hist(Mw,40);n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
area = sum(n1) * (xout(3)-xout(2));
n1 = n1/area;
[xx,yy]=stairs(xout,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(xout,n1,'k');
set(gca,'TickDir','out');box on;
xlabel('Moment magnitude M_w');
if(isyn > 0);
  yLimits = get(gca,'YLim');
  plot([Mw_t Mw_t],[0 yLimits(2)],'-k');
end;
clear xout n1;

nx = 1;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

nx = 1;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc2,spw2,sph2] = get_loc(nx,ny,xim,yim,xymarg);

%%
%% Ensemble mean figure slip components
%%
figens=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl1_ens);shading flat;
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 maxnf1],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
ylabel('Along-dip distance (km)');
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl2_ens);shading flat;
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
colorbar;set(gca,'CLim',[0 maxnf2],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. mean slip 2 (m)','FontSize',14);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse','TickDir','out')
if(NPV == 5);
  set(gca,'XTickLabel',[])
else;
  xlabel('Along-strike distance (km)');
end;
ylabel('Along-dip distance (km)');
box on;

if(NPV == 5);
subplot('Position',[loc2(1,3) loc2(2,3) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,Vr_ens);shading flat;
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
v = [floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
%imagesc(ffsrmn-ffdel);hold on;
contour(x_ev,z_ev,ffsrmn-ffdel,v,'-k','LineWidth',2);
contour(x_ev,z_ev,ffsrmn-ffdel,v,'--w','LineWidth',2);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. rup. vel. (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
%%
%% Ensemble mean figure total slip
%%
figenstot=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
if(isyn == 2);
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
else;
  subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
end;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sltot_ens);
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
if(isyn < 2);
vslip = [0:dvslip:maxslip];
contour(x_ev,z_ev,sltot_ens,vslip,'-w');
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
if(isyn == 2);
  set(gca,'XLim',[225 475],'YLim',[0 150],'YDir','reverse')
  colorbar;set(gca,'CLim',[0 maxnft],'FontSize',14)
else;
  %colorbar;set(gca,'CLim',[minlim(3) sqrt(2*maxlim(3)^2)],'FontSize',14)
  colorbar;set(gca,'CLim',[0 maxslip],'FontSize',14)
end;
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Slip magnitude (m)','FontSize',14);
box on;
if(NPV == 5);
subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,Vr_ens);shading flat;
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
v = [floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
%imagesc(ffsrmn-ffdel);hold on;
contour(x_ev,z_ev,ffsrmn-ffdel,v,'-k','LineWidth',2);
contour(x_ev,z_ev,ffsrmn-ffdel,v,'--w','LineWidth',2);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. rup. vel. (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
%%
%% Ensemble mean rake angle
%%
if(ialf == 1);
figalf=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
if(isyn == 2);
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
else;
  subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
end;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,alf_ens+alfmin);
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
colorbar;set(gca,'CLim',[minrake maxrake],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Rake angle (m)','FontSize',14);
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,alf_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnfalf maxnfalf],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Slip std. dev. (m)','FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
%%
%% Ensemble mean quiver plot
%%
%% Slip 1 is rake+45 degree; rotate to right by amount -rake(2)
%% Slip 2 is rake-45 degree; rotate to right by amount -rake(2)
[xsl,ysl] = rot2d(sl2_ens,sl1_ens,-rake(2));

figrake=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
ysl = -ysl;
%quiverc(x_ev,z_ev,xsl,ysl);
scale = 1;
[minlim(3) sqrt(2*maxlim(3)^2)]
hh=quiverwcolorbar(x_ev,z_ev,xsl,ysl,scale,'bounds',[minlim(3) sqrt(2*maxlim(3)^2)]);
set(hh,'LineWidth',1)
set(gca,'XLim',[0-deltr Rmx+deltr],'YLim',[0-2*deltz Zmx+deltz],'YDir','reverse')
set(gca,'TickDir','out')
%colorbar;set(gca,'CLim',[minslip maxslip],'FontSize',14)
colorbar;set(gca,'CLim',[0 maxslip],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)



figtmp=figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 14])
hq1 = quiver(x_ev,z_ev,xsl,ysl,1.5);
%quiver plots
figrake2=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 14])
%left version (regular)
subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)

%get the line position (first handle)
hkid = get(hq1,'children');
X = get(hkid(1),'XData');
Y = get(hkid(1),'YData');
close(figtmp);

jj = 1;
for ii = 1:3:length(X)-1
  len(jj) = sqrt((Y(ii+1)-Y(ii))^2+(X(ii+1)-X(ii))^2); %get the angle
  jj = jj + 1;
end;
len2 = round(len/max(len)*127)+1;
cmap = jet(128); %colormap, 116 because angles goes up to 115 degrees

jj = 1;
for ii = 1:3:length(X)-1
  %headWidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2);
  headwidth = 60*len2(jj);
  headlength = 60*len2(jj);
  ah = annotation('arrow',...
      'Color', cmap(len2(jj),:),'LineWidth',1.5,...
      'HeadStyle','plain','HeadLength',headlength,'HeadWidth',headwidth);
  set(ah,'parent',gca);
  set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
  jj = jj + 1;
end
set(hh,'LineWidth',1)
set(gca,'XLim',[0-deltr Rmx+deltr],'YLim',[0-2*deltz Zmx+deltz],'YDir','reverse')
set(gca,'TickDir','out')
%colorbar;set(gca,'CLim',[minslip maxslip],'FontSize',14)
colorbar;set(gca,'CLim',[0 maxslip],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)



%%
%% Ensemble mean quiver plot true parameters
%%
%% Slip 1 is rake+45 degree; rotate to right by amount -rake(2)
%% Slip 2 is rake-45 degree; rotate to right by amount -rake(2)
if(isyn >= 1);
  clear xsl ysl;
  [xsl,ysl] = rot2d(sl2tru,sl1tru,-rake(2));

  figraketru=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

  subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  ysl = -ysl;
  %quiverc(x_ev,z_ev,xsl,ysl);
  scale = 1;
  hh=quiverwcolorbar(x_ev,z_ev,xsl,ysl,scale,'bounds',[minlim(3) maxlim(3)]);
  set(hh,'LineWidth',1)
  set(gca,'XLim',[0-deltr Rmx+deltr],'YLim',[0-2*deltz Zmx+deltz],'YDir','reverse')
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minslip maxslip],'FontSize',14)
  colorbar;set(gca,'CLim',[0 maxslip],'FontSize',14)
  cb = colorbar('peer',gca,'FontSize',14);
  set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
  xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
  box on;
end;

%%
%% HPD figure
%%
fighpd=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl1_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 maxnf1-minnf1],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl2_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
colorbar;set(gca,'CLim',[0 maxnf2-minnf2],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (m)','FontSize',14);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse','TickDir','out')
if(NPV == 5);
  set(gca,'XTickLabel',[])
else;
  xlabel('Along-strike distance (km)');
end;
ylabel('Along-dip distance (km)');
box on;

if(NPV == 5);
subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,Vr_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[0 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
%%
%% Total slip 95% HPD
%%
fighpdtot=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sltot_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
%colorbar;set(gca,'CLim',[0 maxnft-minnft],'FontSize',14)
colorbar;set(gca,'CLim',[0 maxslip],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Slip std. dev. (m)','FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
if(NPV == 5);
subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,Vr_hpd);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[0 1],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% Credibility (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
%%
%% True figure
%%
if(isyn >= 1);
figtru=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl1tru);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 maxnf1],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'True slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sl2tru);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
colorbar;set(gca,'CLim',[0 maxnf2],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'True slip 2 (m)','FontSize',14);
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse','TickDir','out')
if(NPV == 5);
  set(gca,'XTickLabel',[])
else;
  xlabel('Along-strike distance (km)');
end;
ylabel('Along-dip distance (km)');
box on;

if(NPV == 5);
subplot('Position',[loc2(1,3) loc2(2,3) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,1./Vrtru);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'True rup. vel. (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
end;
%%
%% True figure total slip
%%
if(isyn > 0);
figtrutot=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
if(isyn == 2);
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
else;
  subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
end;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,sltottru);
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
if(isyn < 2);
vslip=[0:dvslip:maxslip];
contour(x_ev,z_ev,sltottru,vslip,'-w');
end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
set(get(cb,'ylabel'),'String', 'Ens. mean slip 1 (m)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
if(isyn == 2);
  set(gca,'XLim',[225 475],'YLim',[0 150],'YDir','reverse')
  colorbar;set(gca,'CLim',[0 maxnft],'FontSize',14)
else;
  %colorbar;set(gca,'CLim',[minlim(3) sqrt(2*maxlim(3)^2)],'FontSize',14)
  colorbar;set(gca,'CLim',[0 maxnft],'FontSize',14)
end;
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Slip magnitude (m)','FontSize',14);
box on;
if(NPV == 5);
subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,1./Vrtru);shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)
for i=1:length(xplt);
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
end;
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'True rup. vel. (km/s)','FontSize',14);
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;
end;
end;

%%
%% Slice figure slip
%%
figslice=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
if(isyn == 2);
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
  nx = length(xplt);
  ny = 1;
  xim = 0.01;
  yim = 0.15/ny;
  xymarg = [0.1 0.04 0.04 0.14];
else;
  nx = length(xplt);
  ny = 3;
  xim = 0.01;
  yim = 0.15/ny;
  xymarg = [0.07 0.04 0.04 0.14];
end;
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([sl1_ens(1,xplt(i));sl1_ens(:,xplt(i));sl1_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs([sl1tru(1,xplt(i));sl1tru(:,xplt(i));sl1tru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfsl1(1,xplt(i),1);nfsl1(:,xplt(i),1);nfsl1(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfsl1(1,xplt(i),2);nfsl1(:,xplt(i),2);nfsl1(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[0 maxnf1]);
  if(isyn == 2);
    set(gca,'YDir','reverse','YLim',[0 150],'XLim',[0 50]);
  end;
  xlabel('Slip 1 (m)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along-dip distance (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:dvslip2:150]);
  text(1700,110,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
if(isyn < 2);
for i=1:length(xplt);
  subplot('Position',[loc(1,i+length(xplt)) loc(2,i+length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([sl2_ens(1,xplt(i));sl2_ens(:,xplt(i));sl2_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs([sl2tru(1,xplt(i));sl2tru(:,xplt(i));sl2tru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfsl2(1,xplt(i),1);nfsl2(:,xplt(i),1);nfsl2(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfsl2(1,xplt(i),2);nfsl2(:,xplt(i),2);nfsl2(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[0 maxnf2]);
  xlabel('Slip 2 (m)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along-dip distance (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:dvslip2:150]);
box on;
end;
if(NPV == 5);
for i=1:length(xplt);
  subplot('Position',[loc(1,i+2*length(xplt)) loc(2,i+2*length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([Vr_ens(1,xplt(i));Vr_ens(:,xplt(i));Vr_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs(1./[Vrtru(1,xplt(i));Vrtru(:,xplt(i));Vrtru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfVr(1,xplt(i),1);nfVr(:,xplt(i),1);nfVr(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfVr(1,xplt(i),2);nfVr(:,xplt(i),2);nfVr(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[minnf3 maxnf3]);
  xlabel('Rup. vel. (km/s)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along-dip distance (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:.4:4]);
box on;
end;
end;
end;
%%
%% Slice figure total slip
%%
figslicetot=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8])
  nx = length(xplt);
  ny = 2;
  xim = 0.01;
  yim = 0.15/ny;
  xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([sltot_ens(1,xplt(i));sltot_ens(:,xplt(i));sltot_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs([sltottru(1,xplt(i));sltottru(:,xplt(i));sltottru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfsltot(1,xplt(i),1);nfsltot(:,xplt(i),1);nfsltot(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfsltot(1,xplt(i),2);nfsltot(:,xplt(i),2);nfsltot(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[0 maxnftsl]);
  %set(gca,'YDir','reverse','YLim',[0 150],'XLim',[0 maxnftsl]);
  if(isyn == 2);
    %set(gca,'YDir','reverse','YLim',[0 150],'XLim',[minnf1 maxnf1]);
    set(gca,'YDir','reverse','YLim',[0 150],'XLim',[0 50]);
  end;
  xlabel('Slip 1 (m)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along-dip distance (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:dvslip2:150]);
  text(1700,110,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
if(NPV == 5);
for i=1:length(xplt);
  subplot('Position',[loc(1,i+length(xplt)) loc(2,i+length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([Vr_ens(1,xplt(i));Vr_ens(:,xplt(i));Vr_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs(1./[Vrtru(1,xplt(i));Vrtru(:,xplt(i));Vrtru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfVr(1,xplt(i),1);nfVr(:,xplt(i),1);nfVr(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfVr(1,xplt(i),2);nfVr(:,xplt(i),2);nfVr(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[minnf3 maxnf3]);
  xlabel('Rup. vel. (km/s)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along-dip distance (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:.4:4]);
box on;
end;
end;

if(idatfit == 1);

if(I_WP == 1);
  %%
  %% Selected Data plots:
  %%
  figdatsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 2;
  ny = 2;
  xim = 0.01/nx;
  yim = 0.1/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
  idx2 = [];
  for j=1:length(idxdat);
    i = idxdat(j);
    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    idx2 = [idx2,[istart:iend]];
    sigma = std(dat(1,istart:iend)-rep(3,istart:iend));
    plot(deltt*[0:NTSMP(i)-1],dat(1,istart:iend),'b','LineWidth',1.5);
    plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend)-2*sigma,'k','LineWidth',1);
    plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend)+2*sigma,'k','LineWidth',1);
    plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend),'--r','LineWidth',1.5);
    plot([0 deltt*(NTSMP(i)-1)],[0 0],'--k','LineWidth',1);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(j == 1 | j == 3);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  for j=1:length(idxdat);
    i = idxdat(j);
    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
    set(gca,'YLim',[min(dat(1,idx2)) max(dat(1,idx2))],'XLim',[0 deltt*(NTSMP(i)-1)]);
  end;

  %%
  %% Full Data plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 5;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  jj = 1;
  for i = 1:NSTN;
    if(i == (nx*ny)+1);
      figdat2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    plot(deltt*[0:NTSMP(i)-1],dat(1,istart:iend),'b','LineWidth',1.5);
    plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend),'--r','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(dat) max(dat)],'XLim',[0 deltt*NTSMP(i)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:300:1600]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj + 1;
  end;
  %%
  %%  RESIDUALS
  %%
  figres = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figres2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    plot(deltt*[0:NTSMP(i)-1],resraw(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resraw) max(resraw)],'XLim',[0 deltt*NTSMP(i)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj + 1;
  end;
  if(IAR == 1);
  %%
  %%  TOTAL RESIDUALS
  %%
  figrestot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figres2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    plot(deltt*[0:NTSMP(i)-1],resstd(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resstd) max(resstd)],'XLim',[0 deltt*NTSMP(i)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj+1;
  end;
  %%
  %%  AR Predictions
  %%
  figarpred = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figarpred2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    plot(deltt*[0:NTSMP(i)-1],repar(3,istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(repar(3,:)) max(repar(3,:))],'XLim',[0 NTSMP(i)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj+1;
  end;
  end;
  if(ICOV >= 3);
  %%
  %%  TOTAL RESIDUALS
  %%
  figrestot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figrestot2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    L3 = chol(F(i).Cd); 

    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    plot(deltt*[0:NTSMP(i)-1],resstd(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resstd) max(resstd)],'XLim',[0 deltt*NTSMP(i)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj+1;
  end;
  end;
  %%
  %%  Autocovariance
  %%
  figacovr = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figacovr2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    axx = xcorr(resraw(istart:iend),'coeff');
    plot(deltt*[-NTSMP(i)+1:NTSMP(i)-1],axx,'k','LineWidth',1.5);
    plot(deltt*[-max(NTSMP(:)),max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axx;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    xlabel('lag (s)');
    set(gca,'XTick',[-1600:200:1600]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj+1;
  end;
  if(IAR == 1 | ICOV >= 3);
  %%
  %%  TOTAL Autocovariance
  %%
  figacovtot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figacoctot2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    axxtot = xcorr(resstd(istart:iend),'coeff');
    plot(deltt*[-NTSMP(i)+1:NTSMP(i)-1],axxtot,'k','LineWidth',1.5);
    plot([-deltt*max(NTSMP(:)),deltt*max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axxtot;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    xlabel('lag (s)');
    set(gca,'XTick',[-1600:200:1600]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jj = jj+1;
  end;
  end;

  figaxxsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
  idx2 = [];
  for j=1:length(idxdat);
    i = idxdat(j);
    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;

    axx = xcorr(resraw(istart:iend),'coeff');
    plot(deltt*[-NTSMP(i)+1:NTSMP(i)-1],axx,'k','LineWidth',1.5,'color', [0.6 0.6 0.6]);
    clear axx;

    axxtot = xcorr(resstd(istart:iend),'coeff');
    plot(deltt*[-NTSMP(i)+1:NTSMP(i)-1],axxtot,'k','LineWidth',1.5);
    plot(deltt*[-max(NTSMP(:)),max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axxtot;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    set(gca,'XTick',[-1600:400:1600]);
    if(j == 3 | j == 4);
      xlabel('lag (s)');
    else
      set(gca,'XTickLabel',[]);
    end;
    if(j == 1 | j == 3);
      ylabel('');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;

  %%
  %%  RESIDUAL HISTOGRAMS
  %%
  x = -6.2:.4:6.2;
  xx = -6.25:.01:6.25;
  figreshist = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figacoctot2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    [n1,xout] = hist(resraw(istart:iend),x);
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    stairs(xout,n1,'k');
    nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
    plot(xx,nd,'--k')
    set(gca,'XLim',[-6.5 6.5]);
    box on;
    if(jj >= (ny-1)*nx);
      xlabel('Res. (std. dev.)');
    end;
    jj = jj + 1;
  end;

  if(ICOV >= 3| IAR == 1);
  %%
  %%  TOTAL RESIDUAL HISTOGRAMS
  %%
  figrestothist = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jj = 1;
  for i=1:NSTN;
    if(i == (nx*ny)+1);
      figacoctot2 = figure();hold on;box on;
      jj = 1;
    end;
    subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    [n1,xout] = hist(resstd(istart:iend),x);
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    stairs(xout,n1,'k');
    nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
    plot(xx,nd,'--k')
    set(gca,'XLim',[-6.5 6.5]);
    box on;
    if(i >= (ny-1)*nx);
      xlabel('Res. (std. dev.)');
    end;
    jj = jj + 1;
  end;

  figressel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
  idx2 = [];
  for j=1:length(idxdat);
    i = idxdat(j);
    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;

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
    if(j == 3 | j == 4);
      xlabel('Residuals (std. dev.)');
    else
      set(gca,'XTickLabel',[]);
    end;
    if(j == 1 | j == 3);
      ylabel('');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  end;
end;
end;
if(isave == 1)
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext2),'-dpng');
  print(figmap,'-r250',strcat(plotfilemap,plotext2),'-dpng');
  print(fighyp,'-r250',strcat(plotfilehyp,plotext2),'-dpng');
  if(NMISC == 3);
    print(figrupvel,'-r250',strcat(plotfilerupvel,plotext2),'-dpng');
  end;
  print(figmaptot,'-r250',strcat(plotfilemaptot,plotext2),'-dpng');
  print(figens,'-painters','-r250',strcat(plotfileens,plotext2),'-dpng');
  print(figenstot,'-painters','-r250',strcat(plotfileenstot,plotext2),'-dpng');
  print(figrake2,'-painters','-r250',strcat(plotfilerake,plotext2),'-dpng');
  print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext2),'-dpng');
  print(figMw,'-painters','-r250',strcat(plotfileMw,plotext2),'-dpng');
  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext2),'-dpng');
  print(fighpdtot,'-painters','-r250',strcat(plotfilehpdtot,plotext2),'-dpng');
  if(ialf == 1);
    print(figalf,'-painters','-r250',strcat(plotfilealf,plotext2),'-dpng');
  end;
  if(isyn >= 1);
    print(figtru,'-painters','-r250',strcat(plotfiletru,plotext2),'-dpng');
    print(figtrutot,'-painters','-r250',strcat(plotfiletrutot,plotext2),'-dpng');
    print(figraketru,'-painters','-r250',strcat(plotfileraketru,plotext2),'-dpng');
  end;
  %print(figslice,'-painters','-r250',strcat(plotfileslice,plotext2),'-dpng');
  print(figslicetot,'-painters','-r250',strcat(plotfileslice,plotext2),'-dpng');
  print(figk,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
  if(I_WP == 1);
  print(figdat,'-painters','-r250',strcat(plotfiledata,plotext2),'-dpng');
  print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext2),'-dpng');
  print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext2),'-dpng');
  print(figres,'-painters','-r250',strcat(plotfileres,plotext2),'-dpng');
  print(figreshist,'-painters','-r250',strcat(plotfilereshist,plotext2),'-dpng');
  if(IAR == 1 | ICOV > 1);
    print(figrestot,'-painters','-r250',strcat(plotfilerestot,plotext2),'-dpng');
    print(figrestothist,'-painters','-r250',strcat(plotfilerestothist,plotext2),'-dpng');
    print(figacovr,'-painters','-r250',strcat(plotfileacovr,plotext2),'-dpng');
    print(figacovtot,'-painters','-r250',strcat(plotfileacovtot,plotext2),'-dpng');
    print(figressel,'-painters','-r250',strcat(plotfileressel,plotext2),'-dpng');
  end
  if(IAR == 1);
    print(figarpred,'-painters','-r250',strcat(plotfilearpred,plotext2),'-dpng');
  end
  end
  if(ieps == 1);
    print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-depsc');
    print(figmap,'-painters','-r250',strcat(plotfilemap,plotext3),'-depsc');
    print(fighyp,'-painters','-r250',strcat(plotfilehyp,plotext3),'-depsc');
    if(NMISC == 3);
      print(figrupvel,'-painters','-r250',strcat(plotfilerupvel,plotext3),'-deps');
    end;
    print(figmaptot,'-painters','-r250',strcat(plotfilemaptot,plotext3),'-depsc');
    print(figens,'-painters','-r250',strcat(plotfileens,plotext3),'-depsc');
    print(figenstot,'-painters','-r250',strcat(plotfileenstot,plotext3),'-depsc');
    print(figrake2,'-painters','-r250',strcat(plotfilerake,plotext3),'-depsc');
    print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext3),'-depsc');
    print(figMw,'-painters','-r250',strcat(plotfileMw,plotext3),'-depsc');
    print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext3),'-depsc');
    print(fighpdtot,'-painters','-r250',strcat(plotfilehpdtot,plotext3),'-depsc');
    if(isyn >= 1);
      print(figtru,'-painters','-r250',strcat(plotfiletru,plotext3),'-depsc');
      print(figtrutot,'-painters','-r250',strcat(plotfiletrutot,plotext3),'-depsc');
      print(figraketru,'-painters','-r250',strcat(plotfileraketru,plotext3),'-depsc');
    end;
    %print(figslice,'-painters','-r250',strcat(plotfileslice,plotext3),'-depsc');
    print(figslicetot,'-painters','-r250',strcat(plotfileslice,plotext3),'-depsc');
    print(figk,'-painters','-r250',strcat(plotfilek,plotext3),'-depsc');
    %print(figdat,'-painters','-r250',strcat(plotfiledata,plotext3),'-depsc');
    print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext3),'-depsc');
    print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext3),'-depsc');
    %print(figres,'-painters','-r250',strcat(plotfileres,plotext3),'-depsc');
    if(IAR == 1 | ICOV > 1);
      print(figressel,'-painters','-r250',strcat(plotfileressel,plotext3),'-depsc');
    end
    if(IAR == 1);
      %print(figrestot,'-painters','-r250',strcat(plotfilerestot,plotext3),'-depsc');
      %print(figarpred,'-painters','-r250',strcat(plotfilearpred,plotext3),'-depsc');
      %print(figacovr,'-painters','-r250',strcat(plotfileacovr,plotext3),'-depsc');
      %print(figacovtot,'-painters','-r250',strcat(plotfileacovtot,plotext3),'-depsc');
    end
  end;
end

return;
