function []=ffi_tsunami_plot_rjvoro(filename);

isave   = 1;
ifull   = 1;
ieps    = 0;
idatfit = 1;
isyn    = 0;
inonstat= 1;

gridfact = 10.;

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

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
NTW

NRAN2 = NRAN * gridfact;
NDEP2 = NDEP * gridfact;

%% Prior bounds:
minlim(1:2) = [ 0.,  0.];
maxlim(1:2) = [Rmx, Zmx];

%gaugelatlon(1,1:2) = [];

%hyp_loc = [120.,180.];
if(filebase(1:5) == 'maule')
  %% Maule
  idxdat = [11:6:NSTN];
  xplt = [8 15 25 35 40];
  dvslip = 5;
elseif(filebase(1:5) == 'tohok');
  %% Tohoku
  %idxdat = [floor(NSTN/5):floor(NSTN/5):NSTN]
  idxdat = [1:15];
  %idxdat = [   1, 4, 7,    12,14,  16];
  %ymin   = [-3.1,-6,-3.5,-1.5,-1,-0.4];
  %ymax   = [ 2.6, 6, 6.0, 4.2, 2, 1.0];
  if(NRAN <= 8)
    %hyp_loc = [ 55.,202.];
    xplt = [2 4 6];
  else
    %hyp_loc = [120.,180.];
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

fid = fopen('plotpar.txt');
%% Discard first 8 entries, then read NHIST
for i=1:12;
  s = fgets(fid);
  if(i == 9);NHIST = sscanf(s, '%d');end; % convert to number
  if(i ==11);minlim(3:4) = sscanf(s, '%f');end; % convert to number
  if(i ==12);maxlim(3:4) = sscanf(s, '%f');end; % convert to number
end;
minlim(NPV) = Vrmin;
maxlim(NPV) = Vrmax;
maxpert = maxlim-minlim;

%% Sub fault size in range and depth
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
%deltrG = double(Rmx)/double(NRAN2);
%deltzG = double(Zmx)/double(NDEP2);
%% Newton meter scaling for moment

deltrn = double(deltr)/double(Rmx);
deltzn = double(deltz)/double(Zmx);

% Output files
mapfile        = strcat(filebase,'_map.dat');
covfile        = strcat(filebase,'_covmat.mat');
sdfile         = strcat(filebase,'_nonstat_sd.dat');
reparfile      = strcat(filebase,'_replicaar.dat');
repfile        = strcat(filebase,'_replica.dat');
vrmxfile       = strcat(filebase,'_t_Vrmx.txt');
ffdelfile      = strcat(filebase,'_tdelay.txt');
plotfilek      = strcat(filebase,'_khist.');
plotfilelogL   = strcat(filebase,'_logL.');
plotfilemap1   = strcat(filebase,'_map1.');
plotfilemap2   = strcat(filebase,'_map2.');
plotfilehyp    = strcat(filebase,'_hyp.');
plotfilerupvel = strcat(filebase,'_rupvel.');
plotfilemaptot = strcat(filebase,'_maptot.');
plotfileens1   = strcat(filebase,'_ens1.');
plotfileens2   = strcat(filebase,'_ens2.');
%plotfileensG   = strcat(filebase,'_ensG.');
plotfileenstot = strcat(filebase,'_enstot.');
plotfilerake   = strcat(filebase,'_rake.');
plotfileraketru= strcat(filebase,'_raketrue.');
plotfilehpd    = strcat(filebase,'_hpd.');
plotfilehpd66  = strcat(filebase,'_hpd66.');
plotfilehpd66tot= strcat(filebase,'_hpd66tot.');
plotfiletru    = strcat(filebase,'_tru.');
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
  Vrtru = load('true_slowr.txt');
  sl1tru2 = zeros(NDEP,NRAN,NTW);
  for itw = 1:NTW;
    sl1tru2(:,:,itw) = sl1tru((itw-1)*NDEP+1:itw*NDEP,:);
  end;
  sl1tru_tot = sum(sl1tru2,3);
end;
if(idatfit == 1);
  dat = h5read(datfile,'/Observed_data/displacements');
  replin = h5read(datfile,'/Sensitivity_kernel/synthetic_displacements');
  size(replin)
  if(isyn == 1);
    dat = dat';
    replin = replin';
  end;
  NTSMP = cast(h5read(datfile,'/Observed_data/Ntraces'),'like',1);
  %NTSMP = [420,160,216,380,220,220,220,380,220,176,179,171,126,126,160,217,208];
  NDAT = length(dat);
  rep=dlmread(repfile);
  if(IAR == 1);repar=dlmread(reparfile);end;
  %% For rupture contours:
  ffsrmn=load(vrmxfile);
  ffdel=load(ffdelfile);
end;
if(ICOV >= 2);
  load(covfile);
  if(inonstat == 1);sd = dlmread(sdfile);end;
  if(ICOV == 3);
    xi = AS(:,NVMX*NPV+1+4+NMISC+NSTN:NVMX*NPV+NSTN+NSTN+4+NMISC);
  end;
end;
%%
%% Compute residual errors:
%%
if(idatfit == 1);
  res = dat(1,:)-rep(3,:);
end;
%%
%%  MAP model from file:
%%
map=dlmread(mapfile);
kmap = map(1,1);
mapvoro = map(2:end-3,:);
maphyp = map(end-2,1:2);

for istn = 1:NSTN;
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  res(istart:iend) = res(istart:iend)-mean(res(istart:iend));
  resraw(istart:iend) = res(istart:iend)/std(res(istart:iend));
  if(ICOV >= 2);
    L3 = chol(F(istn).Cd);
    resstd(istart:iend) = inv(L3')*res(istart:iend)';
    if(ICOV == 3);
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

x_evn = [deltrn/2.:deltrn:1.-deltrn/2.];
z_evn = [deltzn/2.:deltzn:1.-deltzn/2.];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';

%xG_ev = [deltrG/2.:deltrG:Rmx-deltrG/2.];
%zG_ev = [deltzG/2.:deltzG:Zmx-deltzG/2.];
%zG_ev = zG_ev';

minlimmisc = [ hyp_loc(1)-dhyp, hyp_loc(2)-dhyp, Vrmin];
maxlimmisc = [ hyp_loc(1)+dhyp, hyp_loc(2)+dhyp, Vrmax];

minlimar = -0.5;
maxlimar =  1.0;

dx_ev = x_ev(2)-x_ev(1);

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
  if(idx);
  [logLmx(ii),ilogLmx(ii)]=max(AS(idx,1));
  numk(ii) = length(idx);
  logLmxidx(ii) = idx(ilogLmx(ii));
  %disp([i,logLmx(ii),ilogLmx(ii),numk(ii)]);
  ii = ii + 1;
  end;
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
%stairs([1:length(idx)],AS(idx,4),'k')
stairs(AS(:,4),'k')
ylabel('No. nodes in partition');
xlabel('rjMCMC step');
set(gca,'XLim',[0 size(AS,1)])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(AS(:,4),[1:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. nodes');
ylabel('Probability');
set(gca,'XLim',[40.5 NVMX+0.5]);
set(gca,'YLim',[0.   1.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Count populated nodes:
%%
if(ifull == 1);
  t1_a = tic;
  kv = AS(:,4);
  gv = zeros(sum(kv),1);
  kvpar = zeros(size(kv,1),NPV);
  mv = AS(:,5:end-43);
  igv = 1;
  for i=1:size(AS,1);
    for j=1:kv(i);
      for itw=1:NTW;
        if(mv(i,(j-1)*NPV+2+itw) > -99.);
          kvpar(i,itw+2) = kvpar(i,itw+2) + 1;
          if(j>2);gv(igv) = gv(igv) + 1;end;
        end;
      end;
      if(NPV > 2+NTW);
      if(mv(i,(j-1)*NPV+2+NTW+1) > -99.);
        kvpar(i,2+NTW+1) = kvpar(i,2+NTW+1) + 1;
        if(j>2);gv(igv) = gv(igv) + 1;end;
      end;end;
      %% Position parameter always present:
      kvpar(i,1:2) = kvpar(i,1:2) + 1;
      igv = igv + 1;
    end;
  end;
  gv(find(gv==0))=[];
  nparv = sum(kvpar,2);
  t2 = toc(t1_a);
  disp('Done counting nodes. Time:');disp(t2);
  %%
  %% Save node file
  %%
  save tmp_nodes.mat kvpar gv nparv;
else
  %%
  %% Load pre-computed node file
  %%
  kv = AS(:,4);
  load tmp_nodes.mat;
end; %% end ifull
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% No. parameters PLOT
%%
fig_No_par=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
plot([1:length(nparv)],nparv,'k')
ylabel('No. parameters');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(nparv)],'TickDir','out');
set(gca,'YLim',[min(nparv)-1 max(nparv+1)],'TickDir','out');

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(nparv,[0:max(nparv)]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim - (lim(3)-lim(2))/2;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. parameters');
ylabel('Probability');
set(gca,'XLim',[min(nparv)-1 max(nparv)+1],'TickDir','out');
set(gca,'YLim',[ 0. 0.5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT per dimension
%%
fig_kperd=figure;
figw = 12;
figh = 6;
set(fig_kperd,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = NPV;
ny = 2;
xim = 0.01;
yim = 0.08;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for ipar=1:NPV;
  h1 = subplot('Position',[loc(1,ipar) loc(2,ipar) spw sph]);
  hold on; box off;
  set(gca,'FontSize',12);
  h2 = subplot('Position',[loc(1,ipar+NPV) loc(2,ipar+NPV) spw sph]);
  hold on; box off;
  set(gca,'FontSize',12);
  xbound = round(mean(kvpar));

  subplot(h1);
  set(gca,'FontSize',14);
  plot([1:thinstep:length(kvpar(:,ipar))],kvpar(1:thinstep:end,ipar),'k')
  if(i==1);ylabel('No. nodes');
  else;set(gca,'YTickLabel',[]);end;
  xlabel('rjMCMC step');
  set(gca,'YLim',[xbound(ipar)-10 xbound(ipar)+10],'TickDir','out');
  set(gca,'XLim',[0 length(k)],'TickDir','out');
  box on;

  subplot(h2);
  set(gca,'FontSize',14);
  [n,lim]=hist(kvpar(:,ipar),[0:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;

  if(ipar == 1);xlabel('No. nodes with X');end;
  if(ipar == 2);xlabel('No. nodes with Y');end;
  if(ipar > 2 & ipar < 2+NTW+1);xlabel('No. nodes with TW');end;
  if(ipar == 2+NTW+1);xlabel('No. nodes with Vrup');end;
  if(ipar == 1);ylabel('Probability');
  else;set(gca,'YTickLabel',[]);end;
%  if(isyn == 1);
%    plot([sum(voroidxtru(:,ipar+1)) sum(voroidxtru(:,ipar+1))],[0 1],'-w');
%    plot([sum(voroidxtru(:,ipar+1)) sum(voroidxtru(:,ipar+1))],[0 1],'--k');
%  end;
  set(gca,'XLim',[xbound(ipar)-10 xbound(ipar)+10],'TickDir','out');
  set(gca,'YLim',[ 0. 1.0]);
  box on;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%
figlogL=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)

%for i=1:max(AS(:,end));
i=1;
  %idx = find(AS(:,end)==i);
  %plot(AS(idx(1:thinstep:end),1),'k');
  plot(AS(:,1),'k');
  clear idx;
%   plot([1:thinstep:length(AS(:,1))],AS(1:thinstep:end,1),'k');
%end;

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
%% Acceptance rate plot
%%
nx = 5;
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
  plot(AS(idx,4+NVMX*NPV+NMISC+NSTN+NSTN+i),'k');
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
NMISC
NPV
hyp = AS(:,NVMX*NPV+1+4:NVMX*NPV+4+NMISC);

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
%[n,c]=hist3(hyp(:,1:2),{lim1,lim2});
%n = n';
%imagesc(c{1},c{2},n);
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
%% Plot MAP
for iz = 1:NDEP;
  dz = z_ev(iz)-z;
  for ir = 1:NRAN;
    if(mat_excl(iz,ir) == 1);
      dx = x_ev(ir)-x;
      d  = sqrt(dx.*dx + dz.*dz);
      [dmin,iv] = min(d);
      for itw = 1:NTW;
        if(mapvoro(iv,2+itw)<-99.);
          sl1_evmap(iz,ir,itw) = 0.;
        else;
          sl1_evmap(iz,ir,itw) = mapvoro(iv,2+itw);
        end;
      end;
    end;
  end;
end;

if(NPV > 2+NTW);
  kvrup = length(find(mapvoro(:,NPV)>-99.));
  mapvoro_vr = mapvoro(:,[1:2,NPV]);
  mapvoro_vr(find(mapvoro_vr(:,3)<-99),:) = [];
  x_vr  = mapvoro_vr(:,1);
  z_vr  = mapvoro_vr(:,2);
  for iz = 1:NDEP;
    dz = z_ev(iz)-z_vr;
    for ir = 1:NRAN;
      if(mat_excl(iz,ir) == 1);
        dx = x_ev(ir)-x_vr;
        d  = sqrt(dx.*dx + dz.*dz);
        [dmin,iv] = min(d);
        Vr_evmap(iz,ir) = 1./mapvoro_vr(iv,3);
      end;
    end;
  end;
end;
%%
%% Load previously computed and converted ensembles and histograms
%%
sl1_ens = zeros(NDEP,NRAN,NTW);
%sl1G_ens = zeros(NDEP2,NRAN2);
Vr_ens = zeros(NDEP,NRAN);

sl1_hst2 = zeros(NDEP,NRAN,NHIST,NTW);
Vr_hst2 = zeros(NDEP,NRAN,NHIST);

ffi_tsunami_convert_histfiles(filebase);
load hists.mat
nx = 2;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

save tmp
%%
%%
%%  MAP PLOTS
%%
%%
figmap1=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
axis equal;
imagesc(x_ev,z_ev,sum(sl1_evmap,3));shading flat;
plot(x,z,'.w','Markersize',7);
plot(x,z,'.k','Markersize',5);
[vx,vz]=voronoi(x,z);plot(vx,vz,'w','LineWidth',2);
for i=1:NDEP;plot([x_ev(1) x_ev(end)],[z_ev(i) z_ev(i)],'-w');end;
for i=1:NRAN;plot([x_ev(i) x_ev(i)],[z_ev(1) z_ev(end)],'-w');end;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
%plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'TickDir','out')
colormap(darkb2r(-1,6.5));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
cbfreeze(cb);
%cblabel('Map amplitude (m)');
xlabel('Easting (km)');
ylabel('Northing (km)');
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
freezeColors;

if(NPV > 2+NTW);
  figmap2=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(x_ev,z_ev,Vr_evmap);shading flat;
  plot(x_vr,z_vr,'.w','Markersize',7)
  plot(x_vr,z_vr,'.k','Markersize',5)
  [vx,vz]=voronoi(x_vr,z_vr);plot(vx,vz,'w','LineWidth',2);
  for i=1:NDEP;plot([x_ev(1) x_ev(end)],[z_ev(i) z_ev(i)],'-w');end;
  for i=1:NRAN;plot([x_ev(i) x_ev(i)],[z_ev(1) z_ev(end)],'-w');end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  v=[floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
  %imagesc(ffsrmn-ffdel);hold on;
  contour(x_ev,z_ev,ffsrmn-ffdel,v,'-w');
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  set(gca,'TickDir','out')
  %colormap(jet);
  %colormap(darkb2r(0,3.5));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'MAP rup. vel. (km/s)','FontSize',14);
  %set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
  set(gca,'CLim',[0.5 3.5],'FontSize',14)
  xlabel('Easting (km)');ylabel('Northing (km)');
end;

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
save('nodes.mat','nodes')
disp('saved nodes')

%%
%% Nodal density
%%
fignodes=figure();hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
nx = 1;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
axis equal;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
nodes(:,1) = Rmx - nodes(:,1);
[hdens]=cloudPlot(Rmx-nodes(:,1),nodes(:,2),[0 Rmx 0 Zmx],true,[601 301]);
%colormap( 1-hot );
colormap( 1-gray(256) );
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
set(gca,'CLim',[0 1.5],'FontSize',14)
set(get(cb,'ylabel'),'String', 'Node density','FontSize',14);
xlabel('Easting (km)');ylabel('Northing (km)');
box on;
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx]);

nx = 4;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

%%
%% Compute 95% HPD map.
%%
sl1_hpd = zeros(size(sl1_ens));
sl1_tot_hpd = zeros(size(Vr_ens));
Vr_hpd = zeros(size(Vr_ens));
disp('Starting HPDs');

%% Compute CIs:
for itw = 1:NTW;
for ir = 1:NRAN;
  if(rem(ir,50) == 0);disp(ir);end;
  for iz = 1:NDEP;
    %% 95% HPD
    [nfsl1(iz,ir,itw,:)] = hpd2(sl1_hst2(iz,ir,:,itw),bins(1,:),95);
    sl1_hpd(iz,ir,itw) = abs(nfsl1(iz,ir,itw,2)-nfsl1(iz,ir,itw,1));
    if(NPV == 2+NTW+1);
      [nfVr(iz,ir,:)] = hpd2(Vr_hst2(iz,ir,:),bins(2,:),95);
      Vr_hpd(iz,ir) = abs(nfVr(iz,ir,2)-nfVr(iz,ir,1));
    end;
    %[nfsltot(iz,ir,:)] = hpd2(sltot_hst2(iz,ir,:),bins(4,:),95);
    %sltot_hpd(iz,ir) = abs(nfsltot(iz,ir,2)-nfsltot(iz,ir,1));
  end;
end;
end;
%% Compute CIs for total slip:
for ir = 1:NRAN;
  if(rem(ir,50) == 0);disp(ir);end;
  for iz = 1:NDEP;
    %% 95% HPD
    [nfsl1_tot(iz,ir,:)] = hpd2(sl1_tot_hst2(iz,ir,:),bins(1,:),95);
    sl1_tot_hpd(iz,ir) = abs(nfsl1_tot(iz,ir,2)-nfsl1_tot(iz,ir,1));
  end;
end;

disp('Done HPDs');
minslip = min(min(min(sl1_ens)))
maxslip = max(max(max(sl1_ens)))
%minslipG = min(min(sl1G_ens))
%maxslipG = max(max(sl1G_ens))
minnf1 = min(min(min(nfsl1)))
maxnf1 = max(max(max(nfsl1)))
if(NPV > NTW+2);
  minnf3 = min(min(min(nfVr)))
  maxnf3 = max(max(max(nfVr)))
end;


%%
%% Ensemble mean figure slip components
%%
figens1=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

for itw=1:NTW;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(x_ev,z_ev,sl1_ens(:,:,itw));shading flat;
  for i=1:length(xplt);
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
  end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
  set(gca,'CLim',[-1 6.5],'FontSize',14)
  colormap(darkb2r(-1,6.5));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  cbfreeze(cb);
  %cblabel('Posterior mean (m)','FontSize',14);
  set(get(cb,'ylabel'),'String', 'Posterior mean (m)','FontSize',14);
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  freezeColors;
end;

if(NPV > 2+NTW);
  figens2=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  axis equal;
  imagesc(x_ev,z_ev,Vr_ens);shading flat;
  for i=1:length(xplt);
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
  end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  v=[floor(min(min(ffsrmn-ffdel))):10:ceil(max(max(ffsrmn-ffdel)))];
  %imagesc(ffsrmn-ffdel);hold on;
  contour(x_ev,z_ev,ffsrmn-ffdel,v,'-k','LineWidth',2);
  contour(x_ev,z_ev,ffsrmn-ffdel,v,'--w','LineWidth',2);
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  set(gca,'TickDir','out')
  %colormap(jet);
  %set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
  set(gca,'CLim',[0.5 3.5],'FontSize',14)
  %colormap(darkb2r(0,3.5));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'Ens. rup. vel. (km/s)','FontSize',14);
  xlabel('Easting (km)');ylabel('Northing (km)');
  box on;
end;
%%
%% Ensemble mean figure including Gaussians
%%
%figensG=figure();hold on;box on;
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
%
%subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
%set(gca,'FontSize',14,'layer','top','LineWidth',1)
%axis equal;
%imagesc(xG_ev,zG_ev,sl1G_ens);shading flat;
%for i=1:length(xplt);
%  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
%  plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
%end;
%plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
%%plot(maphyp(1),maphyp(2),'*w','Markersize',16)
%set(gca,'TickDir','out')
%%colorbar;set(gca,'CLim',[minslipG maxslipG],'FontSize',14)
%set(gca,'CLim',[-1.5 4],'FontSize',14)
%colormap(darkb2r(-1.5,4));
%cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
%set(get(cb,'ylabel'),'String', 'Posterior mean amplitude (m)','FontSize',14);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%box on;
%set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')


%%
%% HPD figure
%%
fighpd=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
axis equal;
imagesc(x_ev,z_ev,sl1_tot_hpd(:,:));shading flat;
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
%plot(maphyp(1),maphyp(2),'*w','Markersize',16)
set(gca,'TickDir','out')
colormap(flipud(gray));
set(gca,'CLim',[0 0.5],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
cbfreeze(cb);
%cblabel('95% Credibility (m)','FontSize',14);
xlabel('Easting (km)');ylabel('Northing (km)');
box on;
set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
freezeColors;

if(NPV == NTW+2+1);
  subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  axis equal;
  imagesc(x_ev,z_ev,Vr_hpd);shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[0 maxnf3-minnf3],'FontSize',14)
  colormap(flipud(gray));
  set(gca,'CLim',[0 1],'FontSize',14)
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', '95% Credibility (km/s)','FontSize',14);
  xlabel('Easting (km)');ylabel('Northing (km)');
  box on;
end;

%%
%% True figure
%%
if(isyn >= 1);
figtru=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])

for itw=1:NTW;
  subplot('Position',[loc(1,itw) loc(2,itw) spw sph]);hold on;box on;
  axis equal;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  imagesc(x_ev,z_ev,sl1tru2(:,:,itw));shading flat;
  for i=1:length(xplt);
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
  end;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  %plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  set(gca,'TickDir','out')
  %colorbar;set(gca,'CLim',[minnf1 maxnf1],'FontSize',14)
  set(gca,'CLim',[-1 6.5],'FontSize',14)
  colormap(darkb2r(-1,6.5));
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  cbfreeze(cb);
  %cblabel('True (m)','FontSize',14);
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  box on;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  freezeColors;
end;

if(NPV > 2+NTW);
  figtruslow=figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 11])
  subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  axis equal;
  imagesc(x_ev,z_ev,1./Vrtru);shading flat;
  plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
  plot(maphyp(1),maphyp(2),'*w','Markersize',16)
  for i=1:length(xplt);
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'-k')
    plot([x_ev(xplt(i)) x_ev(xplt(i))],[z_ev(1) z_ev(end)],'--w')
  end;
  set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx],'YDir','reverse')
  set(gca,'TickDir','out')
  set(gca,'CLim',[minnf3 maxnf3],'FontSize',14)
  cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
  set(get(cb,'ylabel'),'String', 'True rup. vel. (km/s)','FontSize',14);
  xlabel('Easting (km)');ylabel('Northing (km)');
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
  ny = 2;
  xim = 0.01;
  yim = 0.15/ny;
  xymarg = [0.1 0.04 0.04 0.14];
else;
  nx = length(xplt);
  ny = 2;
  xim = 0.01;
  yim = 0.15/ny;
  xymarg = [0.07 0.04 0.04 0.14];
end;
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([sl1_tot_ens(1,xplt(i));sl1_tot_ens(:,xplt(i));sl1_tot_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs([sl1tru_tot(1,xplt(i));sl1tru_tot(:,xplt(i));sl1tru_tot(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfsl1_tot(1,xplt(i),1);nfsl1_tot(:,xplt(i),1);nfsl1_tot(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfsl1_tot(1,xplt(i),2);nfsl1_tot(:,xplt(i),2);nfsl1_tot(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[-3 10]);
  if(isyn == 2);
    set(gca,'YDir','reverse','YLim',[0 150],'XLim',[0 50]);
  end;
  xlabel('Amplitude (m)');
  set(gca,'YTick',[0:50:Zmx]);
  if(i == 1);
    ylabel('Northing (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[-15:1:15],'XLim',[-3,10]);
  text(1700,110,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
if(NPV > 2+NTW);
for i=1:length(xplt);
  subplot('Position',[loc(1,i+length(xplt)) loc(2,i+length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs([Vr_ens(1,xplt(i));Vr_ens(:,xplt(i));Vr_ens(end,xplt(i))],[0;z_ev;Zmx],'-b','LineWidth',2);
  if(isyn>0);stairs(1./[Vrtru(1,xplt(i));Vrtru(:,xplt(i));Vrtru(end,xplt(i))],[0;z_ev;Zmx],'--r','LineWidth',2);end;
  stairs([nfVr(1,xplt(i),1);nfVr(:,xplt(i),1);nfVr(end,xplt(i),1)],[0;z_ev;Zmx],'--k','LineWidth',2);
  stairs([nfVr(1,xplt(i),2);nfVr(:,xplt(i),2);nfVr(end,xplt(i),2)],[0;z_ev;Zmx],'--k','LineWidth',2);
  %set(gca,'YDir','reverse','YLim',[0 Zmx],'XLim',[minnf3 maxnf3]);
  set(gca,'YDir','reverse','YLim',[0 Zmx]);
  xlabel('Rup. vel. (km/s)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Along strike (km)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[0:.4:4],'XLim',[0.5,3.5]);
  box on;
end;
end;
if(idatfit == 1);
  %%
  %% Selected Data plots:
  %%
  figdatsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 7])
  nx = 4;
  ny = 4;
  xim = 0.1/nx;
  yim = 0.1/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
  idx2 = [];
  deltt = double(deltt);
  for j=1:length(idxdat);
    i = idxdat(j);
    disp([i,j]);
    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:i));
    istart = iend-NTSMP(i)+1;
    idx2 = [idx2,[istart:iend]];
    sigma = std(dat(1,istart:iend)-rep(3,istart:iend));
    disp([istart,iend,deltt]);
    p1=plot(deltt*[0:NTSMP(i)-1],dat(1,istart:iend),'b','LineWidth',1.5);
    p2=plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend)-2*sigma,'--k','LineWidth',1);
    plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend)+2*sigma,'--k','LineWidth',1);
    p3=plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend),'-r','LineWidth',1.5);
    %plot(deltt*[0:NTSMP(i)-1],replin(1,istart:iend),'--y','LineWidth',1.5);
    plot([0 deltt*(NTSMP(i)-1)],[0 0],'--k','LineWidth',1);
    xlabel('time (s)');
    set(gca,'XTick',[0:400:7200]);
    if(j == 1 | j == 3);
      ylabel('Amplitude');
    else;
      %set(gca,'YTickLabel',[]);
    end;
    %set(gca,'YLim',[ymin(j),ymax(j)],'XLim',[0 deltt*NTSMP(i)]);
    set(gca,'XLim',[0 deltt*NTSMP(i)]);
    %set(gca,'YLim',[min([min(dat(1,istart:iend)),min(rep(3,istart:iend)),...
    %    min(replin(1,istart:iend)),min(rep(3,istart:iend)-2*sigma)]) ...
    %    max([max(dat(1,istart:iend)),max(rep(3,istart:iend)),...
    %    max(1,replin(istart:iend)),max(rep(3,istart:iend)+2*sigma)])],'XLim',[0 deltt*NTSMP(i)]);
    box on;
  end;
%  for j=1:length(idxdat);
%    i = idxdat(j);
%    subplot('Position',[loc_s(1,j) loc_s(2,j) spw_s sph_s]);hold on;box on;
%    set(gca,'YLim',[min(dat(1,idx2)) max(dat(1,idx2))],'XLim',[0 deltt*(NTSMP(i)-1)]);
%  end;
  legend([p1,p2,p3],'observed','95% CI','predicted');

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
    %plot(deltt*[0:NTSMP(i)-1],replin(1,istart:iend),'--y','LineWidth',1.5);
    plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min([min(dat(1,istart:iend)),min(rep(3,istart:iend)),...
        min(replin(1,istart:iend))]) max([max(dat(1,istart:iend)),...
        max(rep(3,istart:iend)),max(replin(1,istart:iend))])],'XLim',[0 deltt*NTSMP(i)]);
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
  legend('obs','pred','lin');
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
  if(ICOV >= 2);
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
  if(IAR == 1 | ICOV >= 2);
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

  if(ICOV >= 2| IAR == 1);
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
if(isave == 1)
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext2),'-dpng');
  print(figmap1,'-r250',strcat(plotfilemap1,plotext2),'-dpng');
  print(figmap2,'-r250',strcat(plotfilemap2,plotext2),'-dpng');
%  if(NMISC == 3);
%    print(figrupvel,'-r250',strcat(plotfilerupvel,plotext2),'-dpng');
%  end;
  print(figens1,'-painters','-r250',strcat(plotfileens1,plotext2),'-dpng');
  print(figens2,'-painters','-r250',strcat(plotfileens2,plotext2),'-dpng');
%  print(figensG,'-painters','-r250',strcat(plotfileensG,plotext2),'-dpng');
  print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext2),'-dpng');
  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext2),'-dpng');
  if(isyn >= 1);
    print(figtru,'-painters','-r250',strcat(plotfiletru,plotext2),'-dpng');
  end;
  print(figslice,'-painters','-r250',strcat(plotfileslice,plotext2),'-dpng');
  print(figk,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
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
  if(ieps == 1);
    print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-depsc');
    print(figmap,'-painters','-r250',strcat(plotfilemap,plotext3),'-depsc');
%    if(NMISC == 3);
%      print(figrupvel,'-painters','-r250',strcat(plotfilerupvel,plotext3),'-deps');
%    end;
    print(figens,'-painters','-r250',strcat(plotfileens,plotext3),'-depsc');
    print(figensG,'-painters','-r250',strcat(plotfileensG,plotext3),'-dpng');
    print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext3),'-depsc');
    print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext3),'-depsc');
    if(isyn >= 1);
      print(figtru,'-painters','-r250',strcat(plotfiletru,plotext3),'-depsc');
    end;
    print(figslice,'-painters','-r250',strcat(plotfileslice,plotext3),'-depsc');
    print(figk,'-painters','-r250',strcat(plotfilek,plotext3),'-depsc');
    print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext3),'-depsc');
    print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext3),'-depsc');
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
