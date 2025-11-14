%function []=cmt_plot_marginals(filename);

%filename = 'maule_data_sample.mat';
%filename = 'bo10_data_sample.mat';
%filename = 'ec10_data_sample.mat';
filename = 'ma10_data_sample.mat';
%filename = 'mi10_data_sample.mat';
%filename = 'am16_data_sample.mat';
%filename = 'cr17_data_sample.mat';
%filename = 'hi17_data_sample.mat';
%filename = 'kk17_data_sample.mat';
%filename = 'td17_data_sample.mat';
%filename = 'ak17_6.3_sample.mat';
%filename = 'ak17_5.7_sample.mat';
%filename = 'an17_data_sample.mat';

isave   = 1;
ieps    = 1;
idatfit = 0;
isyn    = 0;
inonstat= 1;
FNT = 21;

%% Was sampling done with double log prior on M0?
i_dlogprior = 1;

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,ILOC,IAR,IEXCHANGE,IDBLCPL,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay,minMw,minstr,mindip,minrak,...
 maxMw,maxstr,maxdip,maxrak]=cmt_read_parfile(parfile);
[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

IAR
%cmt_plot_cor(filename,NCMT,NLOC);
IDBLCPL
NCMT  = 5;
if(IDBLCPL == 1);
  NCMT2 = 4;
else
  NCMT2 = 5;
end;
NLOC = 4;
NPV = NCMT2+NLOC;
deltt = 4.;

size(minlat)
%% Prior bounds:
if(IDBLCPL == 1);
  minlim(1:NVMX,1:NCMT2) = [minMw,minstr,mindip,minrak];
  maxlim(1:NVMX,1:NCMT2) = [maxMw,maxstr,maxdip,maxrak];
else;
  minlim(1:NVMX,1:NCMT2) = -6;
  maxlim(1:NVMX,1:NCMT2) =  6;
end;

minlim(1:NVMX,NCMT2+1:NPV) = [minlat',minlon',mindepth',mindelay'];
maxlim(1:NVMX,NCMT2+1:NPV) = [maxlat',maxlon',maxdepth',maxdelay'];

%minlim
%maxlim
%NSTN

if(filebase(1:5) == 'maule')
  %% Maule
  idxdat = [15:16:NSTN];
elseif(filebase(1:5) == 'tohok');
  %% Tohoku
  idxdat = [10:round(NSTN/4):NSTN];
  %idxdat = [5:5:NSTN];
elseif(filebase(1:4) == 'hg12');
  %% Haida Gwaii 2012
  idxdat = [10:round(NSTN/4):NSTN];
  %idxdat = [5:5:NSTN];
elseif(filebase(1:4) == 'bo13');
  %% Bohol 2013
  idxdat = [4:round(NSTN/4):NSTN];
  %idxdat = [5:5:NSTN];
elseif(filebase(1:4) == 'am16');
  %% Amberley 2016
  idxdat = [10:round(NSTN/4):NSTN];
  %idxdat = [5:5:NSTN];
elseif(isyn >= 1);
  %% Sim:
  idxdat = [11:6:NSTN];
end;

% Output files
mapfile        = strcat(filebase,'_map.dat');
txtsamplefile  = strcat(filebase,'_sample.txt');
repfile        = strcat(filebase,'_replica.dat');
reparfile      = strcat(filebase,'_replicaar.dat');
covfile        = strcat(filebase,'_noise.mat');
sdfile         = strcat(filebase,'_nonstat_sd.dat');

plotfilelogL   = strcat(filebase,'_logL');
plotfilehyp    = strcat(filebase,'_hyp');
plotfilemarg   = strcat(filebase,'_1Dmarg');
plotfilemargnod= strcat(filebase,'_1Dmargnodal');
plotfileMw     = strcat(filebase,'_Mw');
plotfileMwtot  = strcat(filebase,'_Mwtot');
plotfiledata   = strcat(filebase,'_data');
plotfiledatasel= strcat(filebase,'_datasel');
plotfileaxxsel= strcat(filebase,'_autocovsel');
plotfileressel= strcat(filebase,'_ressel');
plotfileres    = strcat(filebase,'_residuals');
plotfilerestot = strcat(filebase,'_totalres');
plotfilereshist    = strcat(filebase,'_reshist');
plotfilerestothist = strcat(filebase,'_totalreshist');
plotfileacovr  = strcat(filebase,'_autocovraw');
plotfileacovtot= strcat(filebase,'_autocovtot');
plotfilek      = strcat(filebase,'_khist');
plotfilebb     = strcat(filebase,'_bb');
plotfilebb2    = strcat(filebase,'_bb2');
plotfilepctDC    = strcat(filebase,'_pctDC');

plotext1    = '.fig';
plotext2    = '.png';
plotext3    = '.eps';

%%
%% Load sample file
%%
NFPMX = NVMX*NPV;
load(filename);
AS = A;
clear A;
ms = AS(:,5:NFPMX+4);
NSMP  = size(ms,1);

%%
%% Load true model
%%
if(idatfit == 1);
  dat = dobs';
  rep=dlmread(repfile);
  if(IAR == 1);repar=dlmread(reparfile);end;
end;
if(ICOV == 2);
  load(covfile);
  if(inonstat == 1);sd = dlmread(sdfile);end;
  if(ICOV == 3);
    xi = AS(:,NVMX*NPV+1+4+NSTN:NVMX*NPV+NSTN+NSTN+4);
  end;
end;
if(ICOV == 4);
  xi = AS(:,NVMX*NPV+1+4+NSTN:NVMX*NPV+NSTN+NSTN+4);
end;
%%
%%  MAP model from file:
%%
map=dlmread(mapfile);
kmap = map(1,1);
mapvoro = map(2:end-3,:);
[map_logL,idxmap] = max(AS(:,1));

%%
%% Compute residual errors:
%%
if(idatfit == 1);
res = dat(1,:)-rep(3,:);
for istn = 1:NSTN;
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  res(istart:iend) = res(istart:iend)-mean(res(istart:iend));
  resraw(istart:iend) = res(istart:iend)/std(res(istart:iend));
  if(ICOV == 2);
    L3 = chol(n(istn).Cd);
    resstd(istart:iend) = inv(L3')*res(istart:iend)';
    if(ICOV == 3);
      resstd(istart:iend) = resstd(istart:iend)/sqrt(xi(istn));
    end;
  else
    if(inonstat == 0);
      resstd(istart:iend) = res(istart:iend)./sd(istn,1:NTSMP(istn));
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
minlimar = -0.5;
maxlimar =  1.0;

thinstep = 1;
logLmin = min(AS(:,1))-(max(AS(:,1))-min(AS(:,1)))/10;
logLmax = max(AS(:,1))+(max(AS(:,1))-min(AS(:,1)))/10;

k    = AS(:,4);
logL = AS(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figk=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
idx = find(AS(:,end)==1);
stairs([1:length(idx)],AS(idx,4),'k')
ylabel('No. nodes in partition');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(idx)])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim]=hist(AS(:,4),[1:1:NVMX]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
n1 = n1/sum(n1);
lim = lim-0.5;
[xx,yy]=stairs(lim,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n1,'k');
clear n1 lim;
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
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)

iproc=1;
idx = find(AS(:,end)==iproc);
plot(AS(idx(1:thinstep:end),1),'k');
clear idx;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(AS(:,1))])
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim]=hist(AS(:,1),100);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
n1 = n1/sum(n1);
[xx,yy]=stairs(n1,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n1,lim,'k');
clear n1 lim;
xlabel('Probability');
set(gca,'YLim',[logLmin logLmax])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Acceptance rate plot
%%
nx = 3;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
 
figaccept=figure;hold on;box on;
title('Acceptance Rate');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
idx = find(AS(:,end)==1);
jv = 1;
for iv=1:NPV;
  subplot('Position',[loc(1,jv) loc(2,jv) spw sph]);hold on;box on;
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  plot(AS(idx,4+NVMX*NPV+NSTN+NSTN+iv),'k');
  plot([1:length(idx)],0.2*ones(size(idx)),'--k');
  plot([1:length(idx)],0.3*ones(size(idx)),'--k');

  if(jv==1 | jv==4);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  if(jv>3);xlabel('rjMCMC step');else;set(gca,'XTickLabel',[]);end;
  set(gca,'XLim',[0 length(idx)],'YLim',[0 0.5])
  jv=jv+1;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Plot CMT position density
%%
voro = zeros(NSMP,NPV,NVMX);
nodes = zeros(3,NVMX*NSMP);
for ivo = 1:NVMX;
  idx(ivo) = (ivo-1)*NPV+1;
end;
idxk = [1:NSMP];
for ismp = 1:length(idxk);
  for ivo = 1:NVMX;
    voro(ismp,1:NPV,ivo) = squeeze(ms(idxk(ismp),idx(ivo):idx(ivo)+NPV-1));
  end;
  nodes(:,(ismp-1)*NVMX+1:ismp*NVMX) = squeeze(voro(ismp,6:8,:));
end;
save('cmt_pos.txt','nodes','-ascii');
clear idxk;

%%
%% Nodal density
%%
fignodes=figure();hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 6])
nx = 1;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
%[hdens]=cloudPlot(nodes(2,:),nodes(1,:),[min(minlon) max(maxlon) min(minlat) max(maxlat)],true,[601 601]);
[hdens]=cloudPlot(nodes(2,:),nodes(1,:),[min(nodes(2,:)) max(nodes(2,:)) min(nodes(1,:)) max(nodes(1,:))],true,[201 201]);
%colormap( 1-hot );
colormap( 1-gray(256) );
set(gca,'TickDir','out');
set(gca,'XLim',[min(minlon) max(maxlon)],'YLim',[min(minlat) max(maxlat)]);
set(gca,'CLim',[0 1.5],'FontSize',FNT)
xlabel('Along-strike distance (km)');ylabel('Along-dip distance (km)');
box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Centroid plot
%%
if(ILOC == 1);
hyp = AS(:,NCMT2+1+4:NCMT2+4+NLOC);

fighyp = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10])
aspect = 110.574/(111.320*cos(abs(minlim(NCMT2+1))*pi/180.));
maxpert = maxlim-minlim;
spw = 0.44;
sph = spw*(maxpert(NCMT2+1)/maxpert(NCMT2+2))*aspect;
spw2 = 0.15;
sph2 = 0.15;

loc = [0.1,  0.1+0.01+spw,  0.1, 0.1+0.01+spw;
       0.14+sph+.01, 0.14+sph+.01, 0.14, 0.14 ];

subplot('Position',[loc(1,1) loc(2,1) spw sph2]);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim1]=hist(hyp(:,2),[minlim(NCMT2+2):(maxlim(NCMT2+2)-minlim(NCMT2+2))/200:maxlim(NCMT2+2)]);n1 = [0, n1, 0];lim1 = [lim1(1) lim1 lim1(end)];
n1 = n1/sum(n1);
lim1 = lim1-(lim1(3)-lim1(2))/2.;
[xx,yy]=stairs(lim1,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim1,n1,'k');
clear n1;
set(gca,'XLim',[minlim(NCMT2+2) maxlim(NCMT2+2)]);
set(gca,'YLim',[0 0.08],'XTickLabel',[],'TickDir','out');

subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim2]=hist(hyp(:,1),[minlim(NCMT2+1):(maxlim(NCMT2+1)-minlim(NCMT2+1))/200:maxlim(NCMT2+1)]);n1 = [0, n1, 0];lim2 = [lim2(1) lim2 lim2(end)];
n1 = n1/sum(n1);
lim2 = lim2-(lim2(3)-lim2(2))/2.;
[xx,yy]=stairs(n1,lim2,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n1,lim2,'k');
clear n1;
set(gca,'TickDir','out');
set(gca,'XLim',[0 0.08],'YLim',[minlim(NCMT2+1) maxlim(NCMT2+1)]);
set(gca,'YTickLabel',[]);

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[edgesX2,edgesY2,N] = ndhist(hyp(:,2),hyp(:,1));
sanePColor(edgesX2,edgesY2,N);
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
set(gca,'TickDir','out');
set(gca,'XLim',[min(minlon) max(maxlon)],'YLim',[min(minlat) max(maxlat)]);
set(gca,'TickDir','out');
colormap(flipud(gray));
clear n1 lim1 lim2;
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

  alpha = AS(:,NVMX*NPV+1+4:NVMX*NPV+NSTN+4);
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figar1 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figar2 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    [n1,lim]=hist(alpha(:,istn),30);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
    n1 = n1/sum(n1);
    lim = lim-(lim(3)-lim(2))/2.;
    [xx,yy]=stairs(lim,n1,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n1,'k');
    clear n1 lim;

    set(gca,'YLim',[0 .3],'XLim',[minlimar maxlimar]);
    xlabel('AR(1)');
    set(gca,'XTick',[-1:.2:1]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;

  figarb = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figarb1 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figarb2 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)

    idxar = ones(size(alpha(:,istn)));
    idx = find(alpha(:,istn) < -0.5);
    idxar(idx) = 0;
    lim = [-1:1:2];
    [n1]=hist(idxar,lim);%n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
    n1 = n1/sum(n1);
    lim = lim - (lim(2)-lim(1))/2;
    [xx,yy]=stairs(lim,n1,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(lim,n1,'k');
    clear n1 lim;
    if(i==1);ylabel('Prob. AR(1) off/on');end;
    xlabel('');
    if(i>1);set(gca,'YTickLabel',[]);end;
    xticklabel = ({'off','on'});
    set(gca,'XTick',[0,1],'XTickLabel',xticklabel);
    set(gca,'XLim',[-1 2])
    set(gca,'YLim',[0 1.1])
    jjstn = jjstn+1;
  end;
end;

%%
%% Compute moment:
%%
disp('Moment & Percent DC:');
figMw = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4.5])
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
nx = NVMX;
ny = 1;
xim = 0.05/nx;
yim = 0.10/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
%M = zeros(NSMP,NCMT2,NVMX);
MM = zeros(NSMP,3,3,NVMX);
for ivo=1:NVMX;
  if(IDBLCPL == 0);
    if(i_dlogprior == 0)
      Mtmp1(:,:)=AS(:,(ivo-1)*NPV+5:ivo*NPV-NLOC+4)*10.^22;
    else
      % Sampling of the CMT parameters is in a "folded" log space
      Mtmp1(:,:)=sign(AS(:,(ivo-1)*NPV+5:ivo*NPV-NLOC+4)).*10.^( ...
                 abs(AS(:,(ivo-1)*NPV+5:ivo*NPV-NLOC+4)));
    end;
    Mtmp2 = -Mtmp1(:,1)-Mtmp1(:,2);
    M(:,:,ivo)=[Mtmp1(:,1:2),Mtmp2,Mtmp1(:,3:5)];
    MM(:,1,1,ivo) = M(:,1,ivo);
    MM(:,2,2,ivo) = M(:,2,ivo);
    MM(:,3,3,ivo) = M(:,3,ivo);
    
    MT(:,1,1,ivo) = M(:,1,ivo);
    MT(:,2,2,ivo) = M(:,2,ivo);
    MT(:,3,3,ivo) = M(:,3,ivo);
    MT(:,1,2,ivo) = M(:,4,ivo);
    MT(:,1,3,ivo) = M(:,5,ivo);
    MT(:,2,3,ivo) = M(:,6,ivo);
    MT(:,2,1,ivo) = M(:,4,ivo);
    MT(:,3,1,ivo) = M(:,5,ivo);
    MT(:,3,2,ivo) = M(:,6,ivo);
    
    M2=M(:,:,ivo).^2;
    M0(:,ivo) = sqrt(0.5*(M2(:,1)+M2(:,2)+M2(:,3))+M2(:,4)+M2(:,5)+M2(:,6));
    Mw(:,ivo) = 2./3.*(log10(M0(:,ivo))-9.10);  %% 

  else
    M(:,:,1) = AS(:,5:8);
    Mw = AS(:,(ivo-1)*NPV+5);
    M0 = 10.^(1.5*(Mw+6.07));
  end

  if(IDBLCPL == 0);
    for ismp=1:NSMP;
      detM(ismp,ivo) = det(squeeze(MM(ismp,:,:,ivo)));
      trM(ismp,ivo) = trace(squeeze(MM(ismp,:,:,ivo)));
    end;
    if(isyn > 0);
      disp('true Moment:');
      Mtmp = -Mtru(1)-Mtru(2);
      Mtru=[Mtru(1:2),Mtmp,Mtru(3:5)];
      Mtru=Mtru*10^22;
      M2=Mtru.^2;
      M0tru = sqrt(0.5*(M2(1)+M2(2)+M2(3))+M2(4)+M2(5)+M2(6));
      Mwtru = 2./3.*(log10(M0)-9.10);  %% In Nm, see Aki Richards 2002 and 
                                       %% footnote in Shearer 2009 pg 284
    end;
  end;
  subplot('Position',[loc_s(1,ivo) loc_s(2,ivo) spw_s sph_s]);hold on;box on;
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  [n1,xout]=hist(Mw(:,ivo),40);n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
  area = sum(n1) * (xout(3)-xout(2));
  n1 = n1/area;
  [xx,yy]=stairs(xout,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(xout,n1,'k');
  set(gca,'TickDir','out','YTickLabel',[]);box on;
  xlabel('Moment magnitude M_w');
  if(isyn > 0);
    yLimits = get(gca,'YLim');
    plot([Mw_t Mw_t],[0 yLimits(2)],'-k');
  end;
  clear xout n1 Mtmp1 Mtmp2;
end;
for ivo=1:NVMX;
  if(IDBLCPL == 0);
    %% Compute % DC
    for ismp=1:NSMP;
      e(:,ismp) = eig(squeeze(MT(ismp,:,:,ivo)));
      F = min(abs(e(:,ismp)))/max(abs(e(:,ismp)));
      pctDC(ismp) = 100*(1-2*F); 
    end;
    figpctDC = figure();hold on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    [n1,xout]=hist(pctDC,40);n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    [xx,yy]=stairs(xout,n1,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(xout,n1,'k');
    set(gca,'TickDir','out','YTickLabel',[]);box on;
    xlabel('Percent DC');

    figeig = figure();hold on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    for ie=1:3;
      subplot(1,3,ie);hold on;box on;
      [n1,xout]=hist(e(ie,:),40);n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
      area = sum(n1) * (xout(3)-xout(2));
      n1 = n1/area;
      [xx,yy]=stairs(xout,n1,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(xout,n1,'k');
      set(gca,'TickDir','out','YTickLabel',[]);box on;
      lab=strcat('e_',num2str(ie));
      xlabel(lab);
    end;

  end;
end;
%save tmp.mat detM trM MM;

figMwtot = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4])
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
M0tot = sum(M0,2);
Mwtot = 2./3.*(log10(M0tot)-9.10);  %% In Nm, see Aki Richards 2002 and 
[n1,xout]=hist(Mwtot,40);n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
area = sum(n1) * (xout(3)-xout(2));
n1 = n1/area;
[xx,yy]=stairs(xout,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(xout,n1,'k');
set(gca,'TickDir','out','YTickLabel',[]);box on;
xlabel('Total M_w (sum of point sources)');

%%
%% 1D Marginal plots:
%%
nx = 5;
ny = 2;
ML = .02; MR = .02; MB = .0; MT = .02;
SP = .02; PAD = 0; PB = 0.1; FNT = 14;
bins = [40,40,40,40,40,60,60,30,60];
for ivo=1:NVMX;
  figmarg(ivo) = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8])
  jloc = 1;
  if(IDBLCPL == 0);
    label = ({'M_{rr}','M_{tt}','M_{rt}','M_{rp}','M_{tp}','Latitude (^o)','Longitude (^o)','Depth (km)','Delay time (s)'});
  else;
    label = ({'M_w','Strike (^o)','Dip (^o)','Rake (^o)','Latitude (^o)','Longitude (^o)','Depth (km)','Delay time (s)'});
  end;
  for jpar=1:NPV;
    subaxis(ny,nx,jloc,'Spacing',SP,'Padding',PAD,'Paddingbottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
    [n1,xout] = hist(voro(:,jpar,ivo),bins(jpar));
    n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    [xx,yy]=stairs(xout,n1);
    patch(xx,yy,[0.7,0.7,0.7]);
    stairs(xout,n1,'k','LineWidth',1.5);
    set(gca,'YTickLabel',[]);
    xlabel(label(jpar));
    clear n1 xout;
    if(jpar == NCMT2);jloc = nx+1;end;
    jloc = jloc + 1;
  end;
end;

%%
%% Beach Ball plots:
%%
figbb = figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4])
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
nx = NVMX;
ny = 1;
xim = 0.05/nx;
yim = 0.10/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
for ivo=1:NVMX;
  subplot('Position',[loc(1,ivo) loc(2,ivo) spw sph]);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  if(IDBLCPL == 0)
    bb(M(idxmap,:,ivo), 0, 0, 1, 0, 'b',1,1);
  else;
    bb(M(idxmap,2:4,ivo), 0, 0, 1, 0, 'b',1,1);
  end;
  axis equal;
  set(gca,'XTick',[],'YTick',[]);
end;

figbb2=figure();
for ivo=1:NVMX;
  subplot('Position',[loc(1,ivo) loc(2,ivo) spw sph]);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  NMT = 2000;
  FBB= zeros(511,511);
  B1 = zeros(511,511);
  tic;
  for iM=10:10:NSMP;
%    iM
    if(rem(iM,10000) == 0);
        time = toc;
        disp('Focal Mechanism, time = ');
        disp([iM,time])
        tic;
    end
    if(IDBLCPL == 0)
      %% Fuzzy beachball:
      [xVec,yVec,B1]=focalmech(M(iM,:,ivo),0,0,1,0);
      %% Summing 
      FBB = FBB + B1;
    else
      [b(iM).xx,b(iM).yy,b(iM).X,b(iM).Y]=bb(M(iM,2:4,ivo), 0, 0, 1, 0, 'b',0,alph);
    end;
  end;
  %% Normalize by number of models (used stride of 10 on line 686)
  FBB = FBB/(NSMP/10);
  %% Need to mask outside of circle:
  msk = zeros(size(FBB));
  msk(7:505,7:505)=filldisc(249);
  
  figbb3=figure();hold on;
  pcolor(xVec,yVec,FBB.*msk);
  u = .5*cos(0:0.02:2*pi); w = .5*sin(0:0.02:2*pi);
  plot(u,w,'-k')
  colormap(flipud(gray));shading flat;
  axis equal;axis off;
  set(gca,'XTick',[],'YTick',[]);
end;
%%
%% Get nodal planes:
%%
for ivo=1:NVMX;
  if(IDBLCPL == 0);
    for jM=1:NSMP;
      %% ang -> s1, d1, r1, s2, d2, r2
      [ang(jM,1),ang(jM,2),ang(jM,3)] = ...
        mij2sdr(M(jM,1,ivo),M(jM,2,ivo),M(jM,3,ivo),M(jM,4,ivo),M(jM,5,ivo),M(jM,6,ivo));
      [ang(jM,4),ang(jM,5),ang(jM,6)] = AuxPlane(ang(jM,1),ang(jM,2),ang(jM,3));
    end;
  else
    ang = M(:,2:4,ivo);
  end;
%  idxstr = find(ang(:,1)>300.);
  %%
  %% 1D Marginal plots nodal planes:
  %%
  figmargnod(ivo) = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 3;
  ny = 2;
  xim = 0.06/nx;
  yim = 0.20/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc_s,spw_s,sph_s] = get_loc(nx,ny,xim,yim,xymarg);
  jloc = 1;
  label = ({'Strike (^o)','Dip (^o)','Rake (^o)','Aux Strike (^o)','Aux Dip (^o)','Aux Rake (^o)'});
  if(IDBLCPL == 1);
    NANG = 3;
  else;
    NANG = 6;
  end;
  for jpar=1:NANG;
    subplot('Position',[loc_s(1,jloc) loc_s(2,jloc) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)

%    [n1,xout] = hist(ang(idxstr,jpar),50);
    [n1,xout] = hist(ang(:,jpar),50);
    n1 = [0, n1, 0];xout = [xout(1) xout xout(end)];
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    [xx,yy]=stairs(xout,n1);
    patch(xx,yy,[0.7,0.7,0.7]);
    stairs(xout,n1,'k','LineWidth',1.5);
    set(gca,'YTickLabel',[]);
    xlabel(label(jpar));
    clear n1 xout;
    jloc = jloc + 1;
  end;
  %clear ang;
end;

if(idatfit == 1);
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
  for jstn=1:length(idxdat);
    istn = idxdat(jstn);
    subplot('Position',[loc_s(1,jstn) loc_s(2,jstn) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    idx2 = [idx2,[istart:iend]];
    sigma = std(dat(1,istart:iend)-rep(3,istart:iend));
    plot(deltt*[0:NTSMP(istn)-1],dat(1,istart:iend),'b','LineWidth',1.5);
    plot(deltt*[0:NTSMP(istn)-1],rep(3,istart:iend)-2*sigma,'k','LineWidth',1);
    plot(deltt*[0:NTSMP(istn)-1],rep(3,istart:iend)+2*sigma,'k','LineWidth',1);
    plot(deltt*[0:NTSMP(istn)-1],rep(3,istart:iend),'--r','LineWidth',1.5);
    plot([0 deltt*(NTSMP(istn)-1)],[0 0],'--k','LineWidth',1);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(jstn == 1 | jstn == 3);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  for jstn=1:length(idxdat);
    istn = idxdat(jstn);
    subplot('Position',[loc_s(1,jstn) loc_s(2,jstn) spw_s sph_s]);hold on;box on;
    set(gca,'YLim',[min(dat(1,idx2)) max(dat(1,idx2))],'XLim',[0 deltt*(NTSMP(istn)-1)]);
  end;

  %%
  %% Full Data plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 6;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  jjstn = 1;
  for istn = 1:NSTN;
    if(istn == (nx*ny)+1);
      figdat2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figdat3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],dat(1,istart:iend),'b','LineWidth',1.5);
    plot(deltt*[0:NTSMP(istn)-1],rep(3,istart:iend),'--r','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(dat) max(dat)],'XLim',[0 deltt*NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:300:1600]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn + 1;
  end;
  %%
  %%  RESIDUALS
  %%
  figres = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figres2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figres3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],resraw(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resraw) max(resraw)],'XLim',[0 deltt*NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn + 1;
  end;
  if(IAR == 1);
  %%
  %%  TOTAL RESIDUALS
  %%
  figrestot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figrestot2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figrestot3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],resstd(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resstd) max(resstd)],'XLim',[0 deltt*NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;
  %%
  %%  AR Predictions
  %%
  figarpred = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figarpred2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figarpred3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],repar(3,istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(repar(3,:)) max(repar(3,:))],'XLim',[0 NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;
  end;
  if(ICOV == 2);
  %%
  %%  TOTAL RESIDUALS
  %%
  figrestot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figrestot2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figrestot3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)

    L3 = chol(n(istn).Cd); 

    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],resstd(istart:iend),'k','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(resstd) max(resstd)],'XLim',[0 deltt*NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;
  end;
  %%
  %%  Autocovariance
  %%
  figacovr = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figacovr2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figacovr3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    axx = xcorr(resraw(istart:iend),'coeff');
    plot(deltt*[-NTSMP(istn)+1:NTSMP(istn)-1],axx,'k','LineWidth',1.5);
    plot(deltt*[-max(NTSMP(:)),max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axx;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    xlabel('lag (s)');
    set(gca,'XTick',[-1600:200:1600]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;
  if(IAR == 1 | ICOV == 2);
  %%
  %%  TOTAL Autocovariance
  %%
  figacovtot = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figacovtot2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figacovtot3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    axxtot = xcorr(resstd(istart:iend),'coeff');
    plot(deltt*[-NTSMP(istn)+1:NTSMP(istn)-1],axxtot,'k','LineWidth',1.5);
    plot([-deltt*max(NTSMP(:)),deltt*max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axxtot;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    xlabel('lag (s)');
    set(gca,'XTick',[-1600:200:1600]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn+1;
  end;
  end;

  figaxxsel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
  idx2 = [];
  for jstn=1:length(idxdat);
    istn = idxdat(jstn);
    subplot('Position',[loc_s(1,jstn) loc_s(2,jstn) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;

    axx = xcorr(resraw(istart:iend),'coeff');
    plot(deltt*[-NTSMP(istn)+1:NTSMP(istn)-1],axx,'k','LineWidth',1.5,'color', [0.6 0.6 0.6]);
    clear axx;

    axxtot = xcorr(resstd(istart:iend),'coeff');
    plot(deltt*[-NTSMP(istn)+1:NTSMP(istn)-1],axxtot,'k','LineWidth',1.5);
    plot(deltt*[-max(NTSMP(:)),max(NTSMP(:))],[0,0],'--k','LineWidth',1);
    clear axxtot;

    set(gca,'YLim',[-.4 1.1],'XLim',[-deltt*max(NTSMP(:)) deltt*max(NTSMP(:))]);
    set(gca,'XTick',[-1600:400:1600]);
    if(jstn == 3 | jstn == 4);
      xlabel('lag (s)');
    else
      set(gca,'XTickLabel',[]);
    end;
    if(jstn == 1 | jstn == 3);
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
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figreshist2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figreshist3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    [n1,xout] = hist(resraw(istart:iend),x);
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    stairs(xout,n1,'k');
    nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
    plot(xx,nd,'--k')
    set(gca,'XLim',[-6.5 6.5]);
    box on;
    if(jjstn >= (ny-1)*nx);
      xlabel('Res. (std. dev.)');
    end;
    jjstn = jjstn + 1;
  end;

  if(ICOV == 2| IAR == 1);
  %%
  %%  TOTAL RESIDUAL HISTOGRAMS
  %%
  figrestothist = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  jjstn = 1;
  for istn=1:NSTN;
    if(istn == (nx*ny)+1);
      figrestothist2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figrestothist3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    [n1,xout] = hist(resstd(istart:iend),x);
    area = sum(n1) * (xout(3)-xout(2));
    n1 = n1/area;
    stairs(xout,n1,'k');
    nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
    plot(xx,nd,'--k')
    set(gca,'XLim',[-6.5 6.5]);
    box on;
    if(istn >= (ny-1)*nx);
      xlabel('Res. (std. dev.)');
    end;
    jjstn = jjstn + 1;
  end;

  figressel = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
  idx2 = [];
  for jstn=1:length(idxdat);
    istn = idxdat(jstn);
    subplot('Position',[loc_s(1,jstn) loc_s(2,jstn) spw_s sph_s]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;

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
    if(jstn == 3 | jstn == 4);
      xlabel('Residuals (std. dev.)');
    else
      set(gca,'XTickLabel',[]);
    end;
    if(jstn == 1 | jstn == 3);
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
  if(ILOC == 1);
    print(fighyp,'-r250',strcat(plotfilehyp,plotext2),'-dpng');
  end;
  for ivo=1:NVMX;
    print(figmarg(ivo),'-painters','-r250',strcat(plotfilemarg,num2str(ivo),plotext2),'-dpng');
    print(figmargnod(ivo),'-painters','-r250',strcat(plotfilemargnod,num2str(ivo),plotext2),'-dpng');
  end;
  print(figMw,'-painters','-r250',strcat(plotfileMw,plotext2),'-dpng');
  print(figMwtot,'-painters','-r250',strcat(plotfileMwtot,plotext2),'-dpng');
  print(figk,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
  print(figbb2,'-painters','-r250',strcat(plotfilebb,plotext2),'-dpng');
  print(figbb3,'-painters','-r250',strcat(plotfilebb2,plotext2),'-dpng');
  print(figpctDC,'-painters','-r250',strcat(plotfilepctDC,plotext2),'-dpng');
  if(idatfit == 1);
    print(figdat,'-painters','-r250',strcat(plotfiledata,plotext2),'-dpng');
    print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext2),'-dpng');
    print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext2),'-dpng');
    print(figres,'-painters','-r250',strcat(plotfileres,plotext2),'-dpng');
    print(figreshist,'-painters','-r250',strcat(plotfilereshist,plotext2),'-dpng');
    if(ICOV ~= 4);    
    if(IAR == 1 | ICOV > 1);
      print(figrestot,'-painters','-r250',strcat(plotfilerestot,plotext2),'-dpng');
      print(figrestothist,'-painters','-r250',strcat(plotfilerestothist,plotext2),'-dpng');
      print(figacovr,'-painters','-r250',strcat(plotfileacovr,plotext2),'-dpng');
      print(figacovtot,'-painters','-r250',strcat(plotfileacovtot,plotext2),'-dpng');
      print(figressel,'-painters','-r250',strcat(plotfileressel,plotext3),'-dpng');
    end;end;
    if(IAR == 1);
      print(figarpred,'-painters','-r250',strcat(plotfilearpred,plotext2),'-dpng');
    end
  end;
  if(ieps == 1);
    print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-deps');
    if(ILOC == 1);
      print(fighyp,'-r250',strcat(plotfilehyp,plotext3),'-deps');
    end;
    for ivo=1:NVMX;
      print(figmarg(ivo),'-painters','-r250',strcat(plotfilemarg,num2str(ivo),plotext3),'-deps');
      print(figmargnod(ivo),'-painters','-r250',strcat(plotfilemargnod,num2str(ivo),plotext3),'-deps');
    end;
    print(figMw,'-painters','-r250',strcat(plotfileMw,plotext3),'-deps');
    print(figMwtot,'-painters','-r250',strcat(plotfileMwtot,plotext3),'-deps');
    print(figk,'-painters','-r250',strcat(plotfilek,plotext3),'-deps');
    print(figbb,'-painters','-r250',strcat(plotfilebb,plotext3),'-deps');
    print(figbb2,'-painters','-r250',strcat(plotfilebb2,plotext3),'-deps');
    print(figpctDC,'-painters','-r250',strcat(plotfilepctDC,plotext3),'-deps');
    if(idatfit == 1);
      print(figdat,'-painters','-r250',strcat(plotfiledata,plotext3),'-deps');
      print(figdatsel,'-painters','-r250',strcat(plotfiledatasel,plotext3),'-deps');
      print(figaxxsel,'-painters','-r250',strcat(plotfileaxxsel,plotext3),'-deps');
      print(figres,'-painters','-r250',strcat(plotfileres,plotext3),'-deps');
      print(figreshist,'-painters','-r250',strcat(plotfilereshist,plotext3),'-deps');
      if(ICOV ~= 4);
      if(IAR == 1 | ICOV > 1);
        print(figrestot,'-painters','-r250',strcat(plotfilerestot,plotext3),'-deps');
        print(figrestothist,'-painters','-r250',strcat(plotfilerestothist,plotext3),'-deps');
        print(figacovr,'-painters','-r250',strcat(plotfileacovr,plotext3),'-deps');
        print(figacovtot,'-painters','-r250',strcat(plotfileacovtot,plotext3),'-deps');
        print(figressel,'-painters','-r250',strcat(plotfileressel,plotext3),'-deps');
      end;end;
    end;
    if(IAR == 1);
      print(figarpred,'-painters','-r250',strcat(plotfilearpred,plotext3),'-deps');
    end
  end;
end

%return;
