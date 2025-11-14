function []=plot_rjvoro_ensemble(filename);

%%
%% 3 Source 3 Freq
%%
Rmx = 4000;
Zmx = 300;
%NDEP = 1001;
NDEP = 1001;
NRAN = 401;
NHIST = 150;

ISMOOTHW = 1;
NRSMOOTH = 9;
NZSMOOTH = 17;
%%
%% 2 Source 2 Freq
%%
%Rmx = 4000;
%Zmx = 250;
%NDEP = 500;
%NRAN = 250;

NVMX  = 20;
NVMXW = 20;
ISED   = 1;
IWATER = 0;
%% Output files
filebase       = strrep(filename,'sample.mat','');
samplew        = strcat(filebase,'samplew.mat');
plotfilek      = strcat(filebase,'khist.');
plotfilekw     = strcat(filebase,'khistw.');
plotfilelogL   = strcat(filebase,'logL.');
plotfilemap    = strcat(filebase,'map.');
plotfileens    = strcat(filebase,'ens.');
plotfileensw   = strcat(filebase,'ensw.');
plotfilehpd    = strcat(filebase,'hpd.');
plotfilehpd66  = strcat(filebase,'hpd66.');
plotfiletru    = strcat(filebase,'tru.');
plotfileslice  = strcat(filebase,'slice.');
plotfileslicew = strcat(filebase,'slicew.');
plotfileslicew2= strcat(filebase,'slicew2.');
plotfilenodes  = strcat(filebase,'nodes.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';
ctru = load('test_envc.txt');
rtru = load('test_envr.txt');
atru = load('test_enva.txt');

%% Downsample by factor 2 in depth to save compute time:
%ctru = ctru(1:2:end,:);
%rtru = rtru(1:2:end,:);
%atru = atru(1:2:end,:);

isave = 1;
NPV = 5;
NFPMX = NVMX*NPV;
NSMP2 = 2000;
if(ISED == 1);
  load(filename);AS = A;
  clear A;
  ms = AS(:,5:NFPMX+4);
  NSMP  = size(ms,1)
  %NBURNIN = round(NSMP/4);
  NBURNIN = 1;
  voro = zeros(NVMX,NPV,NSMP2);
end;


if(IWATER == 1);
  load(samplew);AW = A;
  clear A;
  NPVW = 3;
  NFPMXW = NVMXW*NPVW;
  mw = AW(:,5:NFPMXW+4);
  if(ISED == 0);
    NSMP  = size(mw,1)
    %NBURNIN = round(NSMP/4);
    NBURNIN = 1;
  end;
  vorow = zeros(NVMXW,NPVW,NSMP2);
  %% Prior bounds:
  minlimw = [ 0.,  0., 1470.];
  maxlimw = [Rmx, Zmx, 1550.];
end;
x_evn = [0:1/(NRAN-1):1];
z_evn = [0:1/(NDEP-1):1];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';

%% Prior bounds:
minlim = [ 0.,  0., 1500., 1200.,0.00];
maxlim = [Rmx, Zmx, 1800., 1800.,0.02];

xplt = [30 100 200 300 380];

%% Experiment geometry:
%sources = [0 50; 500 50; 1000 50; 2000 50; 3000 50];
sources = [4000 50; 3500 50; 3000 50; 2000 50; 1000 50];
dep1 = 10.25;
dep2 = 134.25;
arraydep = [dep1:(dep2-dep1)/3:dep2];
%arrayrange = (Rmx-20)*ones(4,1);
arrayrange = 20*ones(4,1);

%% load bottom depth:
idx_bott = load('rjmh_pecan_bathy.dat');
%idx_bott = round(idx_bott/2);

for ir=1:NRAN;
  bott(ir) = z_ev(idx_bott(ir));
end;
bott = [bott, bott(end)]+(z_ev(2)-z_ev(1))/2;
dx_ev = x_ev(2)-x_ev(1);
xbott = x_ev-dx_ev/2;
xbott = [xbott x_ev(end)+dx_ev];

thinstep = 1;
if(ISED == 1)
  logLmin = min(AS(:,1))-(max(AS(:,1))-min(AS(:,1)))/10;
  logLmax = max(AS(:,1))+(max(AS(:,1))-min(AS(:,1)))/10;
else
  logLmin = min(AW(:,1))-(max(AW(:,1))-min(AW(:,1)))/10;
  logLmax = max(AW(:,1))+(max(AW(:,1))-min(AW(:,1)))/10;
end

%% Discard burnin:
NSMP = NSMP-NBURNIN;
if(ISED == 1)
  ms = ms(NBURNIN:end,:);
  AS = AS(NBURNIN:end,:);
  idx=randperm(NSMP);
  ms = ms(idx(1:NSMP2),:);
  k    = AS(idx(1:NSMP2),4);
  logL = AS(idx(1:NSMP2),1);
end;
if(IWATER == 1);
  if(ISED == 0);idx=randperm(NSMP);end;
  mw = mw(NBURNIN:end,:);
  AW = AW(NBURNIN:end,:);
  mw = mw(idx(1:NSMP2),:);
  kw   = AW(idx(1:NSMP2),4);
  if(ISED == 0);logL = AW(idx(1:NSMP2),1);end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
if(ISED == 1)
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
[n,lim]=hist(AS(:,4),[0:20]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. nodes sediment partition');
ylabel('Probability');
set(gca,'XLim',[0.5 NVMX+0.5]);
set(gca,'YLim',[0.   1.0]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% KW PLOT
%%
if(IWATER == 1);
  figkw=figure;
  subplot(2,1,1);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  clear idx;
  idx = find(AW(:,end-1)==1);
  stairs([1:length(idx)],AW(idx,4),'k')
  ylabel('No. nodes water-colum partition');
  xlabel('rjMCMC step');
  set(gca,'XLim',[0 length(idx)])

  subplot(2,1,2);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  [n,lim]=hist(AW(:,4),[0:20]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim-0.5;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;
  xlabel('No. nodes water-column partition');
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

%for i=1:max(AS(:,end));
i=1;
  if(ISED == 1)
    idx = find(AS(:,end)==i);
    plot(AS(idx(1:thinstep:end),1),'k');
  else;
    idx = find(AW(:,end-1)==i);
    plot(AW(idx(1:thinstep:end),1),'k');
  end;
  clear idx;
%   plot([1:thinstep:length(AS(:,1))],AS(1:thinstep:end,1),'k');
%end;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(AS(:,1))])
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
if(ISED == 1);
  [n,lim]=hist(AS(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
else;
  [n,lim]=hist(AW(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
end;
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
if(ISED == 1)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
idx = find(AS(:,end)==1);
j = 1;
for i=1:NPV;
  subplot('Position',[loc(1,j) loc(2,j) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(AS(idx,end-10+i),'k');
  plot([1:length(idx)],0.2*ones(size(idx)),'--k');
  plot([1:length(idx)],0.3*ones(size(idx)),'--k');

  if(j==1 | j==4);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  if(j>3);xlabel('rjMCMC step');else;set(gca,'XTickLabel',[]);end;
  set(gca,'XLim',[0 length(idx)],'YLim',[0 0.5])
  j=j+1;if(j==3);j=j+1;end;
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Acceptance rate plot water-column
%%
if(IWATER == 1);
nx = 3;
ny = 1;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
figacceptw=figure;hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
idx = find(AW(:,end-1)==1);
for i=1:NPVW;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(AW(idx,end-8+i),'k');
  plot([1:length(idx)],0.2*ones(size(idx)),'--k');
  plot([1:length(idx)],0.3*ones(size(idx)),'--k');

  if(i==1);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  xlabel('rjMCMC step');
  set(gca,'XLim',[0 length(idx)],'YLim',[0 0.5])
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Plot and write MAP model to file:
%%
%[xmap,imap]=max(A(:,1))
%map = A(imap,5:NFPMX+4);
%kmap = A(imap,4);
%for ivo = 1:kmap;
%  mapvoro(ivo,1:5) = map((ivo-1)*NPV+1:ivo*NPV);
%end;
%save('rjmh_pecan_map.dat','kmap','-ascii')
%save('rjmh_pecan_map.dat','mapvoro','-ascii','-append')
%disp('Updated MAP file');
map=dlmread('rjmh_pecan_map.dat');
kmap = map(1,1);
mapvoro = map(2:end,:);

x  = mapvoro(:,1);
z  = mapvoro(:,2);
c_evmap(:,:) = ctru;
r_evmap(:,:) = rtru;
a_evmap(:,:) = atru;
if(ISED == 1);
  xn = mapvoro(:,1)/Rmx;
  zn = mapvoro(:,2)/Zmx;
  c  = mapvoro(:,3);
  r  = mapvoro(:,4);
  a  = mapvoro(:,5);
  %% Plot MAP
  for iz = min(idx_bott):NDEP;
    dz = z_evn(iz)-zn;
    for ir = 1:NRAN;
      if(idx_bott(ir) < iz);
        dx = x_evn(ir)-xn;
        d  = sqrt(dx.*dx + dz.*dz);
        [dmin,iv] = min(d);
        c_evmap(iz,ir) = mapvoro(iv,3);
        r_evmap(iz,ir) = mapvoro(iv,4);
        a_evmap(iz,ir) = mapvoro(iv,5);
      end;
    end;
  end;
end;
if(IWATER == 1);
  mapw  = dlmread('rjmh_pecan_mapw.dat');
  kmapw = mapw(1,1);
  mapvorow = mapw(2:end,:);
  xw  = mapvorow(:,1);
  zw  = mapvorow(:,2);
  xnw = mapvorow(:,1)/Rmx;
  znw = mapvorow(:,2)/Zmx;
  cw  = mapvorow(:,3);
  for iz = 1:max(idx_bott);
    dzw = z_evn(iz)-znw;
    for ir = 1:NRAN;
      if(idx_bott(ir) >= iz);
        dxw = x_evn(ir)-xnw;
        dw  = sqrt(dxw.*dxw + dzw.*dzw);
        [dmin,iv] = min(dw);
        c_evmap(iz,ir) = mapvorow(iv,3);
  end;end;end;
  if(ISMOOTHW == 1);[c_evmap]=smooth_cwater(c_evmap,NRAN,idx_bott,NRSMOOTH,NZSMOOTH);end;
end;

nx = 1;
ny = 4;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

%%
%%
%%  MAP PLOTS
%%
%%
figmap=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_evmap));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
for i=1:kmapw;
  plot(Rmx-mapw(i+1,1),mapw(i+1,2),'.w','Markersize',7)
  plot(Rmx-mapw(i+1,1),mapw(i+1,2),'.k','Markersize',5)
end;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
if(IWATER == 1);colorbar;set(gca,'CLim',[minlimw(3) maxlimw(3)],'FontSize',14);
else;colorbar;set(gca,'CLim',[minlim(3) maxlim(3)],'FontSize',14);end;
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_evmap));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
for i=1:kmap;
  plot(Rmx-map(i+1,1),map(i+1,2),'.w','Markersize',7)
  plot(Rmx-map(i+1,1),map(i+1,2),'.k','Markersize',5)
end;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(r_evmap));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
for i=1:kmap;
  plot(Rmx-map(i+1,1),map(i+1,2),'.w','Markersize',7)
  plot(Rmx-map(i+1,1),map(i+1,2),'.k','Markersize',5)
end;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'Density (kg/m^3)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');

subplot('Position',[loc(1,4) loc(2,4) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(a_evmap));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
for i=1:kmap;
  plot(Rmx-map(i+1,1),map(i+1,2),'.w','Markersize',7)
  plot(Rmx-map(i+1,1),map(i+1,2),'.k','Markersize',5)
end;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse');
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');

clear AS;
t1 = tic;
if(ISED == 1);
  voro = zeros(NVMX,5,NSMP2);
  nodes = zeros(NVMX*NSMP2+NVMXW*NSMP2,2);
  for ivo = 1:NVMX;
    idx(ivo) = (ivo-1)*NPV+1;
  end;
  idxk = [1:NSMP2];
  %idxk = find(k == 5);
  for ismp = 1:length(idxk);
    for ivo = 1:NVMX;
      voro(ivo,1:NPV,ismp) = ms(idxk(ismp),idx(ivo):idx(ivo)+NPV-1);
    end;
    nodes((ismp-1)*NVMX+1:ismp*NVMX,:) = voro(:,1:2,ismp);
  end;
  clear idxk;
end;
if(IWATER == 1);
  vorow = zeros(NVMXW,3,NSMP2);
  for ivo = 1:NVMXW;
    idxw(ivo) = (ivo-1)*NPVW+1;
  end;
  idxk = [1:NSMP2];
  %idxk = find(kw == 10);
  for ismp = 1:length(idxk);
    for ivo = 1:NVMXW;
      vorow(ivo,1:NPVW,ismp) = mw(idxk(ismp),idxw(ivo):idxw(ivo)+NPVW-1);
    end;
    nodes((ismp-1)*NVMXW+1+NVMXW*NSMP2:(ismp*NVMXW)+NVMXW*NSMP2,:) = vorow(:,1:2,ismp);
  end;
  clear idxk;
end;
nodes(find(nodes(:,1) == 0),:)=[];

%%
%% Nodal density
%%
fignodes=figure();hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5])
set(gca,'FontSize',14,'layer','top','LineWidth',1)
nodes(:,1) = Rmx - nodes(:,1);
[hdens]=cloudPlot(nodes(:,1),nodes(:,2),[0 4000 0 300],true,[201 501]);
%colormap( 1-hot );
colormap( 1-gray(256) );
stairs(xbott,fliplr(bott),'-k','LineWidth',1);
plot(sources(:,1),sources(:,2),'*k','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vk','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','k')
plot(arrayrange,arraydep,'--k')
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YLim',[0 4000],'YLim',[0 300]);
set(gca,'CLim',[0 1.5],'FontSize',14)
xlabel('Range (m)');ylabel('Depth (m)');
box on;

c_ens = zeros(NDEP,NRAN);%c_ens2 = zeros(NDEP,NRAN,NSMP2);
r_ens = zeros(NDEP,NRAN);%r_ens2 = zeros(NDEP,NRAN,NSMP2);
a_ens = zeros(NDEP,NRAN);%a_ens2 = zeros(NDEP,NRAN,NSMP2);

c_hst = zeros(NDEP*NHIST,NRAN);c_hst2 = zeros(NDEP,NRAN,NHIST);
r_hst = zeros(NDEP*NHIST,NRAN);r_hst2 = zeros(NDEP,NRAN,NHIST);
a_hst = zeros(NDEP*NHIST,NRAN);a_hst2 = zeros(NDEP,NRAN,NHIST);

%%
%% Load previously computed and converted ensembles and histograms
%%
load hists.mat
if(ISMOOTHW == 1);
  [c_ens]=smooth_cwater(c_ens,NRAN,idx_bott,NRSMOOTH,NZSMOOTH);
  for ihist = 1:size(c_hst2,3);
    [c_hst2(:,:,ihist)]=smooth_cwater(c_hst2(:,:,ihist),NRAN,idx_bott,NRSMOOTH,NZSMOOTH);
  end;
end;

%for ismp = 1:NSMP2;
%  if(rem(ismp,20) == 0);
%    t2 = toc(t1);
%    disp([num2str(ismp),',   ',num2str(NSMP2),',  time: ',num2str(t2)]);
%    t1 = tic;
%  end;
%  
%  c_ev(:,:) = ctru;r_ev(:,:) = rtru;a_ev(:,:) = atru;
%  if(ISED == 1);
%    c = zeros(k(ismp));r = zeros(k(ismp));a = zeros(k(ismp));
%    xn = zeros(k(ismp));zn = zeros(k(ismp));
%    xn = voro(1:k(ismp),1,ismp)/Rmx;
%    zn = voro(1:k(ismp),2,ismp)/Zmx;
%    c = voro(1:k(ismp),3,ismp); 
%    r = voro(1:k(ismp),4,ismp); 
%    a = voro(1:k(ismp),5,ismp); 
%    for ir = 1:NRAN;
%      dx = x_evn(ir)-xn;
%      dxsq = dx.*dx;
%      for iz = idx_bott(ir):NDEP;
%        if(idx_bott(ir) < iz);
%          dz = z_evn(iz)-zn;
%          dzsq = dz.*dz;
%          d  = sqrt(dxsq + dzsq);
%          [dmin,iv] = min(d);
%          c_ev(iz,ir) = c(iv);
%          r_ev(iz,ir) = r(iv);
%          a_ev(iz,ir) = a(iv);
%        end;
%      end;
%    end;
%  end;
%  if(IWATER == 1);
%    xnw = zeros(kw(ismp));znw = zeros(kw(ismp));
%    cw = zeros(kw(ismp));
%    xnw = vorow(1:kw(ismp),1,ismp)/Rmx;
%    znw = vorow(1:kw(ismp),2,ismp)/Zmx;
%    cw  = vorow(1:kw(ismp),3,ismp);
%    for ir = 1:NRAN;
%      dxw = x_evn(ir)-xnw;
%      for iz = 1:idx_bott(ir);
%        if(idx_bott(ir) >= iz);
%          dzw = z_evn(iz)-znw;
%          dw  = sqrt(dxw.*dxw + dzw.*dzw);
%          [dmin,iv] = min(dw);
%          c_ev(iz,ir) = cw(iv);
%        end;
%      end;
%    end;
%    if(ISMOOTHW == 1);[c_ev]=smooth_cwater(c_ev,NRAN,idx_bott,NRSMOOTH,NZSMOOTH);end;
%  end;
%  c_ens = c_ens+c_ev; c_ens2(:,:,ismp) = c_ev;
%  r_ens = r_ens+r_ev; r_ens2(:,:,ismp) = r_ev;
%  a_ens = a_ens+a_ev; a_ens2(:,:,ismp) = a_ev;
%end;
%c_ens = c_ens/NSMP2;
%r_ens = r_ens/NSMP2;
%a_ens = a_ens/NSMP2;

%%
%% Compute 95% HPD map.
%%
c_hpdw = zeros(size(c_ens));
r_hpdw = zeros(size(r_ens));
a_hpdw = zeros(size(a_ens));
c_hpdw66 = zeros(size(c_ens));
r_hpdw66 = zeros(size(r_ens));
a_hpdw66 = zeros(size(a_ens));
%for ir = 1:NRAN;
%%  for iz = idx_bott(ir):NDEP;
%  for iz = 1:NDEP;
%    [nfc(iz,ir,:)] = hpd(c_ens2(iz,ir,:),100,95);
%    c_hpdw(iz,ir) = abs(nfc(iz,ir,2)-nfc(iz,ir,1));
%    [nfr(iz,ir,:)] = hpd(r_ens2(iz,ir,:),100,95);
%    r_hpdw(iz,ir) = abs(nfr(iz,ir,2)-nfr(iz,ir,1));
%    [nfa(iz,ir,:)] = hpd(a_ens2(iz,ir,:),100,95);
%    a_hpdw(iz,ir) = abs(nfa(iz,ir,2)-nfa(iz,ir,1));
%  end;
%end;
disp('Starting HPDs');
for ir = 1:NRAN;
  if(rem(ir,50) == 0);disp(ir);end;
  for iz = 1:NDEP;
    %% 95% HPD
    [nfc(iz,ir,:)] = hpd2(c_hst2(iz,ir,:),bins(1,:),95);
    c_hpdw(iz,ir) = abs(nfc(iz,ir,2)-nfc(iz,ir,1));
    [nfr(iz,ir,:)] = hpd2(r_hst2(iz,ir,:),bins(2,:),95);
    r_hpdw(iz,ir) = abs(nfr(iz,ir,2)-nfr(iz,ir,1));
    [nfa(iz,ir,:)] = hpd2(a_hst2(iz,ir,:),bins(3,:),95);
    a_hpdw(iz,ir) = abs(nfa(iz,ir,2)-nfa(iz,ir,1));
    %% 66% HPD
    [nfc66(iz,ir,:)] = hpd2(c_hst2(iz,ir,:),bins(1,:),66);
    c_hpdw66(iz,ir) = abs(nfc66(iz,ir,2)-nfc66(iz,ir,1));
    [nfr66(iz,ir,:)] = hpd2(r_hst2(iz,ir,:),bins(2,:),66);
    r_hpdw66(iz,ir) = abs(nfr66(iz,ir,2)-nfr66(iz,ir,1));
    [nfa66(iz,ir,:)] = hpd2(a_hst2(iz,ir,:),bins(3,:),66);
    a_hpdw66(iz,ir) = abs(nfa66(iz,ir,2)-nfa66(iz,ir,1));
  end;
end;
disp('Done HPDs');

nx = 1;
ny = 3;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

nx = 1;
ny = 4;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc2,spw2,sph2] = get_loc(nx,ny,xim,yim,xymarg);

%%
%% Ensemble mean figure
%%
figens=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_ens));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out','XTickLabel',[]);
if(IWATER == 1);colorbar;set(gca,'CLim',[minlimw(3) maxlimw(3)],'FontSize',14);
else;colorbar;set(gca,'CLim',[minlim(3) maxlim(3)],'FontSize',14);end;
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_ens));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,3) loc2(2,3) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(r_ens));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Density (kg/m^3)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,4) loc2(2,4) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(a_ens));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

%%
%% Water ensemble mean & true figure
%%
if(IWATER == 1);
figensw=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_ens));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out','XTickLabel',[]);
colorbar;set(gca,'CLim',[minlimw(3) maxlimw(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(ctru));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlimw(3) maxlimw(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;
end;

%%
%% HPD figure
%%
fighpd=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_hpdw));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 300],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% HPD (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(r_hpdw));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 500],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% HPD (kg/m^3)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(a_hpdw));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[0 .022],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', '95% HPD (dB/m)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

%%
%% 66% HPD figure
%%
fighpd66=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_hpdw66));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out','XTickLabel',[]);
colorbar;set(gca,'CLim',[0 10],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Std. dev. (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(c_hpdw66));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 150],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Std. dev. (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,3) loc2(2,3) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(r_hpdw66));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 300],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Std. dev. (kg/m^3)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,4) loc2(2,4) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(a_hpdw66));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[0 .015],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Std. dev. (dB/m)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

%%
%% True figure
%%
figtru=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])

subplot('Position',[loc2(1,1) loc2(2,1) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(ctru));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlimw(3) maxlimw(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,2) loc2(2,2) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(ctru));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,3) loc2(2,3) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(rtru));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Density (kg/m^3)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

subplot('Position',[loc2(1,4) loc2(2,4) spw2 sph2]);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
imagesc(x_ev,z_ev,fliplr(atru));shading flat;
stairs(xbott,fliplr(bott),'-w','LineWidth',1);
plot(sources(:,1),sources(:,2),'*w','MarkerSize',8,'LineWidth',1)
plot(arrayrange,arraydep,'vw','MarkerSize',8,'LineWidth',1,'MarkerFaceColor','w')
plot(arrayrange,arraydep,'--w')
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)],'FontSize',14)
cb = colorbar('peer',gca,'FontSize',14);
set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)','FontSize',14);
xlabel('Range (m)');ylabel('Depth (m)');
box on;

%save tmp.mat x_ev z_ev c_ens ctru xplt;

figslice=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
nx = length(xplt);
ny = 3;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.15/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(c_ens(:,xplt(i)),z_ev,'-b');
  plot(ctru(:,xplt(i)),z_ev,'--r');
  plot(nfc(:,xplt(i),1),z_ev,'--k');
  plot(nfc(:,xplt(i),2),z_ev,'--k');
  set(gca,'YDir','reverse','YLim',[100 300],'XLim',[minlim(3) maxlim(3)]);
  xlabel('Velocity (m/s)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[minlim(3):(maxlim(3)-minlim(3))/5:maxlim(3)]);
  text(1700,110,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
for i=1:length(xplt);
  subplot('Position',[loc(1,i+length(xplt)) loc(2,i+length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(r_ens(:,xplt(i)),z_ev,'-b');
  plot(rtru(:,xplt(i)),z_ev,'--r');
  plot(nfr(:,xplt(i),1),z_ev,'--k');
  plot(nfr(:,xplt(i),2),z_ev,'--k');
  set(gca,'YDir','reverse','YLim',[100 300],'XLim',[minlim(4) maxlim(4)]);
  xlabel('Density (kg/m^3)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[minlim(4):(maxlim(4)-minlim(4))/5:maxlim(4)]);
box on;
end;
for i=1:length(xplt);
  subplot('Position',[loc(1,i+2*length(xplt)) loc(2,i+2*length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(a_ens(:,xplt(i)),z_ev,'-b');
  plot(atru(:,xplt(i)),z_ev,'--r');
  plot(nfa(:,xplt(i),1),z_ev,'--k');
  plot(nfa(:,xplt(i),2),z_ev,'--k');
  set(gca,'YDir','reverse','YLim',[100 300],'XLim',[minlim(5) maxlim(5)]);
  xlabel('Attenuation (dB/m)');
  set(gca,'YTick',[0:50:350]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[minlim(5):(maxlim(5)-minlim(5))/5:maxlim(5)]);
box on;
end;

%%
%%  SLICE figure including water
%%
figslicew=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 10])
if(IWATER == 1);
nx = length(xplt);
ny = 1;
xim = 0.01;
yim = 0.15/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(c_ens(:,xplt(i)),z_ev,'-b');
  plot(ctru(:,xplt(i)),z_ev,'--r');
  plot(nfc(:,xplt(i),1),z_ev,'--k');
  plot(nfc(:,xplt(i),2),z_ev,'--k');
    set(gca,'YDir','reverse','YLim',[minlim(2) maxlim(2)],'XLim',[minlimw(3) maxlim(3)]);
  xlabel('Velocity (m/s)');
  set(gca,'YTick',[0:50:250]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[minlim(3):(maxlim(3)-minlim(3))/5:maxlim(3)]);
  text(1700,15,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
end;
%%
%%  SLICE figure including water
%%
figslicew2=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 8])
if(IWATER == 1);
nx = length(xplt);
ny = 1;
xim = 0.01;
yim = 0.15/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:length(xplt);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  plot(c_ens(:,xplt(i)),z_ev,'-b');
  plot(ctru(:,xplt(i)),z_ev,'--r');
  plot(nfc(:,xplt(i),1),z_ev,'--k');
  plot(nfc(:,xplt(i),2),z_ev,'--k');
  set(gca,'YDir','reverse','YLim',[minlim(2) 210],'XLim',[minlimw(3) maxlimw(3)]);
  xlabel('Velocity (m/s)');
  set(gca,'YTick',[0:50:250]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  %set(gca,'XTick',[minlim(3):(maxlim(3)-minlim(3))/5:maxlim(3)]);
  set(gca,'XTick',[1480:20:1560]);
  text(1700,15,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
box on;
end;
end;

if(isave == 1)
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-depsc');
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext2),'-dpng');
  print(figmap,'-painters','-r250',strcat(plotfilemap,plotext3),'-depsc');
  print(figmap,'-r250',strcat(plotfilemap,plotext2),'-dpng');
  print(figens,'-painters','-r250',strcat(plotfileens,plotext3),'-depsc');
  print(figens,'-painters','-r250',strcat(plotfileens,plotext2),'-dpng');
  print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext3),'-depsc');
  print(fignodes,'-painters','-r250',strcat(plotfilenodes,plotext2),'-dpng');
  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext3),'-depsc');
  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext2),'-dpng');
  print(fighpd66,'-painters','-r250',strcat(plotfilehpd66,plotext3),'-depsc');
  print(fighpd66,'-painters','-r250',strcat(plotfilehpd66,plotext2),'-dpng');
  print(figtru,'-painters','-r250',strcat(plotfiletru,plotext3),'-depsc');
  print(figtru,'-painters','-r250',strcat(plotfiletru,plotext2),'-dpng');
  print(figslice,'-painters','-r250',strcat(plotfileslice,plotext3),'-depsc');
  print(figslice,'-painters','-r250',strcat(plotfileslice,plotext2),'-dpng');
  if(ISED == 1);
    print(figk,'-painters','-r250',strcat(plotfilek,plotext3),'-depsc');
    print(figk,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
  end;
  if(IWATER == 1);
    print(figkw,'-painters','-r250',strcat(plotfilekw,plotext3),'-depsc');
    print(figkw,'-painters','-r250',strcat(plotfilekw,plotext2),'-dpng');
    print(figensw,'-painters','-r250',strcat(plotfileensw,plotext3),'-depsc');
    print(figensw,'-painters','-r250',strcat(plotfileensw,plotext2),'-dpng');
    print(figslicew,'-painters','-r250',strcat(plotfileslicew,plotext3),'-depsc');
    print(figslicew,'-painters','-r250',strcat(plotfileslicew,plotext2),'-dpng');
    print(figslicew2,'-painters','-r250',strcat(plotfileslicew2,plotext3),'-depsc');
    print(figslicew2,'-painters','-r250',strcat(plotfileslicew2,plotext2),'-dpng');
  end;
end

return;
%=============================================================================
function [c_ev]=smooth_cwater(c_ev,NRAN,idx_bott,NRSMOOTH,NZSMOOTH)
%
% Brute force nearest neighbour interpolation
%
%=============================================================================

ctmp = c_ev;
halfr = (NRSMOOTH-1)/2;
halfz = (NZSMOOTH-1)/2;
for ir = 1:NRAN;
  lor = ir - halfr;
  hir = ir + halfr;
  if(lor <= 0); lor = 1;end;
  if(hir >= NRAN); hir = NRAN;end;
  for iz = 1:idx_bott(ir);
    loz = iz - halfz;
    hiz = iz + halfz;
    if(loz <= 0); loz = 1;end;
    if(hiz >= idx_bott(lor)); hiz = idx_bott(lor);end
    if(hiz >= idx_bott(hir)); hiz = idx_bott(hir);end;
    c_ev(iz,ir) = sum(sum(ctmp(loz:hiz,lor:hir),2),1)/((hir-lor+1)*(hiz-loz+1));
  end;
end;
return;
