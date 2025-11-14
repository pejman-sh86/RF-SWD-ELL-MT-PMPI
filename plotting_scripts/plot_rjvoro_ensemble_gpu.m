function []=plot_rjvoro_ensemble_gpu(filename);

%%
%% 3 Source 3 Freq
%%
Rmx = 4000;
Zmx = 300;
M1 = 1001;
M2 = 401;

%%
%% 2 Source 2 Freq
%%
%Rmx = 4000;
%Zmx = 250;
%M1 = 500;
%M2 = 250;

NVMX = 20;
%% Output files
filebase      = strrep(filename,'sample.mat','');
plotfilek     = strcat(filebase,'khist.');
plotfilelogL  = strcat(filebase,'logL.');
plotfilemap   = strcat(filebase,'map.');
plotfileens   = strcat(filebase,'ens.');
plotfilehpd   = strcat(filebase,'hpd.');
plotfiletru   = strcat(filebase,'tru.');
plotfileslice = strcat(filebase,'slice.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';
ctru = load('test_envc.txt');
rtru = load('test_envr.txt');
atru = load('test_enva.txt');

load(filename);
NPV = 5;
isave = 1;

NFPMX = NVMX*NPV;
m = A(:,5:NFPMX+4);
NSMP  = size(m,1)
NSMP2 = 4;
%NBURNIN = round(NSMP/4);
NBURNIN = 100;
x_evn = [0:1/(M2-1):1];
z_evn = [0:1/(M1-1):1];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
voro = zeros(NVMX,NPV,NSMP2);

%% Prior bounds:
minlim = [ 0.,  0., 1500., 1200.,0.00];
maxlim = [Zmx, Rmx, 1800., 1800.,0.02];

xplt = [30 100 200 300 380];

%% load bottom depth:
idx_bott = load('rjmh_pecan_bathy.dat');
%% load true environment
ctru = load('test_envc.txt');

thinstep = 1;
logLmin = min(A(:,1))-(max(A(:,1))-min(A(:,1)))/10;
logLmax = max(A(:,1))+(max(A(:,1))-min(A(:,1)))/10;

%% Discard burnin:
NSMP = NSMP-NBURNIN;
m = m(NBURNIN:end,:);
A = A(NBURNIN:end,:);
idx=randperm(NSMP);
m = m(idx(1:NSMP2),:);
k    = A(idx(1:NSMP2),4);
logL = A(idx(1:NSMP2),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figk=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
idx = find(A(:,end)==1);
stairs([1:length(idx)],A(idx,4),'k')
ylabel('No. nodes in partition');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(idx)])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(A(:,4),[0:20]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. interfaces in partition');
ylabel('Probability');
set(gca,'XLim',[0.5 NVMX+0.5]);
set(gca,'YLim',[0.   1.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%
figlogL=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14);

%for i=1:max(A(:,end));
i=1;
idx = find(A(:,end)==i);
   plot(A(idx(1:thinstep:end),1),'k');
   clear idx;
%   plot([1:thinstep:length(A(:,1))],A(1:thinstep:end,1),'k');
%end;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(A(:,1))])
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(A(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
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
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
idx = find(A(:,end)==1);
j = 1;
for i=1:5;
  subplot('Position',[loc(1,j) loc(2,j) spw sph]);hold on;box on;
  set(gca,'FontSize',14);
  plot(A(idx,end-10+i),'k');
  plot([1:length(idx)],0.2*ones(size(idx)),'--k');
  plot([1:length(idx)],0.3*ones(size(idx)),'--k');

  if(j==1 | j==4);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  if(j>3);xlabel('rjMCMC step');else;set(gca,'XTickLabel',[]);end;
  set(gca,'XLim',[0 length(idx)],'YLim',[0 0.5])
  j=j+1;if(j==3);j=j+1;end;
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
mapvoro = map(1:end,:);

x  = mapvoro(:,1);
z  = mapvoro(:,2);
xn = mapvoro(:,1)/Rmx;
zn = mapvoro(:,2)/Zmx;
c  = mapvoro(:,3);
r  = mapvoro(:,4);
a  = mapvoro(:,5);
%% Plot MAP
c_evmap(:,:) = ctru;
r_evmap(:,:) = rtru;
a_evmap(:,:) = atru;
for iz = min(idx_bott):M1;
  dz = z_evn(iz)-zn;
  dzsq = dz.*dz;
  for ir = 1:M2;
    if(idx_bott(ir) < iz);
      dx = x_evn(ir)-xn;
      d  = sqrt(dx.^2 + dzsq);
      [dmin,iv] = min(d);
      c_evmap(iz,ir) = mapvoro(iv,3);
      r_evmap(iz,ir) = mapvoro(iv,4);
      a_evmap(iz,ir) = mapvoro(iv,5);
    end;
  end;
end;

nx = 1;
ny = 3;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
 
figmap=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])

subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,c_evmap);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)])
cb = colorbar('peer',gca);set(get(cb,'ylabel'),'String', 'Velocity (m/s)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,r_evmap);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)])
cb = colorbar('peer',gca);set(get(cb,'ylabel'),'String', 'Density (kg/m^3)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,a_evmap);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse');
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)])
cb = colorbar('peer',gca);set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)');
xlabel('Range (km)');ylabel('Depth (m)');

clear A;
voro = zeros(NVMX,5,NSMP2);
t1 = tic;
for ivo = 1:NVMX;
  idx(ivo) = (ivo-1)*NPV+1;
end;
for ismp = 1:NSMP2;
  for ivo = 1:NVMX;
    voro(ivo,1:NPV,ismp) = m(ismp,idx(ivo):idx(ivo)+NPV-1);
  end;
end;
c_ens = zeros(M1,M2);c_ens2 = zeros(M1,M2,NSMP2);
r_ens = zeros(M1,M2);r_ens2 = zeros(M1,M2,NSMP2);
a_ens = zeros(M1,M2);a_ens2 = zeros(M1,M2,NSMP2);
for ismp = 1:NSMP2;
  if(rem(ismp,2) == 0);
    t2 = toc(t1);
    disp([num2str(ismp),',   ',num2str(NSMP2),',  time: ',num2str(t2)]);
    t1 = tic;
  end;
  xn = gpuArray.zeros(k(ismp));zn = gpuArray.zeros(k(ismp));
  x_evn_d = gpuArray.zeros(M2);z_evn_d = gpuArray.zeros(M1);
  c = gpuArray.zeros(k(ismp));r = gpuArray.zeros(k(ismp));a = gpuArray.zeros(k(ismp));
  c_ev = gpuArray.zeros(M1,M2);r_ev = gpuArray.zeros(M1,M2);a_ev = gpuArray.zeros(M1,M2);

  x_evn_d = x_evn; z_evn_d = z_evn;
  xn = voro(1:k(ismp),1,ismp)/Rmx;
  zn = voro(1:k(ismp),2,ismp)/Zmx;
  c = voro(1:k(ismp),3,ismp); c_ev(:,:) = ctru;
  r = voro(1:k(ismp),4,ismp); r_ev(:,:) = rtru;
  a = voro(1:k(ismp),5,ismp); a_ev(:,:) = atru;

  for ir = 1:M2;
    dx = x_evn_d(ir)-xn;
    dxsq = dx.*dx;
    for iz = idx_bott(ir):M1;
      dz = z_evn_d(iz)-zn;
      dzsq = dz.*dz;
      if(idx_bott(ir) < iz);
        d  = sqrt(dxsq + dzsq);
        [dmin,iv] = min(d);
        c_ev(iz,ir) = c(iv);
        r_ev(iz,ir) = r(iv);
        a_ev(iz,ir) = a(iv);
      end;
    end;
  end;
  c_ens2(:,:,ismp) = gather(c_ev);
  r_ens2(:,:,ismp) = gather(r_ev);
  a_ens2(:,:,ismp) = gather(a_ev);
end;
c_ens = mean(c_ens2,3);
r_ens = mean(r_ens2,3);
a_ens = mean(a_ens2,3);

%% Compute 95% HPD map.
c_hpdw = zeros(size(c_ens));
r_hpdw = zeros(size(r_ens));
a_hpdw = zeros(size(a_ens));
for ir = 1:M2;
  for iz = idx_bott(ir):M1;
    [nfc(iz,ir,:)] = hpd(c_ens2(iz,ir,:),100,95);
    c_hpdw(iz,ir) = abs(nfc(iz,ir,2)-nfc(iz,ir,1));
    [nfr(iz,ir,:)] = hpd(r_ens2(iz,ir,:),100,95);
    r_hpdw(iz,ir) = abs(nfr(iz,ir,2)-nfr(iz,ir,1));
    [nfa(iz,ir,:)] = hpd(a_ens2(iz,ir,:),100,95);
    a_hpdw(iz,ir) = abs(nfa(iz,ir,2)-nfa(iz,ir,1));
  end;
end;

nx = 1;
ny = 3;
if(length(xplt)>5);ny = 2;end;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

%%
%% Ensemble figure
%%
figens=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,c_ens);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,r_ens);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Density (kg/m^3)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
set(gca,'FontSize',14);
imagesc(x_ev,z_ev,a_ens);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)');
xlabel('Range (km)');ylabel('Depth (m)');

%%
%% HPD figure
%%
fighpd=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,c_hpdw);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 100])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', '95% HPD (m/s)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,r_hpdw);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[0 100])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', '95% HPD (kg/m^3)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,a_hpdw);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[0 .02])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', '95% HPD (dB/m)');
xlabel('Range (km)');ylabel('Depth (m)');

%%
%% True figure
%%
figtru=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,ctru);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(3) maxlim(3)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Velocity (m/s)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,rtru);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'XTickLabel',[],'TickDir','out')
colorbar;set(gca,'CLim',[minlim(4) maxlim(4)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Density (kg/m^3)');
xlabel('Range (km)');ylabel('Depth (m)');

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
imagesc(x_ev,z_ev,atru);shading flat;
set(gca,'XLim',[x_ev(1) x_ev(end)],'YLim',[z_ev(1) z_ev(end)],'YDir','reverse')
set(gca,'TickDir','out')
colorbar;set(gca,'CLim',[minlim(5) maxlim(5)])
cb = colorbar('peer',gca);
set(get(cb,'ylabel'),'String', 'Attenuation (dB/m)');
xlabel('Range (km)');ylabel('Depth (m)');

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
  set(gca,'FontSize',14);
  plot(c_ens(:,xplt(i)),z_ev,'-b');
  plot(ctru(:,xplt(i)),z_ev,'--r');
  plot(nfc(:,xplt(i),1),z_ev,'--k');
  plot(nfc(:,xplt(i),2),z_ev,'--k');
  set(gca,'YDir','reverse','YLim',[100 300],'XLim',[minlim(3) maxlim(3)]);
  xlabel('Velocity (m/s)');
  set(gca,'YTick',[0:50:250]);
  if(i == 1);
    ylabel('Depth (m)');
  else;
    set(gca,'YTickLabel',[]);
  end;
  set(gca,'XTick',[minlim(3):(maxlim(3)-minlim(3))/5:maxlim(3)]);
  text(1600,15,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
end;
for i=1:length(xplt);
  subplot('Position',[loc(1,i+length(xplt)) loc(2,i+length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14);
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
  text(1600,15,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
end;
for i=1:length(xplt);
  subplot('Position',[loc(1,i+2*length(xplt)) loc(2,i+2*length(xplt)) spw sph]);hold on;box on;
  set(gca,'FontSize',14);
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
  text(1600,15,['r=' num2str(round(x_ev(xplt(i)))) ' m'],'FontSize',12,'Color',[0,0,0]);
end;

if(isave == 1)
%  print(figk,'-painters','-r250',strcat(plotfilek,plotext3),'-depsc');
  print(figk,'-painters','-r250',strcat(plotfilek,plotext2),'-dpng');
%  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext3),'-depsc');
  print(figlogL,'-painters','-r250',strcat(plotfilelogL,plotext2),'-dpng');
%  print(figmap,'-painters','-r250',strcat(plotfilemap,plotext3),'-depsc');
  print(figmap,'-r250',strcat(plotfilemap,plotext2),'-dpng');
%  print(figens,'-painters','-r250',strcat(plotfileens,plotext3),'-depsc');
  print(figens,'-painters','-r250',strcat(plotfileens,plotext2),'-dpng');
%  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext3),'-depsc');
  print(fighpd,'-painters','-r250',strcat(plotfilehpd,plotext2),'-dpng');
%  print(figtru,'-painters','-r250',strcat(plotfiletru,plotext3),'-depsc');
  print(figtru,'-painters','-r250',strcat(plotfiletru,plotext2),'-dpng');
%  print(figslice,'-painters','-r250',strcat(plotfileslice,plotext3),'-depsc');
  print(figslice,'-painters','-r250',strcat(plotfileslice,plotext2),'-dpng');
end

return;
