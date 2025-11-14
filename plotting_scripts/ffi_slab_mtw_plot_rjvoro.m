%function []=ffi_plot_rjvoro(filename);

filename = 'bohol_data_sample.mat';
isave   = 1;
ieps    = 1;
idatfit = 1;
isyn    = 0;
inonstat= 1;
ialf    = 0;
smax = 90

alfmin = 59;
NSDGPS = 4; %% 3 sd for x, y, z, and 1 for alp (Okada "constant")

filebase    = strrep(filename,'_sample.mat','');
parfile     = strcat(filebase,'_parameter.dat');
edgefile    = strcat(filebase,'_fault_edge.txt');
M0file      = strcat(filebase,'_M0.dat');
datfile     = strcat(filebase,'.hdf5');
datfilecGPS = strcat(filebase,'_cGPS.hdf5');
covfile     = strcat(filebase,'_covmat.mat');

%%
%% Output plot files
%% 
plotfilesltot= strcat(filebase,'_sltot.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

[IMAP,ICOV,I_WP,I_cGPS,I_GPS,NVMX,NPV,NMISC,IVRUP,ILATLON,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile);

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc');
rake   = h5readatt(datfile,'/Sensitivity_kernel','rake_comp');
dhyp   = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NTW    = h5readatt(datfile,'/Sensitivity_kernel','N_time_win');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin  = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax  = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
mu     = h5read(datfile,'/Rigidity/mu');
sfgrid = h5read(datfile,'/Sensitivity_kernel/subfaultgrid');

NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
NTW  = cast(NTW,'like',1);

NSF = NRAN*NDEP;
NFPMX = NVMX*NPV*NTW;

%% Sub fault size in range and depth
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
%% Newton meter scaling for moment
Nmscale = double(mu)*1.e6*deltr*deltz;

%% Padding in lat (ypad) and lon (xpad):
xpad = deltz/lldistkm([mean(sfgrid(2,:)),1],[mean(sfgrid(2,:)),2])/2;
ypad = deltz/lldistkm([mean(sfgrid(2,:)),1],[mean(sfgrid(2,:))+1,1])/2;
%% Prior bounds:
minlim(1:2) = [min(sfgrid(2,:))-ypad,min(sfgrid(1,:))-xpad];
maxlim(1:2) = [max(sfgrid(2,:))+ypad,max(sfgrid(1,:))+xpad];

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

%% Mesh to grid data on:
[xq,yq] = meshgrid(minlim(1):.02:maxlim(1), minlim(2):.02:maxlim(2));
xxq = xq(1,:);
yyq = yq(:,1);
%% Create clipping mask
edge = load(edgefile);
edgefine(:,1) = interp1(edge(:,2),edge(:,1),xxq,'linear','extrap');
edgefine(:,2) = interp1(edge(:,4),edge(:,3),xxq,'linear','extrap');
edgefine(:,1) = edgefine(:,1) - xpad;
edgefine(:,2) = edgefine(:,2) + xpad;
clip = ones(size(xq));
NX = size(xxq,2); NY = size(yyq,1);
for ix=1:NX;
    for iy = 1:NY;  
        if(yq(iy,ix) < edgefine(ix,1) | yq(iy,ix) > edgefine(ix,2));
            clip(iy,ix) = nan;
        end;
    end
end

%% File extensions for output:
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

%%
%% Load sample file
%%
load(filename);AS = A;
clear A;
ms = AS(:,4+NTW:NFPMX+3+NTW);
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
  if(I_cGPS == 1);
    datcGPS = h5read(datfilecGPS,'/Observed_data/displacements');
  end;
  NTSMP = h5read(datfile,'/Observed_data/Ntraces');
  NDAT = length(dat);
end;

%%
%%  MAP model from file:
%%
map_sl1=load('bohol_data_sl1.txt');
map_sl2=load('bohol_data_sl2.txt');
map_sl = sqrt(map_sl1.^2+map_sl2.^2);
map_slb=reshape(map_sl,NSF,NTW);

if(ILATLON == 1);
  x_ev = sfgrid(2,:);
  z_ev = sfgrid(1,:);
else;
  x_ev = linspace(deltr/2,Rmx-deltr/2,NRAN);
  z_ev = linspace(deltz/2,Zmx-deltz/2,NDEP);
end;
minlon = min(z_ev);
maxlon = max(z_ev);
minlat = min(x_ev);
maxlat = max(x_ev);

minlimmisc = [ hyp_loc(1)-dhyp, hyp_loc(2)-dhyp, Vrmin];
maxlimmisc = [ hyp_loc(1)+dhyp, hyp_loc(2)+dhyp, Vrmax];

minlimar = -0.5;
maxlimar =  1.0;

thinstep = 1;
logLmin = min(AS(:,1))-(max(AS(:,1))-min(AS(:,1)))/10;
logLmax = max(AS(:,1))+(max(AS(:,1))-min(AS(:,1)))/10;

%% Discard burnin:
k    = AS(:,4:3+NTW);
logL = AS(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figk=figure;

for itw = 1:NTW;
  subplot(2,NTW,itw);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  stairs(k(:,itw),'k')
  ylabel('No. nodes in partition');
  xlabel('rjMCMC step');
  set(gca,'XLim',[0 length(k)])

  subplot(2,NTW,itw+NTW);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  [n,lim]=hist(k(:,itw),[1:1:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
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
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%
figlogL=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
plot(logL(1:thinstep:end),'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(logL,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
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
ML = .045; MR = .03; MB = .07; MT = .02;
SP = .02; PAD = 0; FNT = 12;

figaccept=figure;hold on;box on;
title('Acceptance Rate');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12])
j = 1;
for ipar=1:NPV;
  subaxis(ny,nx,ipar,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
  plot(AS(:,end-11:end-11+NPV-1),'k');
  plot([1:length(k)],0.2*ones(size(k)),'--k');
  plot([1:length(k)],0.3*ones(size(k)),'--k');

  if(j==1 | j==4);ylabel('Acceptance rate');else;set(gca,'YTickLabel',[]);end;
  if(j>3);xlabel('rjMCMC step');else;set(gca,'XTickLabel',[]);end;
  set(gca,'XLim',[0 length(k)],'YLim',[0 0.5])
  j=j+1;if(j==3);j=j+1;end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
sltot_hpd = zeros(NDEP,NRAN);
nfsltot = zeros(NDEP,NRAN,2);
disp('Starting HPDs');

%% Compute CIs:
for ir = 1:NRAN;
  if(rem(ir,10) == 0);disp(ir);end;
  for iz = 1:NDEP;
    %% 95% HPD
    [nfsltot(iz,ir,:)] = hpd2(sltot_hst2(iz,ir,:),bins(4,:),95);
    sltot_hpd(iz,ir) = abs(nfsltot(iz,ir,2)-nfsltot(iz,ir,1));
  end;
end;
disp('Done HPDs');

%% COLORMAP
%mycol = [ ones(1,3); summer(128)];
mycol = [ ones(1,3); flipud(hot(128))];colormap(mycol);
mycol(2:15,:)=[];

%%
%%
%%  MAP PLOTS total slip
%%
nx = 2;
ny = 1;
ML = .045; MR = .0; MB = .07; MT = .02;
SP = .0; PAD = 0; PR = 0.1; FNT = 14;

figmap = figure; hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 10])
subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'Paddingright',PR,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
axis equal;
ens_slip = sum(sltot_ens,1);
vq = griddata(z_ev,x_ev,ens_slip,yq,xq,'nearest');
vq = vq.*clip;
pcolor(yq(:,1),xq(1,:),vq');shading flat;
xlabel('East (degree)');
ylabel('North (degree)');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
for i=1:length(z_ev);
  plot(z_ev(i),x_ev(i),'.w','Markersize',7)
  plot(z_ev(i),x_ev(i),'.k','Markersize',5)
end;
plot([140,144.5],[37,37],'-k')
plot([140,144.5],[40,40],'-k')
plot([141,141],[36,40.5],'-k')
plot([144,144],[36,40.5],'-k')
colormap(mycol);
set(gca,'CLim',[0 smax],'FontSize',FNT)
cb = colorbar('peer',gca,'FontSize',FNT,'Location','southoutside');
set(get(cb,'ylabel'),'String', 'Slip (m)','FontSize',FNT);
set(gca,'XLim',[minlon-.4,maxlon+.4],'YLim',[minlat-.4,maxlat+.4]);

subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'Paddingright',PR,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
axis equal;
hpd_slip = sltot_hpd(:);
vq = griddata(z_ev,x_ev,hpd_slip,yq,xq,'nearest');
vq = vq.*clip;
pcolor(yq(:,1),xq(1,:),vq');shading flat;
xlabel('East (degree)');
ylabel('North (degree)');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
for i=1:length(z_ev);
  plot(z_ev(i),x_ev(i),'.w','Markersize',9)
  plot(z_ev(i),x_ev(i),'.k','Markersize',5)
end;
plot([140,144.5],[37,37],'-k')
plot([140,144.5],[40,40],'-k')
plot([141,141],[36,40.5],'-k')
plot([144,144],[36,40.5],'-k')
colormap(mycol);
%colormap('default');
set(gca,'CLim',[0 60],'FontSize',FNT)
cb = colorbar('peer',gca,'FontSize',FNT,'Location','southoutside');
set(get(cb,'ylabel'),'String', '95% CI (m)','FontSize',FNT);
set(gca,'XLim',[minlon-.4,maxlon+.4],'YLim',[minlat-.4,maxlat+.4]);
if(isave == 1)
%  print(figmap,'-painters','-r250',strcat(plotfilesltot,plotext2),'-dpng');
%  saveas(figmap,strcat(plotfilesltot,plotext1),'fig');
  if(ieps == 1);
    print(figmap,'-painters','-r250',strcat(plotfilesltot,plotext3),'-depsc');
  end;
end

figmap2=figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
imagesc(reshape(sum(map_slb,2),9,25));
xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');
colormap(mycol);
colorbar;set(gca,'CLim',[0 smax],'FontSize',FNT)
cb = colorbar('peer',gca,'FontSize',FNT);
set(get(cb,'ylabel'),'String', 'MAP slip (m)','FontSize',FNT);

%%
%% Ensemble mean figure slip components
%%
nx = 2;
ny = 1;
ML = .045; MR = .03; MB = .07; MT = .02;
SP = .02; PAD = 0; FNT = 14;

figens=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 14])
subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
ens_slip = sum(sl1_ens,1);
vq = griddata(z_ev,x_ev,ens_slip,yq,xq);
pcolor(yq(:,1),xq(1,:),vq');shading flat;
xlabel('East (degree)');
ylabel('North (degree)');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
for i=1:length(z_ev);
  plot(z_ev(i),x_ev(i),'.w','Markersize',9)
  plot(z_ev(i),x_ev(i),'.k','Markersize',5)
end;
colormap(mycol);
colorbar;set(gca,'CLim',[0 smax],'FontSize',FNT)
cb = colorbar('peer',gca,'FontSize',FNT);
set(get(cb,'ylabel'),'String', 'Mean slip 1 (m)','FontSize',FNT);
set(gca,'XLim',[minlon-.2,maxlon+.2],'YLim',[minlat-.2,maxlat+.2]);

subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on; 
ens_slip = sum(sl2_ens,1);
vq = griddata(z_ev,x_ev,ens_slip,yq,xq);
pcolor(yq(:,1),xq(1,:),vq');shading flat;
xlabel('East (degree)');
ylabel('North (degree)');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
for i=1:length(z_ev);
  plot(z_ev(i),x_ev(i),'.w','Markersize',9)
  plot(z_ev(i),x_ev(i),'.k','Markersize',5)
end;
colormap(mycol);
colorbar;set(gca,'CLim',[0 smax],'FontSize',FNT)
cb = colorbar('peer',gca,'FontSize',FNT);
set(get(cb,'ylabel'),'String', 'Mean slip 2 (m)','FontSize',FNT);
set(gca,'XLim',[minlon-.2,maxlon+.2],'YLim',[minlat-.2,maxlat+.2]);

NZ1 = 1; NZ2 = 9;
NR1 = 1; NR2 = 6;
ffi_slab_mtw_slip_uncertainty(filename,sltot_hst2,bins,NR1,NR2,NZ1,NZ2,ieps,isyn);
NR1 = 7; NR2 = 12;
ffi_slab_mtw_slip_uncertainty(filename,sltot_hst2,bins,NR1,NR2,NZ1,NZ2,ieps,isyn);
NR1 = 13; NR2 = 18;
ffi_slab_mtw_slip_uncertainty(filename,sltot_hst2,bins,NR1,NR2,NZ1,NZ2,ieps,isyn);
NR1 = 19; NR2 = 24;
ffi_slab_mtw_slip_uncertainty(filename,sltot_hst2,bins,NR1,NR2,NZ1,NZ2,ieps,isyn);

%return;
