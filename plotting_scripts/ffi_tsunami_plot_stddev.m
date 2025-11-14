function []=ffi_tsunami_plot_stddev(file);

filebase = strrep(file,'_sample.mat','')
datfile  = strcat(filebase,'.hdf5');
parfile  = strcat(filebase,'_parameter.dat');

%logLmin = 3800;
%logLmax = 4000;
logLmin = -1280;
logLmax = -1160;
kmin = 100;
kmax = 135;

load(file);
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = cast(h5readatt(datfile,'/Observed_data','Sample_rate'),'like',1.);
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc_Cart');
dhyp   = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin  = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax  = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
NTW    = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
NTSMP  = cast(h5read(datfile,'/Observed_data/Ntraces'),'like',1);
NTW  = cast(NTW,'like',1);
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
NSTN = cast(NSTN,'like',1);
NTSMP = cast(NTSMP,'like',1);

minlimmisc = [-dhyp,-dhyp, 1.4];
maxlimmisc = [ dhyp, dhyp, 1.45];

%gauge = [801,802,803,806,804,807,21418,21401,21413,5741,5742,5861,5862,1002,1006];
%time1 = [200, 200, 200, 200, 200, 200,1380,3000,3900, 600, 900, 150, 150, 150, 150];
%time2 = [3600,3600,3600,4600,3600,3600,3600,6800,6700,3600,3600,1200,2000,2000,2000];
%ylimmax = [ 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
%ylimmin = [-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];
gauge = [ 205, 801, 802, 803, 806, 804, 807,21418,21401,21413,5741,5742,5861,5862,1002,1006];
time1 = [ 100, 300, 300, 300, 300, 300, 300, 1380, 3000, 3900, 600, 900, 300, 300, 250, 250];
time2 = [3800,3600,2400,3600,4600,3600,3600, 3600, 6800, 6700,3600,3600,1200,2000,2000,2000];
ylimmax = [ 8, 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
ylimmin = [-2,-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];

figdat = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
nx = 4;
ny = 4;
xim = 0.01;
yim = 0.02;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
ML = .045;
MR = .03;
MB = .07;
MT = .02;
SP = .02;
PAD = 0;
FNT = 12;


jj = 1;
for istn = 1:NSTN;
  subaxis(ny,nx,istn,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;

  [a,b] = hist(A(:,24+istn));
  a=[0,a,0];
  b=[b(1),b,b(end)];
  [xx,yy]=stairs(b,a,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(b,a,'k');
  set(gca,'XLim',[0 0.75])
  
  ylim=get(gca,'ylim');
  xlim=get(gca,'xlim');
  text(xlim(1)+(xlim(2)-xlim(1))/12,ylim(2)-(ylim(2)-ylim(1))/12,num2str(gauge(istn)));
  
  if(istn == 1 | istn == nx+1 | istn == 2*nx+1 | istn == 3*nx+1 );
    ylabel('Probability');
  else;
    set(gca,'YTickLabel',[]);
  end;
  if(istn > nx*(ny-1) );
    xlabel('Standard deviation (m)');
  else;
    set(gca,'XTickLabel',[]);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Hypocentre plot
%%
[amap,imap]=max(A(:,1));
hyp = A(:,6:8);

hyp(:,1) = hyp(:,1)-hyp_loc(1);
hyp(:,2) = hyp(:,2)-hyp_loc(2);

maphyp = hyp(imap,:);
fighyp2 = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
%cloudPlot(hyp(:,1),hyp(:,2),[],false,[100,100])
%axis equal;
n = hist3(hyp(:,1:2),[100,100]); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(hyp(:,1)),max(hyp(:,1)),size(n,1)+1);
yb = linspace(min(hyp(:,2)),max(hyp(:,2)),size(n,1)+1);
h = imagesc(xb,yb,n1);
shading flat 
mycol = [ ones(1,3); jet(128)];
colormap(mycol);
axis('equal');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)

xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[minlimmisc(1) maxlimmisc(1)],'YLim',[minlimmisc(2) maxlimmisc(2)]);
set(gca,'TickDir','out');
clear n lim1 lim2;



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

n = hist3(hyp(:,1:2),[100,100]); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(hyp(:,1)),max(hyp(:,1)),size(n,1)+1);
yb = linspace(min(hyp(:,2)),max(hyp(:,2)),size(n,1)+1);
h = imagesc(xb,yb,n1);
shading flat 
mycol = [ ones(1,3); jet(128)];
colormap(mycol);
axis('equal');
plot(hyp_loc(1),hyp_loc(2),'pw','Markersize',16,'MarkerFaceColor','w')
plot(maphyp(1),maphyp(2),'*w','Markersize',16)

xlabel('Along-strike distance (km)');
ylabel('Along-dip distance (km)');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'XLim',[minlimmisc(1) maxlimmisc(1)],'YLim',[minlimmisc(2) maxlimmisc(2)]);
set(gca,'TickDir','out');
clear n lim1 lim2;

figrupvel = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.5 3])
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim1]=hist(1./hyp(:,3),[minlimmisc(3):(maxlimmisc(3)-minlimmisc(3))/300:maxlimmisc(3)]);
n = [0, n, 0];lim1 = [lim1(1) lim1 lim1(end)];
n = n/sum(n);
lim1 = lim1-(lim1(3)-lim1(2))/2.;
[xx,yy]=stairs(lim1,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim1,n,'k');
xlabel('Rupture velocity (km/s)');
ylabel('Probability (s/km)');
clear n lim1;
set(gca,'XLim',[minlimmisc(3) maxlimmisc(3)]);
set(gca,'YLim',[0 0.08],'TickDir','out');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
logL = A(:,1);
k_sl = A(:,2);

nx=4;
ny=1;
SP = 0.04;
SH = 0.1;
PD = 0;
PDB = 0.0;
PDL = 0.0;
ML = 0.08;
MR = 0.02;
MB = 0.2;
MT = 0.02;

figk_sl1=figure;
subaxis(ny,nx,1,'SpacingHorizontal',SH,'Padding',PD,'Paddingbottom',PDB,...
        'Paddingleft',PD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);hold on;box on;
hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1);
stairs(k_sl(:),'k');
ylabel('No. coefficients');
xlabel('rjMCMC step');
set(gca,'XLim',[0 size(A,1)],'XTick',[0:50000:1000000],'TickDir','out');
set(gca,'YLim',[kmin kmax]);
    
subaxis(ny,nx,2,'SpacingHorizontal',SH,'Padding',PD,'Paddingbottom',PDB,...
        'Paddingleft',PD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);hold on;box on;
hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(k_sl(:),[0:200]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
set(gca,'YTickLabel',[])
set(gca,'YLim',[kmin kmax],'TickDir','out');
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL plot
%%

subaxis(ny,nx,3,'SpacingHorizontal',SH,'Padding',PD,'Paddingbottom',PDB,...
        'Paddingleft',PD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);hold on;box on;
hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
plot(A(:,1),'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(A(:,1))],'XTick',[0:50000:1000000])
set(gca,'YLim',[logLmin logLmax],'TickDir','out')

subaxis(ny,nx,4,'SpacingHorizontal',SH,'Padding',PD,'Paddingbottom',PDB,...
        'Paddingleft',PD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);hold on;box on;
hold on;box on;
set(gca,'FontSize',14,'layer','top','LineWidth',1)
[n,lim]=hist(A(:,1),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
set(gca,'YTickLabel',[])
if(logLmin == logLmax);
   set(gca,'YLim',[0 2])
else;
   set(gca,'YLim',[logLmin logLmax],'TickDir','out')
end;
set(gca,'TickDir','out')


return;