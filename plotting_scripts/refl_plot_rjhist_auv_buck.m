function [] = plot_rjhist_auv(filename);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 1; %% Normalize profile marginals line by line
Tstar   = 1.;
isyn    = 0;
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
isave   = 1;
idatfit = 1;
idathist= 0;
iboudreau= 0;
ispher  = 1;
iar     = 0;
isd     = 0;
ibucking= 3;
ivgsnosh= 0;
if(ibucking == 3 & ivgsnosh == 0);
  NPL = 5;
else;
  NPL = 4;
end;
if(ivgsnosh == 1)
  NMISC = 5;
else;
  NMISC = 4;
end;

filebase    = strrep(filename,'sample.mat','');
datafile    = strrep(filename,'_sample.mat','.txt');
repfile     = strrep(filename,'_sample.mat','_rep.dat');
mapfile     = strrep(filename,'_sample.mat','_map.dat');
profilefile = strcat(filebase,'profile.mat');
resfile     = strcat(filebase,'res.dat');
resarfile   = strcat(filebase,'resar.dat');
repensfile  = strcat(filebase,'rep_ens.dat');

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';
plotfile1   = strcat(filebase,'chains_logl.');
plotfile2   = strcat(filebase,'ar_par.');
plotfile3   = strcat(filebase,'sigma.');
plotfile33  = strcat(filebase,'misc.');
plotfile4   = strcat(filebase,'number_interfaces.');
plotfile5   = strcat(filebase,'logL.');
plotfile6   = strcat(filebase,'transdim_marg_prof.');
plotfile66  = strcat(filebase,'transdim_map_prof.');
plotfile7   = strcat(filebase,'datafit.');
plotfile77  = strcat(filebase,'datafit2.');
plotfile8   = strcat(filebase,'axx.');
plotfile9   = strcat(filebase,'reshist.');
corefile    = 'core.mat';

NAVEF = 5;
IGA = 2
IEPS = 0
frbw = 1.0/15.0;
%frbw = 1.0/20.0;
FBW  = 37.5;
icore   = 0;
ibl     = 0;
nffile = 'nf_band3.mat';
hmin = 0.05;
logLmin = 100;
logLmax = 150;
%logLmin = 600;
%logLmax = 700;
%logLmin = 780;
%logLmax = 880;

%% Site 1  8 m packet
%bands = [400 504  635  800 1008 1270 1600 2016];

%% Site 02  4 m packet
%bands = [300.,400.,504.,635.,800.,1008.,1270,1600];
bands = [315.,400.,500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.];

%% Site 13  8 m packet
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000 6300];

%% NJ Site 05  12 m packet
%bands = [400 504  635  800 1008 1270 1600 2016 2540 3200];

%% Clutter site (Panama City)
%bands = [1601.4, 2017.6, 2542, 3202.8, 4035.2, 5084.1, 6405.5];
%bands = [1677.4, 2017.6, 2542, 3202.8, 4035.2, 5084.1];
%bands = [2017.6, 2542, 3202.8, 4035.2, 5084.1, 5800.];
%bands = [2212.96, 2542.02, 3202.75, 4035.21, 5084.05, 5840.04 ];

%% Site 16  8 m packet
%bands = [504  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  600  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000 6300];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032];
%bands = [504  635  800 1008 1270 1600 2016 2540];

%% Site 07
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250.,1600.,2000.,2500.,3150.];
%bands = [ 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000.];
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000.,6300.,8000.,10000.];
%% Site 16  8 m packet
%bands = [500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.,4000.,5000.];
NBAND = length(bands);
NSD = NBAND;
%% AUV:
%pmin = [1420 1.2 0.0]';
%pmax = [1800 2.4 1.0]';
%% Site07:
%pmin = [1550 1.4 0.0]';
%pmax = [1800 2.3 1.0]';
%% Bucking simulation:
if(ibucking == 1)
  pmin = [0.25 1.0 log(50.00)]';
  pmax = [0.90 2.1 log(2.e6)]';
elseif(ibucking == 2)
  pmin = [0.25 log(2.e7) 0.05]';
  pmax = [0.90 log(8.e8) 0.15]';
elseif(ibucking == 3)
  if(ivgsnosh == 0)
    pmin = [0.20 log(1.e6) log(0.001) log(8.95e-5)]';
    pmax = [0.91 log(5.e9) log(0.147) log(0.4)]';
  else;
    pmin = [0.20 log(2.e7) -10]';
    pmax = [0.90 log(8.e8) -3]';
  end;
elseif(ibucking == 4)
  pmin = [0.20 1.0 log(1.e-17)]';
  pmax = [0.91 2.0 log(1.e-9)]';
end;
if(ibucking > 0)
  miscmin = [3.55e10, 2.34e9,2710.,1022.,0.02]';
  miscmax = [3.65e10, 2.38e9,2740.,1026.,0.12]';
end;
%   NDAVE   = ceil(NPROF/10);
%NDAVE   = 500;
thinstep = 1;

%% Charles' VGS simulations:
  mtru = [0.1500      0.8000     18.0000     -6.6036     -4.1730,...
          0.6000      0.7000     19.0000     -6.6036     -1.0370,...
          1.4000      0.6000     20.0000     -6.6036     -6.2030,...
          4.4450      0.5000     21.0000     -6.6036     -5.2100,...
                      0.4000     22.0000     -3.9380     -5.2100];
  mtrumisc = [35612723556.66, 2355895918.32, 2724.76, 1029.88];
load(filename);
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NBAND);
   armax = 1.*ones(1,NBAND);
else
   order = 1;
end;
%if(ispher == 0)
%   tmp = load(datafile);
%   dobs   = tmp(1:length(bands),:)';
%   angobs = tmp(length(bands)+1,:);
%   rex = tmp(length(bands)+2:end,:)';
%else
%hmax    = 1.5;
%if(idatfit == 1)
   tmp = dlmread(datafile);
   z_t    = tmp(1,1);
   cw     = tmp(2,1);
   rw     = tmp(3,1);
   hmax   = tmp(4,1)+.2;
%   hmax   = 1.5;
   dobs   = tmp(5:length(bands)+4,:)';
   angobs = tmp(length(bands)+5,:);
   NANG   = length(angobs);
   rex    = tmp(length(bands)+6:end-1,:)';
%end

if(icore == 1)
    load(corefile);
end;

NPROF = length(A);
if(NPROF > 100000);
  BURNIN = 20000;
else;
  BURNIN = ceil(NPROF/5);
end;
BURNIN = 1;
A = A(BURNIN:end,:);
%A=A(1:2000,:);
%idx=find(A(:,112)==1);
%A=A(idx,:);

NPROF = length(A);

%%
%% Random permutation of sample is only applied to data fit computation
%%
idxran = randperm(length(A));
%A = A(idxran,:);
disp('Done lodaing data.');

k = A(:,4);
%% Compute prior volume
logP = zeros(size(k));
for i=1:length(k);
  vol = 1.;
  for ik=1:k(i);
    vol = vol * (hmax-(ik-1.)*hmin)*0.489*(pmax(end)-pmin(end));
  end;
  vol = vol * 0.489 * (pmax(end)-pmin(end));  %% This is half-space volume
  logP(i) = -log(vol);
end;

NFP = (k*4)+3;
m = A(:,5:end-5);

%%
%% Save MAP for max(P(k)) to file for starting new inversion
%%
kmin = min(k);
kmax = max(k);
if(kmax > 0);
  for i = 1:kmax-kmin+1;
    ik = kmin+i-1;
    idx = find(A(:,4) == ik);
    nmod_k(i) = length(idx);
    [a,b] = max(A(idx,1));
    mapk(i).par=A(idx(b),5:4+A(idx(b),4)*4+3);
    clear idx;
  end
  [a,b] = max(nmod_k);
  kpeak = b+kmin - 1
  idx = find(A(:,4) == kpeak);
  [a,b] = max(A(idx,1));map=A(idx(b),4:end-6);
  clear idx;
else
  [a,b] = max(A(:,1));map=A(b,4:end-6);
end;
save(mapfile,'map','-ascii');

if(iar == 1)
   (size(A,2)-6-order*NBAND-NMISC+1:size(A,2)-NMISC-6)
   alpha = A(:,end-6-(order*NBAND)-NMISC+1:end-NMISC-6);
end
if(ibucking > 0)
   (size(A,2)-6-NMISC-order*NBAND+1:size(A,2)-6-NMISC)
   bulkrho = A(:,end-6-(NMISC)+1:end-6);
end

logL = A(:,1);
if(imap == 1);
   [logLmap,jmap] = max(logL);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

if(kmax > 0);
  for i = 1:size(A,1);
    idxh = [0:k(i)-1]*4+1;
    h(i) = m(i,idxh(end));
    clear idxh;
  end;
end;
if(isd == 1)
   (size(A,2)-6-(order*NBAND)-NSD-NMISC+1:size(A,2)-6-(order*NBAND)-NMISC)
   sd(:,1:NSD) = A(:,end-6-(order*NBAND)-NSD-NMISC+1:end-6-(order*NBAND)-NMISC);
   disp('Done getting sigma.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig1=figure;
NTH = max(A(:,end));
NPT = max(A(:,end-1));
nx = ceil(sqrt(NTH*NPT));
ny = nx;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ith = 1;
ipt = 1;
jj=0;
for i=1:NTH;
for j=1:NPT;
   jj=jj+1;
   h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
   hold on; box off;
   set(gca,'FontSize',14);

   clear idx1 idx2;
   if(find(A(:,end)==i));
   idx1 = find(A(:,end)==i);
   idx2 = find(A(idx1,end-1)==j);
   [AX,H1,H2] = plotyy([1:length(idx2)],...
                A(idx1(idx2),1),[1:length(idx2)],...
                A(idx1(idx2),4));
   set(AX(1),'YTick',[0:40:1000]);
   if(rem(i-1,nx)==0);
      set(get(AX(1),'Ylabel'),'String','logL') 
      set(AX(1),'YTickLabel',[0:40:1000]);
   else;
      set(AX(1),'YTickLabel',[]);
   end;
   set(AX(2),'YTick',[0:2:20]);
   if(rem(i,nx)==0);
      set(get(AX(2),'Ylabel'),'String','No. interfaces') 
      set(AX(2),'YTickLabel',[0:2:20]);
   else;
      set(AX(2),'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(AX(2),'XTickLabel',[],'YLim',[kmin-1 kmax+1]);
   set(gca,'YLim',[logLmin logLmax]);
   end;
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Posterior CHAIN PLOTS
%%
fig11   = figure;
NTH = max(A(:,end));
NPT = max(A(:,end-1));
nx = ceil(sqrt(NTH*NPT));
ny = nx;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ith = 1;
ipt = 1;
jj=0;
for i=1:NTH;
for j=1:NPT;
   jj=jj+1;
   h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
   hold on; box off;
   set(gca,'FontSize',14);

   clear idx1 idx2;
   if(find(A(:,end)==i));
   idx1 = find(A(:,end)==i);
   idx2 = find(A(idx1,end-1)==j);
   plot([1:length(idx2)],A(idx1(idx2),1)+logP(idx1(idx2)));
   set(gca,'YTick',[0:40:1000]);
   if(rem(i-1,nx)==0);
      set(gca,'YTickLabel',[0:40:1000]);
   else;
      set(gca,'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(gca,'YLim',[logLmin+max(logP) logLmax+min(logP)]);
   end;
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
if(isd == 1)
   fig3=figure;
   nx = NSD;
   ny = 2;
   xim = 0.01;
   yim = 0.06;
   xymarg = [0.1 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NSD;
      h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14);
      h2 = subplot('Position',[loc(1,i+NSD) loc(2,i+NSD) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14);
      subplot(h1);hold on;box on;
      plot([1:thinstep:length(sd(:,i))],sd(1:thinstep:end,i),'k');
      if(i == 1);ylabel('Data error standard deviation');end;
      xlabel('rjMCMC step');
      set(gca,'XLim',[0 length(sd(:,i))])
      set(gca,'YLim',[0.0 0.10],'YTick',[0.0 0.02 0.04 0.06])
      if(i > 1);set(gca,'YTickLabel',[]);end;

      subplot(h2);set(gca,'Layer','top');hold on;
      [n,lim]=hist(sd(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      [xx,yy]=stairs(n,lim,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(n,lim,'k');
      clear n lim;
      if(i == 1);ylabel('Data error standard deviation');end;
%      xlabel('Probability');
      set(gca,'YLim',[0.0 0.10],'YTick',[0.0 0.02 0.04 0.06])
      if(i > 1);set(gca,'YTickLabel',[]);end;
      box on;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Bulk Moduli and Density (pore fluid/grain material) PLOT
%%
if(ibucking > 0)
   fig33=figure;
   nx = ceil(NMISC/2);
   ny = 2;
   xim = 0.01;
   yim = 0.10;
   xymarg = [0.1 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NMISC;
      h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14,'Layer','top');
      miscmead(i) = median(bulkrho(:,i));
      [n,lim]=hist(bulkrho(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      [xx,yy]=stairs(lim,n,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(lim,n,'k');
      if(isyn == 1);
        plot([mtrumisc(i) mtrumisc(i)],[0 .04],'-w');
        plot([mtrumisc(i) mtrumisc(i)],[0 .04],'--k');
      end;
      clear n lim;
      set(gca,'XLim',[miscmin(i) miscmax(i)])
      if(i == 1);xlabel('Ks (mineral grains) [Pa]');end;
      if(i == 2);xlabel('Kf (interstitial fluid) [Pa]');end;
      if(i == 3);xlabel('Density grains (kg/m^3)');end;
      if(i == 4);xlabel('Density fluid (kg/m^3)');end;
      if(i == 5);xlabel('Strain hardening');end;
%      set(gca,'YLim',[0.0 0.06],'YTick',[0.0 0.02 0.04 0.06])
      set(gca,'YTickLabel',[]);
      box on;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ALPHA PLOT
%%
if(iar == 1)
   figw = 12;
   figh = 8;
   fig2=figure('visible','on');
   set(fig2,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = NBAND;
   ny = 3;
   xim = 0.01;
   yim = 0.01;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NBAND
     h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);
     h2 = subplot('Position',[loc(1,i+NBAND) loc(2,i+NBAND) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);

     subplot(h1);hold on;box on;
     [n,lim]=hist(alpha(:,(i-1)*order+j),60);n = [0, n, 0];lim = [lim(1) lim lim(end)];
     n = n/sum(n);
     [xx,yy]=stairs(n,lim,'k');
     patch(xx,yy,[0.8,0.8,0.8]);
     stairs(n,lim,'k');
     clear n lim;
     if(j==order);xlabel('Probability');end;
     if(i==1);ylabel('AR coefficient');end;
     if(i>1);set(gca,'YTickLabel',[]);end;
     if(j<order);set(gca,'XTickLabel',[]);end;
     set(gca,'XLim',[0 0.12],'XAxisLocation','top');
     set(gca,'YLim',[-.6 1]);

     idxar = ones(size(alpha(:,(i-1)*order+j)));
     idx = find(alpha(:,(i-1)*order+j) < -0.5);
     idxar(idx) = 0;
     subplot(h2);hold on;box on;
     lim = [-1:1:2];
     [n]=hist(idxar,lim);%n = [0, n, 0];lim = [lim(1) lim lim(end)];
     n = n/sum(n);
     [xx,yy]=stairs(lim-.5,n,'k');
     patch(xx,yy,[0.8,0.8,0.8]);
     stairs(lim-.5,n,'k');
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
%% K PLOT
%%
fig4=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
ylabel('No. interfaces in partition');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(k)])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(k,[0:20]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. interfaces in partition');
ylabel('Probability');
set(gca,'XLim',[-0.5 20.5]);
set(gca,'YLim',[ 0.   1.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL PLOT
%%
fig5=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14);

for i=1:max(A(:,end));
   idx = find(A(:,end)==i);
   plot(A(idx(1:thinstep:end),1),'k');
   clear idx;
%   plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
end;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(logL)])
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14);
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
%% COMPUTE PROFILE MARGINALS
%%
if(imarg == 1)
   NZ = 300;
   NZI= 500;
   nsmooth = ceil(NZ/80.);
   NC = 300;
   NR = 200;
   NA = 200;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
   if(ibucking == 3 & ivgsnosh == 0);
   tlim = pmin(4)+cumsum((pmax(4)-pmin(4))/NA*ones(1,NA));end;
   

   dz  = hmax/(NZ-1);
   z   = cumsum(dz*ones(1,NZ))-dz;
   dzi = hmax/(NZI-1);
   zi  = cumsum(dzi*ones(1,NZI))-dzi;

   h = zeros(1,NZI);
   nlo = zeros(NZ,1);
   nhi = zeros(NZ,1);
   c = zeros(NPROF,NZ);
   r = zeros(NPROF,NZ);
   a = zeros(NPROF,NZ);
   if(ibucking == 3 & ivgsnosh == 0);t = zeros(NPROF,NZ);end;
   disp('Sample size: '),disp(size(m))
   for iprof = 1:NPROF

     if(rem(iprof,5000)==0)
        fprintf(1,'%8i',iprof)
     end
     clear idxh idxc idxr idxa prof;
     %% Find index for current model
     if(k(iprof) > 0)
        idxh = (([1:k(iprof)]-1)*NPL)+1;
        idxc = (([1:k(iprof)]-1)*NPL)+2;
        idxr = (([1:k(iprof)]-1)*NPL)+3;
        idxa = (([1:k(iprof)]-1)*NPL)+4;
        if(ibucking == 3 & ivgsnosh == 0);idxt = (([1:k(iprof)]-1)*NPL)+5;end;
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+NPL-1];
        idxr = [idxr idxr(end)+NPL-1];
        idxa = [idxa idxa(end)+NPL-1];
        if(ibucking == 3 & ivgsnosh == 0);idxt = [idxt idxt(end)+NPL-1];end;
     else
        idxh = [];
        idxc = [1];
        idxr = [2];
        idxa = [3];
        if(ibucking == 3 & ivgsnosh == 0);idxt = [4];end;
     end

     %% Compute the profile for current model
     if(k(iprof) > 0)
        prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
        prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
        prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
     else
        prof(1,1) = hmax;
     end
     prof(:,2) = m(iprof,idxc);
     prof(:,3) = m(iprof,idxr);
     prof(:,4) = m(iprof,idxa);
     if(ibucking == 3 & ivgsnosh == 0);prof(:,5) = m(iprof,idxt);end;

     for ilay=2:k(iprof)+1  %% k is # layers of current model
        idxzi = round(prof(ilay-1,1)/dzi);
        h(idxzi) = h(idxzi) + 1;
     end;
     c(iprof,:) = prof(1,2);
     r(iprof,:) = prof(1,3);
     a(iprof,:) = prof(1,4);
     if(ibucking == 3 & ivgsnosh == 0);t(iprof,:) = prof(1,5);end;
     for ilay=2:k(iprof)+1  %% k is # layers of current model
        idxz = round(prof(ilay-1,1)/dz);
%        h(idxz)     = h(idxz) + 1;
        c(iprof,idxz:end) = prof(ilay,2);
        r(iprof,idxz:end) = prof(ilay,3);
        a(iprof,idxz:end) = prof(ilay,4);
        if(ibucking == 3 & ivgsnosh == 0);t(iprof,idxz:end) = prof(ilay,5);end;
     end;
   end;
   fprintf(1,'\n')
   disp('Done with profiles.');

   %
   % Compute histograms for each depth
   %
   disp('Starting histograms...');
   for iz=1:NZ
%      Nc(iz,:) = histc_tstar(c(:,iz),logL,clim,Tstar);
%      Nr(iz,:) = histc_tstar(r(:,iz),logL,rlim,Tstar);
%      Na(iz,:) = histc_tstar(a(:,iz),logL,alim,Tstar);
      Nc(iz,:) = histc(c(:,iz),clim);
      Nr(iz,:) = histc(r(:,iz),rlim);
      Na(iz,:) = histc(a(:,iz),alim);
      if(ibucking == 3 & ivgsnosh == 0);Nt(iz,:) = histc(t(:,iz),tlim);end;
      if(iboudreau == 1);
        nlo(iz) = sqrt(1-0.36*log(max(c(:,iz))^2))-0.06;
        nhi(iz) = sqrt(1-1.6*log(min(c(:,iz))^2))+0.2;
      end;
      [nfc(iz,:)] = hpd(c(:,iz),100,95);
      [nfr(iz,:)] = hpd(r(:,iz),100,95);
      [nfa(iz,:)] = hpd(a(:,iz),100,95);
      if(ibucking == 3 & ivgsnosh == 0);[nft(iz,:)] = hpd(t(:,iz),100,95);end;
      meac(iz) = median(c(:,iz));
      mear(iz) = median(r(:,iz));
      meaa(iz) = median(a(:,iz));
      if(ibucking == 3 & ivgsnosh == 0);meat(iz) = median(t(:,iz));end;
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         Na(iz,:) = Na(iz,:)/max(Na(iz,:));
        if(ibucking == 3 & ivgsnosh == 0);Nt(iz,:) = Nt(iz,:)/max(Nt(iz,:));end;
      end;
   elseif(inorm == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      Na = Na/max(max(Na));
        if(ibucking == 3 & ivgsnosh == 0);Nt = Nt/max(max(Nt));end;
   end;
   disp('Done histograms.');

if(ibucking == 3 & ivgsnosh == 0);
save tmp.mat z nfc nfa nfr nft meac mear meaa meat miscmead;
else;
save tmp.mat z nfc nfa nfr meac mear meaa miscmead;
end;

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
if(ibucking == 3 & ivgsnosh == 0);t_mean = mean(t);end;
c_mead = median(c);
r_mead = median(r);
a_mead = median(a);
if(ibucking == 3 & ivgsnosh == 0);t_mead = median(t);end;

[ntmp,idxcmax] = max(Nc,[],2);
[ntmp,idxrmax] = max(Nr,[],2);
[ntmp,idxamax] = max(Na,[],2);
if(ibucking == 3 & ivgsnosh == 0);[ntmp,idxtmax] = max(Nt,[],2);end;
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   r_max(iz) = rlim(idxrmax(iz));
   a_max(iz) = alim(idxamax(iz));
   if(ibucking == 3 & ivgsnosh == 0);t_max(iz) = alim(idxtmax(iz));end;
end;
NAVE = 8;
for i=1:length(c_max)
   if(i <= NAVE)
       NAVE1 = 1;
   else
       NAVE1 = i-NAVE;
   end
   if(i >= length(c_max)-NAVE)
       NAVE2 = length(c_max);
   else
       NAVE2 = i+NAVE;
   end
   c_max_sm(i) = sqrt(mean(c_max(NAVE1:NAVE2).^2));
   r_max_sm(i) = sqrt(mean(r_max(NAVE1:NAVE2).^2));
   a_max_sm(i) = sqrt(mean(a_max(NAVE1:NAVE2).^2));
   if(ibucking == 3 & ivgsnosh == 0);t_max_sm(i) = sqrt(mean(t_max(NAVE1:NAVE2).^2));end;
end

end; % end imarg

if(kmax > 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT MAP PROFILES for each K
%%
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
loc = [0.08,  0.18, 0.4533, 0.7266;...
      0.14, 0.14, 0.14, 0.14];
spw1 = 0.09;spw2 = 0.2633;sph = 0.82;

fig66 = figure;hold on; box on;
set(fig66, 'renderer', 'painters')
h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
hold on; box off;
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;

subplot(h4);hold on;box on;
if(imarg == 1)
   h = h/sum(h);
   [xx,yy]=stairs(h,zi,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,zi,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax],'XLim',[0 0.1]);
   set(gca,'YDir','reverse');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;

subplot(h1)
set(gca,'FontSize',14);
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Porosity');
col = [{'k'},{'b'},{'r'},{'c'},{'g'},{'m'},{'k'},{'b'},{'r'},{'c'},{'g'},{'k'},{'m'},{'b'},{'r'},{'c'}];
col = char(col);
for ik = 1:length(nmod_k);
   if(nmod_k(ik)>0);
   plprof(mapk(ik).par,hmax,NPL,col(ik,:),1,0);
   end;
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[0.25 0.5 0.75]);
set(gca,'XTick',[0.25 0.5 0.75]);
box on;

subplot(h2)
set(gca,'FontSize',14);
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
if(ibucking == 1);
  xlabel('Tortuosity');
elseif(ibucking > 2);
  xlabel('log(Kg)');
end;
for ik = 1:length(nmod_k);
   if(nmod_k(ik)>0);
   plprof(mapk(ik).par,hmax,NPL,col(ik,:),2,0);
   end;
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1:1:20]);
set(gca,'XTick',[1:1:20]);
box on;

subplot(h3)
set(gca,'FontSize',14);
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
if(ibucking == 1);
  xlabel('log Transition freq.');
elseif(ibucking >= 2);
  xlabel('strain hardening idx');
end;
for ik = 1:length(nmod_k);
   if(nmod_k(ik)>0);
   plprof(mapk(ik).par,hmax,NPL,col(ik,:),3,0);
   end;
end;
set(gca,'YTickLabel',[]);
box on;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT PROFILE MARGINALS
%%
%nx = 4;
%ny = 1;
%xim = 0.01;
%yim = 0.06;
%xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
%[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

if(ibucking==3 & ivgsnosh==0);
  loc = [0.08,  0.18, 0.385, 0.59, 0.795;...
         0.14, 0.14, 0.14, 0.14, 0.14];
  spw1 = 0.09;
  spw2 = 0.197;
  sph = 0.82;
else
  loc = [0.08,  0.18, 0.4533, 0.7266;...
         0.14, 0.14, 0.14, 0.14];
  spw1 = 0.09;
  spw2 = 0.2633;
  sph = 0.82;
end;

fig6 = figure;hold on; box on;
set(fig6, 'renderer', 'painters')
h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
hold on; box off;
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;
if(ibucking == 3 & ivgsnosh==0);
h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]);
hold on; box off;end;

subplot(h4);hold on;box on;
if(imarg == 1)
   h = h/sum(h);
   [xx,yy]=stairs(h,zi,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,zi,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax],'XLim',[0 0.1]);
   set(gca,'YDir','reverse');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;

subplot(h1)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
   if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Porosity');
if(isyn == 1)
   plprof(mtru,hmax,NPL,'-w',1,0);
   plprof(mtru,hmax,NPL,'--k',1,0);
end;
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',1,0);
end;
if(imax == 1)
  for i=1:length(z);
     [cmx,j] = max(Nc(i,:));
     c_max(i) = clim(j);
  end;
  plot(c_max,z,'w');
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[0.25 0.5 0.75]);
set(gca,'XTick',[0.25 0.5 0.75]);
box on;

subplot(h2)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
if(iboudreau == 1);plot(nlo,z,'-w');end;
if(iboudreau == 1);plot(nhi,z,'-w');end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
if(ibucking == 1);
  xlabel('Tortuosity');
elseif(ibucking >= 2);
  xlabel('log(Kg)');
end;
if(isyn == 1)
   plprof(mtru,hmax,NPL,'-w',2,0);
   plprof(mtru,hmax,NPL,'--k',2,0);
end
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',2);
end
if(imax == 1)
  for i=1:length(z);
     [rmx,j] = max(Nr(i,:));
     r_max(i) = rlim(j);
  end;
  plot(r_max,z,'w');
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1:1:20]);
set(gca,'XTick',[1:1:20]);
box on;

subplot(h3)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(alim,z,Na);shading flat;
   if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
if(ibucking == 1);
  xlabel('log Transition freq.');
elseif(ibucking == 2);
  xlabel('strain hardening idx');
elseif(ibucking==3 & ivgsnosh==0);
  xlabel('strain hardening idx');
elseif(ibucking==3 & ivgsnosh==1);
  xlabel('log(tau idx)');
end;
if(isyn == 1)
   plprof(mtru,hmax,NPL,'-w',3,0);
   plprof(mtru,hmax,NPL,'--k',3,0);
end;
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',3);
end;
if(imax == 1)
  for i=1:length(z);
     [amx,j] = max(Na(i,:));
     a_max(i) = alim(j);
  end;
  plot(a_max,z,'w');
  save('max_model.mat','z','c_max','r_max','a_max');
end;
set(gca,'YTickLabel',[]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
%cmap(1,:) = [1 1 1];
colormap(cmap);
%colorbar('peer',h1,'location','WestOutside');
box on;
if(ibucking==3 & ivgsnosh==0);
  subplot(h5)
  set(gca,'FontSize',14);
  if(imarg == 1)
    pcolor(tlim,z,Nt);shading flat;
    if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
    if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
  end;
  set(h5,'layer','top')
  set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax]);
  set(gca,'YDir','reverse');
  xlabel('log(tau idx)');
  set(gca,'YTickLabel',[]);
  cmap = colormap(jet);
  colormap(cmap);
  box on;
  if(isyn == 1)
    plprof(mtru,hmax,NPL,'-w',4,0);
    plprof(mtru,hmax,NPL,'--k',4,0);
  end;
end;

if(icore == 1)
   barstep1 = 5;
   barstep1r = 5;
   barstep2 = 10;
   barstep2r = 1;
   barstep3 = 10;

   subplot(h1);
%   stairs([c1(:,2);c1(end,2)],[c1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(c1(:,2),c1(:,1),'w','Linewidth',1);
      errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
                 10.*ones(size(c1(1:barstep1:end,2))),...
                 zeros(size(c1(1:barstep1:end,2))),'w','w');
%      plot(c2(:,2),c2(:,1),'g','Linewidth',1);
%      errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%                10.*ones(size(c2(1:barstep2:end,2))),...
%                 zeros(size(c2(1:barstep2:end,2))),'g','g');
%      plot(c3(:,2),c3(:,1),'r','Linewidth',1);
%      errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%                 10.*ones(size(c3(1:barstep3:end,2))),...
%                 zeros(size(c3(1:barstep3:end,2))),'r','r');
%   end
   subplot(h2);
%   stairs([r1(:,2);r1(end,2)],[r1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(r1(:,2),r1(:,1),'w','Linewidth',1);
      errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
                 2./100.*r1(1:barstep1r:end,2),...
                 zeros(size(r1(1:barstep1r:end,2))),'w','w');
      plot(r2(:,2),r2(:,1),'w','Linewidth',1);
      errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
                 2./100.*r2(1:barstep2r:end,2),...
                 zeros(size(r2(1:barstep2r:end,2))),'w','w');
%      plot(r3(:,2),r3(:,1),'r','Linewidth',1);
%      errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%                 2./100.*r3(1:barstep3:end,2),...
%                 zeros(size(r3(1:barstep3:end,2))),'r','r');
%   end
%   subplot(h3);
%   stairs([a1(:,2);a1(end,2)],[a1(:,1);hmax],'k','Linewidth',1.5);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  COMPUTE DATA FIT
%%
if(idatfit == 1)

   for i=1:length(bands);
      stat(i).idxnan  = find(rex(:,i) == 0);
      stat(i).idxgood = find(rex(:,i) ~= 0);
      stat(i).rnan = rex(:,i);
      stat(i).rnan(stat(i).idxnan) = NaN;
   end;
   size(dobs)
   size(rex)
   dobs = dobs .* rex;
   for iband=1:length(bands);
     if(NAVEF > 1)
       if(IGA == 1)
         flo = bands(iband) - bands(iband)*frbw;
         fhi = bands(iband) + bands(iband)*frbw;
       elseif(IGA == 2)
         flo = bands(iband) / (2.^(1./15.)); % Use 1/3 octave
         fhi = bands(iband) * (2.^(1./15.)); %
       elseif(IGA == 3)
         flo = bands(iband) - bands(iband)*frbw;
         fhi = bands(iband) + bands(iband)*frbw;
       else
         flo = bands(iband) - FBW;
         fhi = bands(iband) + FBW;
       end
       fstep = (fhi-flo)/(NAVEF-1);
       fr((iband-1)*NAVEF+1:iband*NAVEF) = flo + ([1:NAVEF]-1) .* fstep;
     else
       fr = bands(iband);
     end
   end
   NFREQ = length(fr);

   restmp=dlmread(resfile);
   resartmp=dlmread(resarfile);
   reptmp=dlmread(repensfile);
   NDAVE = length(reptmp)/NBAND;

   ipass = zeros(NBAND,1);
   ifail = zeros(NBAND,1);
   for j=1:NDAVE;
     if(rem(j,100)==0)
        fprintf(1,'%8i',j);
     end
     l = idxran(j);
     model = m(l,:);
%     if(ibucking == 1); 
%       cp   = zeros(k(l)+1,NFREQ);
%       rho  = zeros(k(l)+1,1);
%       alfp = zeros(k(l)+1,NFREQ);
%       for ilay = 1:k(l)+1;
%         ipar = (ilay-1)*NPL+2;
%         if(ilay == k(l)+1);ipar = (ilay-1)*NPL+1;end;
%         [cp(ilay,:),rho(ilay,1),alfp(ilay,:)] = Buckingham2004(model(end-4),...
%         model(end-3),model(end-2),model(end-1),...
%         model(ipar),model(ipar+1),exp(model(ipar+2)),fr);
%       end;
%       rho = rho/1000.;
%     elseif(ibucking == 2);
%       cp   = zeros(k(l)+1,NFREQ);
%       cs   = zeros(k(l)+1,NFREQ);
%       rho  = zeros(k(l)+1,1);
%       alfp = zeros(k(l)+1,NFREQ);
%       alfs = zeros(k(l)+1,NFREQ);
%       for ilay = 1:k(l)+1;
%         ipar = (ilay-1)*NPL+2;
%         if(ilay == k(l)+1);ipar = (ilay-1)*NPL+1;end;
%         [cp(ilay,:),alfp(ilay,:),cs(ilay,:),alfs(ilay,:),rho(ilay,1)] = ...
%         GS_Mod(model(end-4),model(end-3),model(end-2),model(end-1),...
%         model(ipar),exp(model(ipar+1)),model(ipar+2),fr);
%       end;
%       rho = rho/1000.;
%     elseif(ibucking == 3);
%       if(ivgsnosh==0);
%         cp   = zeros(k(l)+1,NFREQ);
%         cs   = zeros(k(l)+1,NFREQ);
%         rho  = zeros(k(l)+1,1);
%         alfp = zeros(k(l)+1,NFREQ);
%         alfs = zeros(k(l)+1,NFREQ);
%         for ilay = 1:k(l)+1;
%           ipar = (ilay-1)*NPL+2;
%           if(ilay == k(l)+1);ipar = (ilay-1)*NPL+1;end;
%           [cp(ilay,:),alfp(ilay,:),cs(ilay,:),alfs(ilay,:),rho(ilay,1)] = ...
%           VGSlambda_Mod(model(end-4),model(end-3),model(end-2),model(end-1),...
%           model(ipar),exp(model(ipar+1)),model(ipar+2),exp(model(ipar+3)),fr);
%         end;
%         rho = rho/1000.;
%       else
%         cp   = zeros(k(l)+1,NFREQ);
%         cs   = zeros(k(l)+1,NFREQ);
%         rho  = zeros(k(l)+1,1);
%         alfp = zeros(k(l)+1,NFREQ);
%         alfs = zeros(k(l)+1,NFREQ);
%         for ilay = 1:k(l)+1;
%           ipar = (ilay-1)*NPL+2;
%           if(ilay == k(l)+1);ipar = (ilay-1)*NPL+1;end;
%           [cp(ilay,:),alfp(ilay,:),cs(ilay,:),alfs(ilay,:),rho(ilay,1)] = ...
%           VGSlambda_Mod(model(end-5),model(end-4),model(end-3),model(end-2),...
%           model(ipar),exp(model(ipar+1)),model(end-1),exp(model(ipar+2)),fr);
%         end;
%         rho = rho/1000.;
%       end;
%     end;

%     modh = zeros(k(l),1);
%     modh(1) = model(1);
%     for i = 2:k(l);
%        modh(i) = model((i-1)*NPL+1)-model((i-2)*NPL+1);
%     end;
     for iband=1:length(bands);
%       clear geo refl r_ave;
%       for ifr = 1:NAVEF;
%         ifr2 = (iband-1)*NAVEF+ifr;
%         geo = zeros(k(l)+2,4);
%         geo(1,:) = [NaN, cw, 0., rw];
%         for i=1:k(l); 
%            geo(i+1,:) = [modh(i),cp(i,ifr2),...
%                          alfp(i,ifr2),rho(i,1)];
%         end;
%         geo(k(l)+2,:) = [NaN,cp(k(l)+1,ifr2),...
%                          alfp(k(l)+1,ifr2),rho(k(l)+1,1)];        
%         fr2(ifr) = fr(ifr2);
%         if(ispher == 0)
%            [refl(:,ifr)]=ref_nlay3(angobs,geo,fr(ifr2));
%         else
%            [refl(:,ifr), Rp] = spherical_refl(geo,z_t,fr(ifr2),angobs);
%%         end;
%       end;
%       refl = abs(refl);
%       refl = refl.';
%       %% Do freq average
%       if(NAVEF > 1)
%         if(IGA == 1)
%           df = mean(diff(fr2));
%           f0 = bands(iband);
%           for iang = 1:NANG
%             r_ave(iang) = sum(refl(iang,:)' .* ...
%                         exp(-(fr2'-f0).^2/(frbw*f0)^2)*df)/...
%                         sum(exp(-(fr2'-f0).^2/(frbw*f0)^2)*df);
%           end;
%         else
%           r_ave = sqrt(sum(refl.*refl,2)/NAVEF);
%         end;
%       else
%         r_ave = abs(refl);
%       end;
%       ref(:,iband,j) = r_ave;
%       if(ibl==1);
%         ref(:,iband,j) = -20.*log10(abs(ref(:,iband,j)));
%       end;

       ref(:,iband,j) = reptmp(NBAND*(j-1)+iband,:);
       res(:,iband,j) = restmp(NBAND*(j-1)+iband,:);
%       res(:,iband,j) = ref(:,iband,j)-dobs(:,iband);
       res(:,iband,j) = res(:,iband,j)-mean(res(:,iband,j));

       %% Do runs test
       [h(j),p(j)]=runstest(res(:,iband,j),0.);
       if(h(j) == 1);ifail(iband) = ifail(iband)+1;end;
       if(h(j) == 0);ipass(iband) = ipass(iband)+1;end;

       stat(iband).res(:,j) = res(:,iband,j);
       [stat(iband).axx(:,j),stat(iband).axxlags(:,j)] = ...
            xcorr(stat(iband).res(:,j),'coeff');
     %%if(iband==12);save tmp.mat refl fr r_ave angobs;end;
     end;
   end;

   ipass'
   ifail'
   ipass'+ifail'
   ipass'./(ipass'+ifail')
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT DATA MISFIT
   %%
   nx = 4;
   ny = ceil(length(bands)/nx);
   if(length(bands)>12);ny = 3;end;
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig7 = figure;hold on; box on;
   if(length(bands)>12);fig77 = figure;hold on; box on;end;
   ii = 1;
   diffang = diff(angobs)/2;
   angplt = [angobs(1)-2*diffang(1),angobs-[diffang,diffang(end)],...
             angobs(end)+2*diffang(end)];
   for i=1:length(bands);
      figure(fig7);
      jj = i;
      if(i>12);
        figure(fig77);
        jj=i-12;
      end;
      subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      refmean = mean(ref,3);
%      set(gca,'XLim',[angobs(1)-2 angobs(end)+2]);
%      set(gca,'XLim',[angplt(1) angplt(end)]);
      set(gca,'XLim',[15 62]);
      set(gca,'FontSize',14);
      set(gca,'LineWidth',1);
      if(idathist == 1);
         NV = 800;
         vmin = 0.0;
         vmax = 1.1;
         vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
         for iang=1:NANG
            Nd(:,iang) = histc(ref(iang,i,:),vlim);
         end;
         %
         % Normalize Histograms
         %
         if(inorm == 1)
            for iz=1:NANG
               Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
            end;
         elseif(inorm == 2)
            Nd = Nd/max(max(Nd));
         end;
         Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
         pcolor(angplt,vlim,Nd);shading flat;
         clear Nd,vlim;
         plot(angobs,dobs(:,i),'o','MarkerEdgeColor','w',...
                   'MarkerFaceColor','k',...
                   'MarkerSize',4,'Linewidth',1);
         plot(angobs,dobs(:,i),'-w','Linewidth',1.5);
         plot(angobs,dobs(:,i),'--k','Linewidth',1.5);
      else;
         for j=1:length(angobs);
            y = ref(j,i,:);
            [nf(j,:)] = hpd(y,100,95);
         end;
         for j=1:NDAVE;
            plot(angobs,ref(:,i,j),':r');
         end;
%         plot(angobs,refmean(:,i),'-r');
         plot(angobs,nf(:,1),'-w');
         plot(angobs,nf(:,1),'--k');
         plot(angobs,nf(:,2),'-w');
         plot(angobs,nf(:,2),'--k');
%         plot(angobs,dobs(:,i),'-k','Linewidth',1.5);
         plot(angobs,dobs(:,i),'xk','Linewidth',1.5);
      end;
      if(ibl == 0)
         if(max(max(dobs))> 1.)
            ylim  = [0,1.3];
            ytick = [0,0.3,0.6,0.9,1.2];
            text(angobs(end)-10,1.,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         elseif(max(max(dobs))> 0.64)
            ylim  = [0,1.01];
            ytick = [0,0.3,0.6,0.9];
            text(angobs(end)-10,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         else
            ylim  = [0,0.66];
            ytick = [0,0.2,0.4,0.6];
            text(angobs(end)-10,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         end
         set(gca,'YLim',ylim);
         if ((i == (ii-1)*nx+1))
             set(gca,'YTickLabel',ytick,'YTick',ytick);
             ylabel('Refl. Coeff.');
             ii = ii + 1;
         else
              set(gca,'YTickLabel',[],'YTick',ytick);
         end
      else
         text(48,15,[num2str(bands(i)) ' Hz'],'FontSize',12)
         if(max(max(dobs))< 20)
            ylim  = [0,20];
            ytick = [0,5,10,15];
         else
            ylim  = [10,40];
            ytick = [10,20,30];
         end
         set(gca,'YLim',ylim);
         if ((i == 1) | (i == nx+1) | (i == (2*nx)+1))
             set(gca,'YTickLabel',ytick,'YTick',ytick);
             ylabel('Bottom Loss (dB)');
         else
             set(gca,'YTickLabel',[],'YTick',ytick);
         end
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[10:10:90],'XTick',[10:10:90]);
          xlabel('Angle (deg.)');
      else
          set(gca,'XTickLabel',[],'XTick',[10:10:90]);
      end
   end;
   clear nf;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT Axx
   %%
   nx = 4;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig8=figure();
   ii = 1;
   for i=1:length(bands); 
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      clear axxmean;
      axxmean = mean(stat(i).axx,2);
      plot([stat(i).axxlags(1,1) stat(i).axxlags(end,1)],[0 0],'--k');
      plot(stat(i).axxlags(:,1),axxmean,'.-k');
      for j=1:length(axxmean);
         y = stat(i).axx(j,:);
         [nf(j,:)] = hpd(y,100,98);
         clear y;
      end;
      plot(stat(i).axxlags(:,1),nf(:,1),'--k');
      plot(stat(i).axxlags(:,1),nf(:,2),'--k');
      text(stat(i).axxlags(end,1)-20,0.8,[num2str(bands(i)) ' Hz'],'FontSize',12)
      set(gca,'XLim',[stat(i).axxlags(1,1) stat(i).axxlags(end,1)],'YLim',[-.45 1.1]);
      if ((i == (ii-1)*nx+1))
          set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
          ylabel('Autocorrelation');
          ii = ii+1;
      else
           set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[-20 0 20],'XTick',[-20 0 20]);
          xlabel('Lag');
      else
          set(gca,'XTickLabel',[],'XTick',[-20 0 20]);
      end
      clear nf;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT res hist 
   %%
   fig9=figure();
   x = -16.375:.75:16.375;
   xx = -16.25:.01:16.25;
   ii = 1;
   for i=1:length(bands); 
      resmean = mean(stat(i).res(stat(i).idxgood),2);
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      resmean = resmean - mean(resmean);
      if(isd == 1)
         [n1,xout] = hist(resmean./mean(sd(:,i)),x);
      else
         [n1,xout] = hist(resmean./std(resmean),x);
      end
      narea = sum(n1) * (xout(2)-xout(1));
      n1 = n1/narea;
      stairs(xout,n1,'k');
      nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
      plot(xx,nd,'--k')

      text(1,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12)
      axis([-6 6 0 0.6]);
      set(gca,'XTick',[-4 0 4],'FontSize',14);
      set(gca,'XTickLabel',[],'FontSize',14);
      if ((i == (ii-1)*nx+1))
          set(gca,'YTickLabel',[0 0.2 0.4],'YTick',[0 0.2 0.4]);
          ii = ii+1;
      else
           set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[-4 0 4],'XTick',[-4 0 4]);
          xlabel('Res. (std.dev.)');
      else
          set(gca,'XTickLabel',[],'XTick',[-4 0 4]);
      end
      clear resmean;
   end;
end;

if(isave == 1)
%   saveas(fig1,strcat(plotfile1,plotext1),'fig');
   saveas(fig1,strcat(plotfile1,plotext2),'png');
   if(IEPS == 1);
   print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
   end;
   if(iar == 1)
%     saveas(fig2,strcat(plotfile2,plotext1),'fig');
     saveas(fig2,strcat(plotfile2,plotext2),'png');
   if(IEPS == 1);
     print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
   end;
   end
   if(isd == 1)
%      saveas(fig3,strcat(plotfile3,plotext1),'fig');
      saveas(fig3,strcat(plotfile3,plotext2),'png');
   if(IEPS == 1);
      print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
   end;
   end
   if(ibucking > 0)
%      saveas(fig3,strcat(plotfile3,plotext1),'fig');
      saveas(fig33,strcat(plotfile33,plotext2),'png');
   if(IEPS == 1);
      print(fig33,'-painters','-r250',strcat(plotfile33,plotext3),'-depsc');
   end;
   end
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   saveas(fig4,strcat(plotfile4,plotext2),'png');
   if(IEPS == 1);
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');
   end;
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   saveas(fig5,strcat(plotfile5,plotext2),'png');
   if(IEPS == 1);
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');
   end;
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
   saveas(fig6,strcat(plotfile6,plotext2),'png');
   if(IEPS == 1);
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');
   end;
   if(kmax > 0);saveas(fig66,strcat(plotfile66,plotext2),'png');end;
   if(idatfit == 1);
%      saveas(fig7,strcat(plotfile7,plotext1),'fig');
      saveas(fig7,strcat(plotfile7,plotext2),'png');
   if(IEPS == 1);
      print(fig7,'-painters','-r250',strcat(plotfile7,plotext3),'-depsc');
   end;
     if(length(bands)>12);
%       saveas(fig77,strcat(plotfile77,plotext1),'fig');
       saveas(fig77,strcat(plotfile77,plotext2),'png');
   if(IEPS == 1);
       print(fig77,'-painters','-r250',strcat(plotfile77,plotext3),'-depsc');
   end;
     end;
%      saveas(fig8,strcat(plotfile8,plotext1),'fig');
      saveas(fig8,strcat(plotfile8,plotext2),'png');
   if(IEPS == 1);
      print(fig8,'-painters','-r250',strcat(plotfile8,plotext3),'-depsc');
   end;
%      saveas(fig9,strcat(plotfile9,plotext1),'fig');
      saveas(fig9,strcat(plotfile9,plotext2),'png');
   if(IEPS == 1);
      print(fig9,'-painters','-r250',strcat(plotfile9,plotext3),'-depsc');
   end;
   end;
%   if(imarg == 1)
%      save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');
%   end;
end;

return;
