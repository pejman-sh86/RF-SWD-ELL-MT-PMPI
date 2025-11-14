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
ispher  = 1;
iar     = 0;
isd     = 1;
icore   = 0;

filebase    = strrep(filename,'sample.mat','')
datafile    = strrep(filename,'_sample.mat','.txt')
repfile    = strrep(filename,'_sample.mat','_rep.dat')
mapfile    = strrep(filename,'_sample.mat','_map.dat')
profilefile = strcat(filebase,'profile.mat')
plotext1    = 'fig';
plotext2    = 'png';
plotfile1   = strcat(filebase,'chains_logl.')
plotfile2   = strcat(filebase,'ar_par.')
plotfile3   = strcat(filebase,'sigma.')
plotfile4   = strcat(filebase,'number_interfaces.')
plotfile5   = strcat(filebase,'logL.')
plotfile6   = strcat(filebase,'transdim_marg_prof.')
plotfile66  = strcat(filebase,'transdim_map_prof.')
plotfile7   = strcat(filebase,'datafit.')
plotfile8   = strcat(filebase,'axx.')
plotfile9   = strcat(filebase,'reshist.')
corefile    = 'core.mat';

NMISC = 4;
NAVEF = 5;
IGA = 0
frbw = 1.0/15.;
FBW  = 37.5;
%FBW  = 50;
ibl     = 0;
nffile = 'nf_band3.mat';
hmin = 0.05;
logLmin = 220;
logLmax = 280;
logPmin = 140;
logPmax = 200;

%% AUV single pings
%bands = [925, 1075, 1275]
bands = [1000, 1250, 2400]
%bands = [2400,2850,3200]
%bands = [1000, 2400]
%% Site 02  4 m packet
%bands = [300.,400.,504.,635.,800.,1008.,1270,1600]
%bands = [315.,400.,500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.]

%bands = [400, 600,  800,1000,1200,1400,1600];
%% Site 08
%bands = [630, 900, 1000,1600,2400,3600]
%bands = [2100,2300,2500,2700,2900,3100]
%% Site 07
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.];
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000.];
NBAND = length(bands);
NSD = NBAND;
%% AUV:
pmin = [1.450 1.2 0.0]';
pmax = [1.700 2.1 1.0]';
%% Site07:
%pmin = [1550 1.4 0.0]';
%pmax = [1800 2.3 1.0]';
%   NDAVE   = ceil(NPROF/10);
NDAVE   = 100;
thinstep = 1;

mtru = [0.5, 1480., 1.25, 0.010,...
        1.2, 1530., 1.45, 0.085,...
        2.1, 1580., 1.65, 0.300,...
        2.5, 1660., 1.90, 0.550,...
             1600., 1.75, 0.800];

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
%hmax    = 10.;
if(idatfit == 1)
   tmp = dlmread(datafile);
   z_t    = tmp(1,1);
   cw     = tmp(2,1);
   rw     = tmp(3,1);
   hmax   = tmp(4,1)+.5;
%   hmax    = 10.;
   dobs   = tmp(5:length(bands)+4,:)';
   angobs = tmp(length(bands)+5,:);
   NANG   = length(angobs);
%   rex = tmp(length(bands)+6:end,:)';
end

if(icore == 1)
    load(corefile);
end;

NPROF = length(A);
BURNIN = ceil(NPROF/3);
%BURNIN = 2000;
%BURNIN = 1;
A = A(BURNIN:end,:);
NPROF = length(A);
%NPROF = 1e3;

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
    vol = vol * (hmax-(ik-1.)*hmin)*116.712*(pmax(end)-pmin(end));
  end;
  vol = vol * 116.712 * (pmax(end)-pmin(end));  %% This is half-space volume
  logP(i) = -log(vol);
end;
NFP = (k*4)+3;
m = A(:,5:end-5);

%%
%% Save MAP for max(P(k)) to file for starting new inversion
%%
kmin = min(k);
kmax = max(k);
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
[a,b] = max(A(:,1)+logP);map=A(b,4:end-6);
save(mapfile,'map','-ascii');
clear idx;

if(iar == 1)
   (size(A,2)-6-order*NBAND+1:size(A,2)-6)
   alpha = A(:,end-6-(order*NBAND)+1:end-6);
end

logL = A(:,1);
logPPD = A(:,1)+logP;
if(imap == 1);
   [logmap,jmap] = max(logPPD);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

for i = 1:size(A,1);
   idxh = [0:k(i)-1]*4+1;
%   h(i) = sum(m(i,idxh));
   h(i) = m(i,idxh(end));
   clear idxh;
end
if(isd == 1)
   (size(A,2)-6-(order*NBAND)-NSD+1:size(A,2)-6-(order*NBAND))
   sd(:,1:NSD) = A(:,end-6-(order*NBAND)-NSD-NMISC+1:end-6-(order*NBAND)-NMISC);
   disp('Done getting sigma.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig1=figure;
NTH = max(A(:,end))
NPT = max(A(:,end-1))
nx = ceil(sqrt(NTH*NPT))
ny = nx
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
NTH = max(A(:,end))
NPT = max(A(:,end-1))
nx = ceil(sqrt(NTH*NPT))
ny = nx
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
   set(gca,'YLim',[logPmin logPmax]);
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
      set(gca,'YLim',[0.0 0.06],'YTick',[0.0 0.02 0.04 0.06])
      if(i > 1);set(gca,'YTickLabel',[]);end;

      subplot(h2);set(gca,'Layer','top');hold on;
      [n,lim]=hist(sd(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      [xx,yy]=stairs(n,lim,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(n,lim,'k');
      clear n lim;
      if(i == 1);ylabel('Data error standard deviation');end;
      xlabel('rjMCMC probability');
      set(gca,'YLim',[0.0 0.06],'YTick',[0.0 0.02 0.04 0.06])
      if(i > 1);set(gca,'YTickLabel',[]);end;
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
      for j=1:order
         h1 = subplot('Position',[loc(1,i+(j-1)*NBAND) loc(2,i+(j-1)*NBAND) spw sph]);
         hold on; box off;
         set(gca,'FontSize',12);
%      h2 = subplot('Position',[loc(1,i+order) loc(2,i+order) spw sph]);
%      hold on; box off;
%      set(gca,'FontSize',12);
%      subplot(h1);hold on;box on;
%      plot([1:thinstep:length(alpha(:,i))],alpha(1:thinstep:end,i),'k');
%      if(i==1);ylabel('AR coefficient');end;
%      xlabel('rjMCMC step');
%      if(i>1);set(gca,'YTickLabel',[]);end;
%      set(gca,'XLim',[0 length(alpha)])
%      set(gca,'YLim',[-1 1])
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
         set(gca,'XLim',[0 0.08]);
         set(gca,'YLim',[-.6 1]);
      end;
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
set(gca,'YLim',[min(A(:,1))-abs((max(A(:,1))-min(A(:,1)))/20) max(A(:,1))+abs((max(A(:,1))-min(A(:,1)))/20)])

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
%set(gca,'YLim',[min(A(:,1))-abs((max(A(:,1))-min(A(:,1)))/20) max(A(:,1))+abs((max(A(:,1))-min(A(:,1)))/20)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% COMPUTE PROFILE MARGINALS
%%
NPARL = 4;
if(imarg == 1)
   NZ = 600;
   nsmooth = ceil(NZ/80.);
   NC = 400;
   NR = 400;
   NA = 200;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

   dz = hmax/(NZ-1);
   z = cumsum(dz*ones(1,NZ))-dz;

   h = zeros(1,NZ);
   c = zeros(NPROF,NZ);
   r = zeros(NPROF,NZ);
   a = zeros(NPROF,NZ);
   disp('Sample size: '),disp(size(m))
   for iprof = 1:NPROF

     if(rem(iprof,5000)==0)
        fprintf(1,'%8i',iprof)
     end
     clear idxh idxc idxr idxa prof;
     %% Find index for current model
     if(k(iprof) > 0)
        idxh = (([1:k(iprof)]-1)*NPARL)+1;
        idxc = (([1:k(iprof)]-1)*NPARL)+2;
        idxr = (([1:k(iprof)]-1)*NPARL)+3;
        idxa = (([1:k(iprof)]-1)*NPARL)+4;
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+3];
        idxr = [idxr idxr(end)+3];
        idxa = [idxa idxa(end)+3];
     else
        idxh = [];
        idxc = [1];
        idxr = [2];
        idxa = [3];
     end

     %% Compute the profile for current model
     if(k(iprof) > 0)
%        prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
        prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
        prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
     else
        prof(1,1) = hmax;
     end
     prof(:,2) = m(iprof,idxc);
     prof(:,3) = m(iprof,idxr);
     prof(:,4) = m(iprof,idxa);

     c(iprof,:) = prof(1,2);
     r(iprof,:) = prof(1,3);
     a(iprof,:) = prof(1,4);
     for ilay=2:k(iprof)+1  %% k is # layers of current model
        idxz = round(prof(ilay-1,1)/dz);
        h(idxz)     = h(idxz) + 1;
        c(iprof,idxz:end) = prof(ilay,2);
        r(iprof,idxz:end) = prof(ilay,3);
        a(iprof,idxz:end) = prof(ilay,4);
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
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         Na(iz,:) = Na(iz,:)/max(Na(iz,:));
      end;
   elseif(inorm == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      Na = Na/max(max(Na));
   end;
   disp('Done histograms.');

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
c_mead = median(c);
r_mead = median(r);
a_mead = median(a);

[ntmp,idxcmax] = max(Nc,[],2);
[ntmp,idxrmax] = max(Nr,[],2);
[ntmp,idxamax] = max(Na,[],2);
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   r_max(iz) = rlim(idxrmax(iz));
   a_max(iz) = alim(idxamax(iz));
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
end

end; % end imarg
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
   [xx,yy]=stairs(h,z,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,z,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax]);
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
xlabel('Velocity (m/s)');
col = [{'k'},{'b'},{'r'},{'c'},{'g'},{'k'},{'b'},{'r'},{'c'},{'g'},{'k'},{'b'}];
col = char(col);
%for ik = 1:length(nmod_k);
%   plprof(mapk(ik).par,hmax,col(ik,:),1,0);
%end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
box on;

subplot(h2)
set(gca,'FontSize',14);
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
%for ik = 1:length(nmod_k);
%   plprof(mapk(ik).par,hmax,col(ik,:),2,0);
%end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
box on;

subplot(h3)
set(gca,'FontSize',14);
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Attenuation');
%for ik = 1:length(nmod_k);
%   plprof(mapk(ik).par,hmax,col(ik,:),3,0);
%end;
set(gca,'YTickLabel',[]);
box on;

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

loc = [0.08,  0.18, 0.4533, 0.7266;...
      0.14, 0.14, 0.14, 0.14];
spw1 = 0.09;
spw2 = 0.2633;
sph = 0.82;

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

subplot(h4);hold on;box on;
if(imarg == 1)
   h = h/sum(h);
   [xx,yy]=stairs(h,z,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,z,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax]);
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
xlabel('Velocity (m/s)');
if(isyn == 1)
   plprof(mtru,hmax,'-w',1,0);
   plprof(mtru,hmax,'--k',1,0);
end;
if(imap == 1)
   plprof(mmap,hmax,'--k',1,0);
end;
if(imax == 1)
  for i=1:length(z);
     [cmx,j] = max(Nc(i,:));
     c_max(i) = clim(j);
  end;
  plot(c_max,z,'w');
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
box on;

subplot(h2)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
   plprof(mtru,hmax,'-k',2,0);
   plprof(mtru,hmax,'--w',2,0);
end
if(imap == 1)
   plprof(mmap,hmax,'--k',2);
end
if(imax == 1)
  for i=1:length(z);
     [rmx,j] = max(Nr(i,:));
     r_max(i) = rlim(j);
  end;
  plot(r_max,z,'w');
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
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
xlabel('Attenuation');
if(isyn == 1)
   plprof(mtru,hmax,'-k',3,0);
   plprof(mtru,hmax,'--w',3,0);
end;
if(imap == 1)
   plprof(mmap,hmax,'--k',3);
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

if(icore == 1)
   barstep1 = 5;
   barstep1r = 5;
   barstep2 = 10;
   barstep2r = 1;
   barstep3 = 10;

   subplot(h1);
%   stairs([c1(:,2);c1(end,2)],[c1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(c1(:,2),c1(:,1),'-k','Linewidth',1);
      plot(c1(:,2),c1(:,1),'--w','Linewidth',1);
%      errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
%                 10.*ones(size(c1(1:barstep1:end,2))),...
%                 zeros(size(c1(1:barstep1:end,2))),'w','w');
      plot(c2(:,2),c2(:,1),'-k','Linewidth',1);
      plot(c2(:,2),c2(:,1),'--w','Linewidth',1);
%      errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%                10.*ones(size(c2(1:barstep2:end,2))),...
%                 zeros(size(c2(1:barstep2:end,2))),'g','g');
      plot(c3(:,2),c3(:,1),'-k','Linewidth',1);
      plot(c3(:,2),c3(:,1),'--w','Linewidth',1);
%      errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%                 10.*ones(size(c3(1:barstep3:end,2))),...
%                 zeros(size(c3(1:barstep3:end,2))),'r','r');
%   end
   subplot(h2);
%   stairs([r1(:,2);r1(end,2)],[r1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(r1(:,2),r1(:,1),'-k','Linewidth',1);
      plot(r1(:,2),r1(:,1),'--w','Linewidth',1);
%      errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
%                 2./100.*r1(1:barstep1r:end,2),...
%                 zeros(size(r1(1:barstep1r:end,2))),'w','w');
      plot(r2(:,2),r2(:,1),'-k','Linewidth',1);
      plot(r2(:,2),r2(:,1),'--w','Linewidth',1);
%      errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%                 2./100.*r2(1:barstep2r:end,2),...
%                 zeros(size(r2(1:barstep2r:end,2))),'w','w');
      plot(r3(:,2),r3(:,1),'-k','Linewidth',1);
      plot(r3(:,2),r3(:,1),'--w','Linewidth',1);
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

%   for i=1:length(bands);
%      stat(i).idxnan  = find(rex(:,i) == 0);
%      stat(i).idxgood = find(rex(:,i) ~= 0);
%      stat(i).rnan = rex(:,i);
%      stat(i).rnan(stat(i).idxnan) = NaN;
%   end;
%   dobs = dobs .* rex;

   for j=1:NDAVE;
      if(rem(j,500)==0)
         fprintf(1,'%8i',j);
      end
      l = idxran(j);
      model = m(l,:);

      modh = zeros(k(l),1);
      modh(1) = model(1);
      for i = 2:k(l);
         modh(i) = model((i-1)*4+1)-model((i-2)*4+1);
      end;

      clear geo vp alf1dB refl r_ave;
      geo = zeros(k(l)+2,4);
      geo(1,:) = [NaN, cw, 0., rw];
      for i=1:k(l); 
         geo(i+1,:) = [modh(i),1000.*model((i-1)*4+2),...
                       model((i-1)*4+4),model((i-1)*4+3)];
      end;
      geo(k(l)+2,:) = [NaN, 1000.*model(NFP(l)-2),model(NFP(l)),...
                       model(NFP(l)-1)];
%      k(l) = 2;
%      geo = [ NaN, 1511., 0., 1.029;...
%              1.76997149, 1561.0642, 0.977785179, 1.71335923;...
%              5.07706022, 1527.2039, 0.758054738, 1.46706705;...
%               NaN,       1644.7873, 0.0419632979, 2.07418183];
%      afdep(l) = 0.584093125;

      for iband=1:length(bands);
         if(NAVEF > 1)
            if(IGA == 1)
               flo = bands(iband) - bands(iband)*frbw;
               fhi = bands(iband) + bands(iband)*frbw;
            elseif(IGA == 2)
               flo = bands(iband) / (2.^(1./6.)); % Use 1/3 octave
               fhi = bands(iband) * (2.^(1./6.)); %
            else
               flo = bands(iband) - FBW;
               fhi = bands(iband) + FBW;
            end
            fstep = (fhi-flo)/(NAVEF-1);
            fr = flo + ([1:NAVEF]-1) .* fstep;
         else
            fr = bands(iband);
         end
         if(ispher == 0)
            [refl]=ref_nlay3(angobs,geo,fr);
            refl = refl';
            refl = abs(refl);
         else
            [refl, Rp] = spherical_refl(geo,z_t,fr,angobs);
            refl = abs(refl);
         end;
         %% Do freq average
         if(NAVEF > 1)
            if(IGA == 1)
               refl = abs(refl);
               df = mean(diff(fr));
               f0 = bands(iband);
               for iang = 1:NANG
                  r_ave(iang) = sum(refl(iang,:)' .* ...
                                exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
                                sum(exp(-(fr'-f0).^2/(frbw*f0)^2)*df);
               end;
            else
               r_ave = sqrt(sum(refl.*refl,2)/length(fr));
            end;
         else
            r_ave = abs(refl);
         end;
         ref(:,iband,j) = r_ave;
         if(ibl==1);
            ref(:,iband,j) = -20.*log10(abs(ref(:,iband,j)));
         end;
         ref(:,iband,j) = ref(:,iband,j);
         res(:,iband,j) = ref(:,iband,j)-dobs(:,iband);
         stat(iband).res(:,j) = res(:,iband,j);
         [stat(iband).axx(:,j),stat(iband).axxlags(:,j)] = ...
              xcorr(stat(iband).res(:,j),'coeff');
      end;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT DATA MISFIT
   %%
   nx = 3;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig7 = figure;hold on; box on;
   ii = 1;
   diffang = diff(angobs)/2;
   angplt = [angobs(1)-2*diffang(1),angobs-[diffang,diffang(end)],...
             angobs(end)+2*diffang(end)];
   for i=1:length(bands); 
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      refmean = mean(ref,3);
%      set(gca,'XLim',[angobs(1)-2 angobs(end)+2]);
      set(gca,'XLim',[angplt(1) angplt(end)]);
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
      else;
         for j=1:length(angobs);
            y = ref(j,i,:);
            [nf(j,:)] = hpd(y,100,95);
         end;
         for j=1:NDAVE;
            plot(angobs,ref(:,i,j),':r');
         end;
         plot(angobs,refmean(:,i),'-r');
         plot(angobs,nf(:,1),'-k');
         plot(angobs,nf(:,1),'--w');
         plot(angobs,nf(:,2),'-k');
         plot(angobs,nf(:,2),'--w');
      end;
      plot(angobs,dobs(:,i),'o','MarkerEdgeColor','w',...
                'MarkerFaceColor','k',...
                'MarkerSize',6,'Linewidth',2);
      plot(angobs,dobs(:,i),'-k','Linewidth',2);
      plot(angobs,dobs(:,i),'--w','Linewidth',2);
      if(ibl == 0)
         if(max(max(dobs))< 0.64)
            ylim  = [0,0.66];
            ytick = [0,0.2,0.4,0.6];
            text(55,0.5,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[1,1,1])
         else
            ylim  = [0,1.01];
            ytick = [0,0.3,0.6,0.9];
            text(55,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[1,1,1])
         end
         set(gca,'YLim',ylim,'XLim',[angplt(1) angplt(end)]);
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
          set(gca,'XTickLabel',[30 40 50 60],'XTick',[30 40 50 60]);
          xlabel('Angle (deg.)');
      else
          set(gca,'XTickLabel',[],'XTick',[30 40 50 60]);
      end
   end;
   clear nf;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT Axx
   %%
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
      resmean = mean(stat(i).res,2);
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      resmean = resmean - mean(resmean);
      if(isd == 1)
         [n1,xout] = hist(resmean./mean(sd(:,i)),x);
      else
         [n1,xout] = hist(resmean,x);
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
   if(iar == 1)
%     saveas(fig2,strcat(plotfile2,plotext1),'fig');
     saveas(fig2,strcat(plotfile2,plotext2),'png');
   end
   if(isd == 1)
%      saveas(fig3,strcat(plotfile3,plotext1),'fig');
      saveas(fig3,strcat(plotfile3,plotext2),'png');
   end
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   saveas(fig4,strcat(plotfile4,plotext2),'png');
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   saveas(fig5,strcat(plotfile5,plotext2),'png');
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
   saveas(fig6,strcat(plotfile6,plotext2),'png');
   saveas(fig66,strcat(plotfile66,plotext2),'png');
   if(idatfit == 1);
%      saveas(fig7,strcat(plotfile7,plotext1),'fig');
      saveas(fig7,strcat(plotfile7,plotext2),'png');
%      saveas(fig8,strcat(plotfile8,plotext1),'fig');
      saveas(fig8,strcat(plotfile8,plotext2),'png');
%      saveas(fig9,strcat(plotfile9,plotext1),'fig');
      saveas(fig9,strcat(plotfile9,plotext2),'png');
   end;
%   if(imarg == 1)
%      save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');
%   end;
end;

return;
