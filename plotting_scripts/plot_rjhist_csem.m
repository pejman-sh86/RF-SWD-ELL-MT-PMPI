function [] = plot_rjhist_csem(filename);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 0; %% Normalize profile marginals line by line
isyn    = 1;
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
isave   = 1;
iar     = 1;
isd     = 1;
NPL = 2;

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
plotfile4   = strcat(filebase,'number_interfaces.');
plotfile5   = strcat(filebase,'logL.');
plotfile6   = strcat(filebase,'transdim_marg_prof.');
plotfile66  = strcat(filebase,'transdim_map_prof.');
plotfile7   = strcat(filebase,'datafit.');
plotfile77  = strcat(filebase,'datafit2.');
plotfile8   = strcat(filebase,'axx.');
plotfile9   = strcat(filebase,'reshist.');
corefile    = 'core.mat';

IEPS = 0;
ibl     = 0;
nffile = 'nf_band3.mat';
hmin = 0.05;

NRF = 3; %% No. of receivers
NSD = NRF;

pmin = [log10(.1)]';
pmax = [log10(300)]';

%   NDAVE   = ceil(NPROF/10);
%NDAVE   = 500;
thinstep = 1;

%% Simulation case:
mtru = [ 50.0,  log10(1.),...
        150.0,  log10(2.),...
                log10(4.)];

load(filename);
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NRF);
   armax = 1.*ones(1,NRF);
else
   order = 1;
end;
hmax    = 500.;

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
logLmin = min(A(:,1))-10.;
logLmax = max(A(:,1))+10;

NPROF = length(A);

k = A(:,4);
%% Compute prior volume
logP = zeros(size(k));
for i=1:length(k);
  logP(i) = log(factorial(k(i)-1)*hmax^-(k(i)-1) * prod(pmax-pmin)^-k(i));
end;
logPPD = A(:,1)+logP(:);

NFP = (k*NPL)+NPL-1;
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
  mapk(i).par=A(idx(b),5:4+A(idx(b),4)*NPL+NPL-1);
  clear idx;
end
[a,b] = max(PPD);map=A(b,4:end-6);
save(mapfile,'map','-ascii');

if(iar == 1)
   (size(A,2)-6-order*NRF+1:size(A,2)-6)
   alpha = A(:,end-6-(order*NRF)+1:end-6);
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
    idxh = [0:k(i)-1]*NPL+1;
    h(i) = m(i,idxh(end));
    clear idxh;
  end;
end;
if(isd == 1)
   (size(A,2)-6-(order*NRF)-NSD+1:size(A,2)-6-(order*NRF))
   sd(:,1:NSD) = A(:,end-6-(order*NRF)-NSD+1:end-6-(order*NRF));
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
   set(AX(1),'YTick',[0:40:10000]);
   if(rem(i-1,nx)==0);
      set(get(AX(1),'Ylabel'),'String','logL') 
      set(AX(1),'YTickLabel',[0:40:10000]);
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
   plot([1:length(idx2)],logPPD(idx1(idx2)));
   set(gca,'YTick',[0:10:10000]);
   if(rem(i-1,nx)==0);
      set(gca,'YTickLabel',[0:10:10000]);
   else;
      set(gca,'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   %set(gca,'YLim',[logLmin+max(logP) logLmax+min(logP)]);
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
      set(gca,'YLim',[0.0 1.e-6],'YTick',[0.0 0.02 0.04 0.06])
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
      set(gca,'YLim',[0.0 1.e-6],'YTick',[0.0 0.02 0.04 0.06])
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
   nx = NRF;
   ny = 3;
   xim = 0.01;
   yim = 0.01;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NRF
      for j=1:order
         h1 = subplot('Position',[loc(1,i+(j-1)*NRF) loc(2,i+(j-1)*NRF) spw sph]);
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
if(imarg == 1);
   NZ = 300;
   NZI= 500;
   nsmooth = ceil(NZ/80.);
   NC = 300;
   NR = 200;
   NA = 200;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   if(NPL == 3);
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));end;
   if(NPL == 4);
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));end;
   if(NPL == 5);
   tlim = pmin(4)+cumsum((pmax(4)-pmin(4))/NA*ones(1,NA));end;
   

   dz  = hmax/(NZ-1);
   z   = cumsum(dz*ones(1,NZ))-dz;
   dzi = hmax/(NZI-1);
   zi  = cumsum(dzi*ones(1,NZI))-dzi;

   h = zeros(1,NZI);
   nlo = zeros(NZ,1);
   nhi = zeros(NZ,1);
   c = zeros(NPROF,NZ);
   if(NPL == 3);r = zeros(NPROF,NZ);end;
   if(NPL == 4);a = zeros(NPROF,NZ);end;
   if(NPL == 5);t = zeros(NPROF,NZ);end;
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
        if(NPL == 3);idxr = (([1:k(iprof)]-1)*NPL)+3;end;
        if(NPL == 4);idxa = (([1:k(iprof)]-1)*NPL)+4;end;
        if(NPL == 5);idxt = (([1:k(iprof)]-1)*NPL)+5;end;
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+NPL-1];
        if(NPL == 3);idxr = [idxr idxr(end)+NPL-1];end;
        if(NPL == 4);idxa = [idxa idxa(end)+NPL-1];end;
        if(NPL == 5);idxt = [idxt idxt(end)+NPL-1];end;
     else
        idxh = [];
        idxc = [1];
        if(NPL == 3);idxr = [2];end;
        if(NPL == 4);idxa = [3];end;
        if(NPL == 5);idxt = [4];end;
     end

     %% Compute the profile for current model
     if(k(iprof) > 0)
        %prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
        prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
        prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
     else
        prof(1,1) = hmax;
     end
     prof(:,2) = m(iprof,idxc);
     if(NPL == 3);prof(:,3) = m(iprof,idxr);end;
     if(NPL == 4);prof(:,4) = m(iprof,idxa);end;
     if(NPL == 5);prof(:,5) = m(iprof,idxt);end;

     for ilay=2:k(iprof)+1  %% k is # layers of current model
        idxzi = round(prof(ilay-1,1)/dzi);
        h(idxzi) = h(idxzi) + 1;
     end;
     c(iprof,:) = prof(1,2);
     if(NPL == 3);r(iprof,:) = prof(1,3);end;
     if(NPL == 4);a(iprof,:) = prof(1,4);end;
     if(NPL == 5);t(iprof,:) = prof(1,5);end;
     for ilay=2:k(iprof)+1  %% k is # layers of current model
        idxz = round(prof(ilay-1,1)/dz);
        if(idxz==0);idxz=1;end;
        c(iprof,idxz:end) = prof(ilay,2);
        if(NPL == 3);r(iprof,idxz:end) = prof(ilay,3);end;
        if(NPL == 4);a(iprof,idxz:end) = prof(ilay,4);end;
        if(NPL == 5);t(iprof,idxz:end) = prof(ilay,5);end;
     end;
   end;
   fprintf(1,'\n')
   disp('Done with profiles.');

   %
   % Compute histograms for each depth
   %
   disp('Starting histograms...');
   for iz=1:NZ
      [Nc(iz,:),binsc] = hist(c(:,iz),clim);
      if(NPL == 3);[Nr(iz,:),bincr] = hist(r(:,iz),rlim);end;
      if(NPL == 4);[Na(iz,:),binca] = hist(a(:,iz),alim);end;
      if(NPL == 5);[Nt(iz,:),binct] = hist(t(:,iz),tlim);end;
      [nfc(iz,:)] = hpd(c(:,iz),100,95);
      if(NPL == 3);[nfr(iz,:)] = hpd(r(:,iz),100,95);end;
      if(NPL == 4);[nfa(iz,:)] = hpd(a(:,iz),100,95);end;
      if(NPL == 5);[nft(iz,:)] = hpd(t(:,iz),100,95);end;
      meac(iz) = median(c(:,iz));
      if(NPL == 3);mear(iz) = median(r(:,iz));end;
      if(NPL == 4);meaa(iz) = median(a(:,iz));end;
      if(NPL == 5);meat(iz) = median(t(:,iz));end;
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 0)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         if(NPL == 3);Nr(iz,:) = Nr(iz,:)/trapz(bincr,Nr(iz,:));end;
         if(NPL == 4);Na(iz,:) = Na(iz,:)/trapz(binca,Na(iz,:));end;
         if(NPL == 5);Nt(iz,:) = Nt(iz,:)/trapz(binct,Nt(iz,:));end;
      end;
   end;
   disp('Done histograms.');

c_mean = mean(c);
if(NPL == 3);r_mean = mean(r);end;
if(NPL == 4);a_mean = mean(a);end;
if(NPL == 5);t_mean = mean(t);end;
c_mead = median(c);
if(NPL == 3);r_mead = median(r);end;
if(NPL == 4);a_mead = median(a);end;
if(NPL == 5);t_mead = median(t);end;

[ntmp,idxcmax] = max(Nc,[],2);
if(NPL == 3);[ntmp,idxrmax] = max(Nr,[],2);end;
if(NPL == 4);[ntmp,idxamax] = max(Na,[],2);end;
if(NPL == 5);[ntmp,idxtmax] = max(Nt,[],2);end;
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   if(NPL == 3);r_max(iz) = rlim(idxrmax(iz));end;
   if(NPL == 4);a_max(iz) = alim(idxamax(iz));end;
   if(NPL == 5);t_max(iz) = alim(idxtmax(iz));end;
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
   if(NPL == 3);r_max_sm(i) = sqrt(mean(r_max(NAVE1:NAVE2).^2));end;
   if(NPL == 4);a_max_sm(i) = sqrt(mean(a_max(NAVE1:NAVE2).^2));end;
   if(NPL == 5);t_max_sm(i) = sqrt(mean(t_max(NAVE1:NAVE2).^2));end;
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
xlabel('Vp');
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

if(NPL == 3);
subplot(h2)
set(gca,'FontSize',14);
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Vp/Vs');
for ik = 1:length(nmod_k);
   if(nmod_k(ik)>0);
   plprof(mapk(ik).par,hmax,NPL,col(ik,:),2,0);
   end;
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1:1:20]);
set(gca,'XTick',[1:1:20]);
box on;
end;

if(NPL == 4);
subplot(h3)
set(gca,'FontSize',14);
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('VP/VS ratio');
for ik = 1:length(nmod_k);
   if(nmod_k(ik)>0);
   plprof(mapk(ik).par,hmax,NPL,col(ik,:),3,0);
   end;
end;
set(gca,'YTickLabel',[]);
box on;
end;
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

if(NPL == 5);
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
if(NPL == 3);
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;end
if(NPL == 4);
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;end;
if(NPL == 5);
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
xlabel('Vp');
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
set(gca,'XTickLabel',[pmin(1):(pmax(1)-pmin(1))/5:pmax(1)]);
set(gca,'XTick',[pmin(1):(pmax(1)-pmin(1))/5:pmax(1)]);
box on;

if(NPL == 3);
subplot(h2)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
xlabel('Vp/Vs');
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
set(gca,'XTickLabel',[pmin(2):(pmax(2)-pmin(2))/5:pmax(2)]);
set(gca,'XTick',[pmin(2):(pmax(2)-pmin(2))/5:pmax(2)]);
box on;
end;

if(NPL == 4);
  subplot(h3)
  set(gca,'FontSize',14);
  if(imarg == 1)
    pcolor(alim,z,Na);shading flat;
    if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
    if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
  end;
  set(h3,'layer','top')
  set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
  set(gca,'YDir','reverse');
  xlabel('VP/VS ratio');
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
  cmap = colormap(jet);
  colormap(cmap);
  box on;
end;
if(NPL == 5);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  SAVE PLOTS
%%
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
end;

return;
