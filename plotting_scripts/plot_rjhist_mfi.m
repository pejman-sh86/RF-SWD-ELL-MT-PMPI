function [] = plot_rjhist_mfi2(filename);

imarg   = 1; %% Plot depth-marginal distributions?
inorm1  = 1; %% Normalize profile marginals line by line
inorm2  = 1; %% Normalize profile marginals line by line
isyn    = 0;
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
isave   = 1;   % Save plots?
iar     = 1;
icov    = 0;

filebase    = strrep(filename,'sample.mat','');
%datafile    = 'fort.3';
datafile    = strcat(filebase,'observed.dat');
mapfile    = strrep(filename,'_sample.mat','_map.dat');
repfile     = strcat(filebase,'replica.dat');
resfile     = strcat(filebase,'residuals.dat');
resfilear     = strcat(filebase,'residualsar.dat');
resfileari    = strcat(filebase,'residualsari.dat');
plotfile1   = strcat(filebase,'sigma.');
plotfile2   = strcat(filebase,'arma.');
plotfile3   = strcat(filebase,'sw.');
plotfile4   = strcat(filebase,'number_interfaces.');
plotfile5   = strcat(filebase,'logL.');
plotfile6   = strcat(filebase,'transdim_marg_prof.');
plotfile66  = strcat(filebase,'transdim_map_prof.');
plotfile7   = strcat(filebase,'chains.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

corefile    = 'core.mat';
hmin = 0.5;

%% DELTA SITE: 
hmax    = 62.;
hmax2    = hmax;

icore   = 0;
NBANDS = 9;
bands = [400., 450., 500., 550.,...
         600., 650., 700., 750., 800.];
pmin = [1460  1.3 0.]';
pmax = [1620  2.2 1.0]';
minlimsw = [3.75, 10., 130., 1516., 1510., 1505., 1505.];
maxlimsw = [3.95, 14., 136., 1524., 1518., 1513., 1513.];
thinstep = 1;

load(filename);
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NBANDS);
   armax = 1.00*ones(1,NBANDS);
else
   order = 1;
end;

if(icore == 1)
    load(corefile);
end;

NPROF = length(A);
BURNIN = ceil(NPROF/4);
if(BURNIN >20000);BURNIN = 20000;end;
A = A(BURNIN:end,:);
NPROF = length(A);
logLmin = min(A(:,1))-10;
logLmax = max(A(:,1))+10;

k = A(:,4);
%% Compute prior volume
logP = zeros(size(k));
for i=1:length(k);
  logP(i) = log(factorial(k(i))*hmax^-(k(i)) * prod(pmax-pmin)^-k(i));
end;
logPPD = A(:,1)+logP(:);

NFP = (k*4)+3;
if(icov == 1)
   m = A(:,5:end-7-5-1);
else
   m = A(:,5:end-7-5-order*NBANDS-1);
end

if(iar == 1)
   (size(A,2)-6-order*NBANDS-NBANDS+1:size(A,2)-NBANDS-6)
   alpha = A(:,end-6-(order*NBANDS)-NBANDS+1:end-NBANDS-6);
end
sw = A(:,end-6-(order*NBANDS)-NBANDS-7+1:end-6-(order*NBANDS)-NBANDS);

logL = A(:,1);
if(imap == 1);
   [logLmap,jmap] = max(logL);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

for i = 1:size(A,1);
   sd(i) = m(i,NFP(i)+1);
   idxh = [0:k(i)-1]*4+1;
%   h(i) = sum(m(i,idxh));
   h(i) = m(i,idxh(end));
   clear idxh;
end
disp('Done getting sigma.');

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
%[a,b] = max(A(idx,1));map=A(idx(b),4:end-6);
[a,b] = max(logPPD);map=A(b,4:end-6);
save(mapfile,'map','-ascii');
map
clear idx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig7   = figure;
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
fig77   = figure;
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
   plot([1:length(idx2)],logPPD(idx1(idx2)));
   set(gca,'YTick',[0:40:10000]);
   if(rem(i-1,nx)==0);
      set(gca,'YTickLabel',[0:40:10000]);
   else;
      set(gca,'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(gca,'YLim',[min(logPPD) max(logPPD)]);
   end;
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
figw = 8;
figh = 4;
fig1=figure('visible','on');
set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = 2;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
set(gca,'FontSize',12);
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
set(gca,'FontSize',12);
subplot(h1);hold on;box on;
plot([1:thinstep:length(sd)],sd(1:thinstep:end),'k');
ylabel('Data error standard deviation');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(sd)])
set(gca,'YLim',[0 12])

subplot(h2);hold on;box on;
[n,lim]=hist(sd,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('rjMCMC probability');
set(gca,'YTickLabel',[])
set(gca,'YLim',[0 12])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ALPHA PLOT
%%
if(iar == 1)
   figw = 12;
   figh = 8;
   fig2=figure('visible','on');
   set(fig2,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = NBANDS;
   ny = 3;
   xim = 0.01;
   yim = 0.01;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NBANDS
      for j=1:order
         h1 = subplot('Position',[loc(1,i+(j-1)*NBANDS) loc(2,i+(j-1)*NBANDS) spw sph]);
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
         dx = mean(diff(lim));
         n = n/(sum(n)*dx);
         [xx,yy]=stairs(n,lim,'k');
         patch(xx,yy,[0.8,0.8,0.8]);
         stairs(n,lim,'k');
         clear n lim;
%         if(j==order);xlabel('Probability');end;
         if(i==1);ylabel('AR coefficient');end;
         if(i>1);set(gca,'YTickLabel',[]);end;
         set(gca,'XTickLabel',[]);
         if(j<order);set(gca,'XTickLabel',[]);end;
         text(0.3,-0.48,[num2str(bands(i)) ' Hz'],'FontSize',12)
         set(gca,'XLim',[0 6.]);
         set(gca,'YLim',[-.6 1]);
      end;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Source and Water column  parameters PLOT
%%
figw = 8;
figh = 6;
fig3=figure('visible','on');
set(fig3,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = 4;
ny = 2;
xim = 0.03;
yim = 0.08;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for i=1:size(sw,2)
   if(i<=3)
     h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   else
     h1 = subplot('Position',[loc(1,i+1) loc(2,i+1) spw sph]);
   end
   hold on; box off;
   set(gca,'FontSize',12);

   subplot(h1);hold on;box on;
   lim = [minlimsw(i):(maxlimsw(i)-minlimsw(i))/(100-1):maxlimsw(i)];
   [n]=hist(sw(:,i),lim);n = [0, n, 0];lim = [lim(1) lim lim(end)];
   dx = mean(diff(lim));
   n = n/(sum(n)*dx)*(lim(end)-lim(1));
   [xx,yy]=stairs(lim,n,'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(lim,n,'k');
   clear n lim;
   if(i==1 | i == 4);
     ylabel('Probability');
   end;
   set(gca,'YLim',[0  20])
   if(i==1);xlabel('Range (km)');end;
   if(i==2);xlabel('Source depth (m)');end;
   if(i==3);xlabel('Water depth (m)');end;
   if(i>=4);xlabel(strcat('v(',num2str(i-3),') (m/s)'));end;
%   if(i>1);set(gca,'YTickLabel',[]);end;
   set(gca,'XLim',[minlimsw(i)-((maxlimsw(i)-minlimsw(i))/100) maxlimsw(i)+((maxlimsw(i)-minlimsw(i))/100)])
   if(i>=4);set(gca,'XLim',[minlimsw(i)-((maxlimsw(i)-minlimsw(i))/300) maxlimsw(i)+((maxlimsw(i)-minlimsw(i))/300)]);end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
figw = 5;
figh = 5;
fig4=figure('visible','on');
set(fig4,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',12);
plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
ylabel('No. interfaces');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(k)])
set(gca,'YLim',[0 15.5])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',12);
[n,lim]=hist(k,[0:25]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. interfaces');
ylabel('Probability');
set(gca,'XLim',[-0.5 15.5]);
set(gca,'YLim',[0 .7])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL PLOT
%%
figw = 8;
figh = 4;
fig5=figure('visible','on');
set(fig5,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',12);
plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(logL)])
%set(gca,'YLim',[2200  2300])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',12);
[n,lim]=hist(logL,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
%set(gca,'YLim',[2200 2300])

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
   NA = 400;
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
%        prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
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
        if(round(prof(ilay-1,1)/dz) > 0);
           idxz = round(prof(ilay-1,1)/dz);
        else;
           idxz = 1;
        end;
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
      [Nc(iz,:),binsc] = hist(c(:,iz),clim);
      [Nr(iz,:),binsa] = hist(r(:,iz),rlim);
      [Na(iz,:),binsr] = hist(a(:,iz),alim);
   end;
   %
   % Normalize Histograms (depth by depth)
   %
   if(inorm1 == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
         Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
      end;
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
   set(gca,'Fontsize',14,'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;

subplot(h1)
set(gca,'FontSize',14);
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
col = [{'k'},{'b'},{'r'},{'c'},{'g'},{'k'},{'b'},{'r'},{'c'},{'g'},{'k'},{'b'}];
col = char(col);
%for ik = 1:length(nmod_k);
%   if(length(mapk(ik).par)>1);plprof(mapk(ik).par,hmax2,col(ik,:),1,0);end;
%end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1460:40:1800]);
set(gca,'XTick',[1460:40:1800]);
box on;

subplot(h2)
set(gca,'FontSize',14);
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
%for ik = 1:length(nmod_k);
%   if(length(mapk(ik).par)>1);plprof(mapk(ik).par,hmax2,col(ik,:),2,0);end
%end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
box on;

subplot(h3)
set(gca,'FontSize',14);
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Attenuation');
%for ik = 1:length(nmod_k);
%   if(length(mapk(ik).par)>1);plprof(mapk(ik).par,hmax2,col(ik,:),3,0);end;
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

figw = 9;
figh = 6;
fig6 = figure('visible','on');hold on; box on;
set(fig6,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
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
   [xx,yy]=stairs([0,h],[0,z],'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs([0,h],[0,z],'-k');
   set(gca,'Fontsize',12,'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;


subplot(h1)
set(gca,'FontSize',12);
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
   if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top','CLim',[0 1])
set(gca,'Fontsize',12,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Sound velocity (m/s)');
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',1);
end;
if(imap == 1)
   plprof(mmap,hmax,'--k',1);
end;
if(imax == 1)
  for i=1:length(z);
     [cmx,j] = max(Nc(i,:));
     c_max(i) = clim(j);
  end;
  plot(c_max,z,'w');
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1460:40:1800]);
set(gca,'XTick',[1460:40:1800]);
box on;

subplot(h2)
set(gca,'FontSize',12);
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top','CLim',[0 1])
set(gca,'Fontsize',12,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
   plprof(mtru,hmax,'--k',2);
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
set(gca,'XTickLabel',[1.2:.3:2.2]);
set(gca,'XTick',[1.2:.3:2.2]);
box on;

subplot(h3)
set(gca,'FontSize',12);
if(imarg == 1)
   pcolor(alim,z,Na);shading flat;
   if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top','CLim',[0 1])
set(gca,'Fontsize',12,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)');
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',3);
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
set(gca,'XTickLabel',[0.0 0.5 1.0]);
set(gca,'XTick',[0. 0.5 1.0]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
%cmap(1,:) = [1 1 1];
colormap(cmap);
%colorbar('peer',h1,'location','WestOutside');
box on;

if(icore == 1)
   subplot(h1);
   if(hmax < 100.);
      plot(core(:,2),core(:,1),'-k','Linewidth',2);
      plot(core(:,2),core(:,1),'--w','Linewidth',2);
   else
      plot(scpt_14(:,2),scpt_14(:,1),'-k','Linewidth',2);
      plot(scpt_14(:,2),scpt_14(:,1),'--w','Linewidth',2);
%      plot(scpt_15(:,2),scpt_15(:,1),'-k','Linewidth',2);
%      plot(scpt_15(:,2),scpt_15(:,1),'--w','Linewidth',2);
%      plot(scpt_16(:,2),scpt_16(:,1),'-k','Linewidth',2);
%      plot(scpt_16(:,2),scpt_16(:,1),'--w','Linewidth',2);
      plot(bore(:,2),bore(:,1),'-k','Linewidth',2);
      plot(bore(:,2),bore(:,1),'--w','Linewidth',2);
   end;
end;
if(isave == 1)
   print(fig1,'-painters','-r250',strcat(plotfile1,plotext2),'-dpng');
   print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
   if(iar == 1);
      print(fig2,'-painters','-r250',strcat(plotfile2,plotext2),'-dpng');
      print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
   end;
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext2),'-dpng');
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext2),'-dpng');
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext2),'-dpng');
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext2),'-dpng');
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');
   print(fig7,'-painters','-r250',strcat(plotfile7,plotext2),'-dpng');
   print(fig7,'-painters','-r250',strcat(plotfile7,plotext3),'-depsc');
end;

return;
