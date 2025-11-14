function [] = plot_rjhist(filename);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 1; %% Normalize profile marginals line by line
Tstar   = 1.;
isyn    = 0;
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
imapk   = 0;
isave   = 1;
idatfit = 0;

filebase    = strrep(filename,'sample.mat','')
datafile    = strrep(filename,'_sample.mat','.txt')
repfile    = strrep(filename,'_sample.mat','_rep.dat')
mapfile    = strrep(filename,'_sample.mat','_map.dat')
%datafile    = 'x_sim_1_3_16_core.txt'
%datafile    = 'x_s02_1_3_16_mind.txt'
%datafile    = 'x_s02_1_3_25_3lay.txt'
%datafile    = 'x_s02_1_3_16_8m.txt'
%datafile    = 'x_s02_1_3_25_8m.txt'
%datafile    = 'x_s07_1_1_100_noalf.txt'
plotext1    = 'fig';
plotext2    = 'png';
profilefile = strcat(filebase,'profile.mat')
plotfile1   = strcat(filebase,'afdep.')
plotfile2   = strcat(filebase,'chains.')
plotfile3   = strcat(filebase,'sigma.')
plotfile4   = strcat(filebase,'number_interfaces.')
plotfile5   = strcat(filebase,'logL.')
plotfile6   = strcat(filebase,'transdim_marg_prof.')
plotfile66  = strcat(filebase,'transdim_map_prof.')
plotfile7   = strcat(filebase,'data.')
plotfile8   = strcat(filebase,'axx.')
plotfile9   = strcat(filebase,'reshist.')
corefile    = 'core.mat';

%site = str2num(filebase(4:5))
site = 2;
%frange = str2num(filebase(11:12))
frange = 16
NAVEF = 8;

%% Site 02: 
if(site == 2)
   IGA = 2;
   frbw = 1./15.;
   hmax    = 4.25;
   icore   = 1;
   ibl     = 0;
   if(hmax <= 5)
      %%4m:
      pmin = [1450 1.2 0.0]';
      pmax = [1650 2.0 1.0]';
   else
      %% Site 02 8m:
      pmin = [1450 1.2 0.0]';
      pmax = [1700 2.2 1.0]';
   end
   if(frange <= 16)
      bands = [300., 400., 504., 635., 800., 1008., 1270., 1600.];
   else
      bands = [300., 400., 504., 635., 800., 1008., 1270., 1600., 2016., 2540.];
   end;
%   NDAVE   = ceil(NPROF/10);
   NDAVE   = 5000;
   thinstep = 100;
   NDISP = 0;
elseif(site == 7)
   IGA = 2;
   frbw = 0.;
   hmax    = 5.25;
   icore   = 1;
   ibl     = 0;
%   bands = [ 100., 125., 160., 200., 250., 315., 400.];
%   bands = [ 500., 630., 800.,1000.,1250., 1600.,2000.];
   bands = [2500.,3150.,4000.,5000., 6300.,8000.,10000.];
%   bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., 800.,...
%           1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000., 6300.,8000.,10000.];
   %% Site 07:
   pmin = [1500 1.4 0.0]';
   pmax = [1800 2.2 2.0]';
%   datafile    = 'x_s07_1_1_100_noalf.txt'
%   NDAVE   = ceil(NPROF/10);
   NDAVE   = 1e3;
   thinstep = 1;
   NDISP = 1;
else
   IGA = 0
   frbw = 0.;
   hmax    = 6.25;
   icore   = 0;
   ibl     = 0;
   nffile = 'nf_band3.mat';
   bands = [2100,2300,2500,2700,2900,3100]
%   bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., 800.,...
%           1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000., 6300.,8000.,10000.];
%   bands = [1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000., 6300.,8000.,10000.];
%   bands = [ 100., 125., 160., 200., 250.];
%   bands = [ 315., 400., 500., 630., 800.];
%   bands = [1000.,1250., 1600.,2000.,2500];
%   bands = [3150.,4000.,5000., 6300.,8000.,10000.];
   %% Sim:
   pmin = [1450 1.2 0.0]';
   pmax = [1650 2.2 1.5]';
%   NDAVE   = ceil(NPROF/10);
   NDAVE   = 5e2;
   thinstep = 10;
   NDISP = 1;
end

load(filename);

if(icore == 1)
    load(corefile);
end;

NPROF = length(A);
BURNIN = ceil(NPROF/2);
%BURNIN = 60000;
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

if(isyn == 1)
   mtru = [0.60, 1480, 1.30, 0.10, ...
           1.00, 1550, 1.60, 0.30, ...
           0.60, 1585, 1.78, 0.30, ...
                 1600, 1.80, 0.40 ];
   mtrutmp = mtru;
end;

k = A(:,4);
NFP = (k*4)+3;
if(NDISP > 0)
   afdep = A(:,end-5-NDISP+1:end-5);
end;
m = A(:,5:end-5-NDISP);
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
   h(i) = sum(m(i,idxh));
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
[a,b] = max(A(idx,1));map=A(idx(b),4:end-6);
save(mapfile,'map','-ascii');
clear idx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Alpha freq dependence PLOT
%%
if(NDISP > 0)
   fig1=figure;hold on; box on;
   set(gca,'FontSize',14);
   [n,lim]=hist(afdep,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
   n = n/sum(n);
   [xx,yy]=stairs(lim,n,'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(lim,n,'k');
   clear n lim;
   xlabel('n');
   ylabel('Probability');
   set(gca,'Layer','top');
   set(gca,'XLim',[0.0 2.0]);
   set(gca,'YLim',[0 0.05],'YTick',[0,0.02,0.04]);
   box on;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig2=figure;
nx = ceil(sqrt(max(A(:,end))));
ny = nx;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

NTH = max(A(:,end));
for i=1:NTH;
   h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   hold on; box off;
   set(gca,'FontSize',14);

   if(find(A(:,end)==i));
   [AX,H1,H2] = plotyy([1:length(find(A(:,end)==i))],...
                A(find(A(:,end)==i),1),[1:length(find(A(:,end)==i))],...
                A(find(A(:,end)==i),4));
   set(AX(1),'YTick',[220 260 300 340]);
   if(rem(i-1,nx)==0);
      set(get(AX(1),'Ylabel'),'String','logL')
      set(AX(1),'YTickLabel',[220 260 330 340]);
   else;
      set(AX(1),'YTickLabel',[]);
   end;
   set(AX(2),'YTick',[2 4 6 8 10 12 14]);
   if(rem(i,nx)==0);
      set(get(AX(2),'Ylabel'),'String','No. interfaces')
      set(AX(2),'YTickLabel',[2 4 6 8 10 12 14]);
   else;
      set(AX(2),'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(AX(2),'XTickLabel',[],'YLim',[2  8]);
   set(gca,'YLim',[min(A(:,1)) max(A(:,1))]);
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
fig3=figure;
nx = 2;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
set(gca,'FontSize',14);
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
set(gca,'FontSize',14);
subplot(h1);hold on;box on;
plot([1:thinstep:length(sd)],sd(1:thinstep:end),'k');
ylabel('Data error standard deviation');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(sd)])
%set(gca,'YLim',[0.0 0.06])

subplot(h2);hold on;box on;
[n,lim]=hist(sd,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('rjMCMC probability');
set(gca,'YTickLabel',[])
%set(gca,'YLim',[0.0 0.06])

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
plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(logL)])
%set(gca,'YLim',[460 560])
%set(gca,'YLim',[600 700])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(logL,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
%set(gca,'YLim',[460 560])
%set(gca,'YLim',[600 700])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% COMPUTE PROFILE MARGINALS
%%
NPARL = 4;
if(imarg == 1)
   NZ = 400;
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
        prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
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
if(imapk == 1)
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
for ik = 3:6;
   plprof(mapk(ik).par,hmax,col(ik,:),1,0);
end;
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
for ik = 3:6;
   plprof(mapk(ik).par,hmax,col(ik,:),2,0);
end;
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
for ik = 3:6;
   plprof(mapk(ik).par,hmax,col(ik,:),3,0);
end;
set(gca,'YTickLabel',[]);
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
      plot(c1(:,2),c1(:,1),'w','Linewidth',1);
      plot(c1(:,2),c1(:,1),'--k','Linewidth',1);
      errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
                 10.*ones(size(c1(1:barstep1:end,2))),...
                 zeros(size(c1(1:barstep1:end,2))),'k','k');
      plot(c2(:,2),c2(:,1),'w','Linewidth',1);
      plot(c2(:,2),c2(:,1),'--b','Linewidth',1);
      errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
                10.*ones(size(c2(1:barstep2:end,2))),...
                 zeros(size(c2(1:barstep2:end,2))),'b','b');
      plot(c3(:,2),c3(:,1),'w','Linewidth',1);
      plot(c3(:,2),c3(:,1),'--r','Linewidth',1);
      errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
                 10.*ones(size(c3(1:barstep3:end,2))),...
                 zeros(size(c3(1:barstep3:end,2))),'r','r');
%   end
   subplot(h2);
%   stairs([r1(:,2);r1(end,2)],[r1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(r1(:,2),r1(:,1),'w','Linewidth',1);
      plot(r1(:,2),r1(:,1),'--k','Linewidth',1);
      errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
                 2./100.*r1(1:barstep1r:end,2),...
                 zeros(size(r1(1:barstep1r:end,2))),'k','k');
      plot(r2(:,2),r2(:,1),'w','Linewidth',1);
      plot(r2(:,2),r2(:,1),'--b','Linewidth',1);
      errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
                 2./100.*r2(1:barstep2r:end,2),...
                 zeros(size(r2(1:barstep2r:end,2))),'b','b');
      plot(r3(:,2),r3(:,1),'w','Linewidth',1);
      plot(r3(:,2),r3(:,1),'--r','Linewidth',1);
      errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
                 2./100.*r3(1:barstep3:end,2),...
                 zeros(size(r3(1:barstep3:end,2))),'r','r');
%   end
%   subplot(h3);
%   stairs([a1(:,2);a1(end,2)],[a1(:,1);hmax],'k','Linewidth',1.5);
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
set(gca,'XTickLabel',[1500 1600 1700 1800]);
set(gca,'XTick',[1500 1600 1700 1800]);
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
set(gca,'XTickLabel',[1.6 1.8 2.0]);
set(gca,'XTick',[1.6 1.8 2.0]);
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
%xlabel('Attenuation coeff. \kappa');
xlabel('Attenuation (dB/m/kHz)');
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
      plot(c1(:,2),c1(:,1),'w','Linewidth',1);
      plot(c1(:,2),c1(:,1),'--k','Linewidth',1);
      errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
                 10.*ones(size(c1(1:barstep1:end,2))),...
                 zeros(size(c1(1:barstep1:end,2))),'k','k');
      plot(c2(:,2),c2(:,1),'w','Linewidth',1);
      plot(c2(:,2),c2(:,1),'--b','Linewidth',1);
      errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
                10.*ones(size(c2(1:barstep2:end,2))),...
                 zeros(size(c2(1:barstep2:end,2))),'b','b');
      plot(c3(:,2),c3(:,1),'w','Linewidth',1);
      plot(c3(:,2),c3(:,1),'--r','Linewidth',1);
      errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
                 10.*ones(size(c3(1:barstep3:end,2))),...
                 zeros(size(c3(1:barstep3:end,2))),'r','r');
%   end
   subplot(h2);
%   stairs([r1(:,2);r1(end,2)],[r1(:,1);hmax],'k','Linewidth',1.5);
%   if(site == 2)
      plot(r1(:,2),r1(:,1),'w','Linewidth',1);
      plot(r1(:,2),r1(:,1),'--k','Linewidth',1);
      errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
                 2./100.*r1(1:barstep1r:end,2),...
                 zeros(size(r1(1:barstep1r:end,2))),'k','k');
      plot(r2(:,2),r2(:,1),'w','Linewidth',1);
      plot(r2(:,2),r2(:,1),'--b','Linewidth',1);
      errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
                 2./100.*r2(1:barstep2r:end,2),...
                 zeros(size(r2(1:barstep2r:end,2))),'b','b');
      plot(r3(:,2),r3(:,1),'w','Linewidth',1);
      plot(r3(:,2),r3(:,1),'--r','Linewidth',1);
      errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
                 2./100.*r3(1:barstep3:end,2),...
                 zeros(size(r3(1:barstep3:end,2))),'r','r');
%   end
%   subplot(h3);
%   stairs([a1(:,2);a1(end,2)],[a1(:,1);hmax],'k','Linewidth',1.5);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  COMPUTE DATA FIT
%%
if(idatfit == 1)

%   rep = load(repfile);
   tmp = load(datafile);
   dobs   = tmp(1:length(bands),:)';
   angobs = tmp(length(bands)+1,:);
   rex = tmp(length(bands)+2:end,:)';
   for i=1:length(bands);
      stat(i).idxnan  = find(rex(:,i) == 0);
      stat(i).idxgood = find(rex(:,i) ~= 0);
      stat(i).rnan = rex(:,i);
      stat(i).rnan(stat(i).idxnan) = NaN;
   end;
   dobs = dobs .* rex;

   if(NDISP ~= 0)
      ithresh = floor(length(bands)/NDISP)+1;
      fref   = bands(floor(ithresh/2.))
      fthresh= bands(end);
   end

   vpt = zeros(NDAVE,NAVEF*length(bands));
   alfrdBt = zeros(NDAVE,NAVEF*length(bands));
   frt = zeros(1,NAVEF*length(bands));
   for j=1:NDAVE;
      if(rem(j,500)==0)
         fprintf(1,'%8i',j);
      end
      l = idxran(j);
      model = m(l,:);
      clear geo vp alf1dB refl r_ave;
      geo = zeros(k(l)+2,4);
      geo(1,:) = [NaN, 1511., 0., 1.029];
      for i=1:k(l); 
         geo(i+1,:) = [model((i-1)*4+1),model((i-1)*4+2),...
                       model((i-1)*4+4),model((i-1)*4+3)];
      end;
      geo(k(l)+2,:) = [NaN, model(NFP(l)-2),model(NFP(l)),...
                       model(NFP(l)-1)];
%      k(l) = 2;
%      geo = [ NaN, 1511., 0., 1.029;...
%              0.633004308, 1593.9545,  0.61620242,  1.97112622;...
%              0.144615054, 1601.31862, 0.618422141, 2.01869688;...
%              NaN, 1589.17373, 0.334276375, 2.00589373];
%      afdep(l) = 0.584093125;

      for iband=1:length(bands);
         if(IGA == 1)
            flo = bands(iband) - bands(iband)*frbw;
            fhi = bands(iband) + (iband)*frbw;
         elseif(IGA == 2)
            flo = bands(iband) / (2.^(1./6.)); % Use 1/3 octave
            fhi = bands(iband) * (2.^(1./6.)); %
         else
            flo = bands(iband) - FBW;
            fhi = bands(iband) + FBW;
         end
         fstep = (fhi-flo)/(NAVEF-1);
         fr = flo + ([1:NAVEF]-1) .* fstep;
         if(NDISP > 0)
            for iidx = 1:k(l)+1;
               cref = geo(iidx+1,2);
               alfrdB = geo(iidx+1,3);
               y = afdep(l);
               [alf1dB(:,iidx),vp(:,iidx)] = ...
                           KK_Waters_dBmkHz(y,fref,alfrdB,cref,fr);
            end
%            [vp(1,1), vp(end,1)]
            %% over total bandwidth:
            vpt(j,(iband-1)*NAVEF+1:iband*NAVEF)     = vp(:,1);
            alfrdBt(j,(iband-1)*NAVEF+1:iband*NAVEF) = alf1dB(:,1);
            alfrt(j,(iband-1)*NAVEF+1:iband*NAVEF) = alf1dB(:,1)*(1/1000)/(20*log10(exp(1)));

            frt((iband-1)*NAVEF+1:iband*NAVEF) = fr;
            geo2 = geo;
            for ifr = 1:NAVEF
               for iidx = 1:k(l)+1
                  geo2(iidx+1,2) = vp(ifr,iidx);
                  geo2(iidx+1,3) = alf1dB(ifr,iidx);
               end
               [refl(:,ifr)] = ref_nlay3(angobs,geo2,fr(ifr));
            end
         else
            [refl]=ref_nlay3(angobs,geo,fr);
         end;
         refl = abs(refl);
         refl = refl';
         %% Do freq average
         if(IGA == 1)
            refl = abs(refl);
            nang = length(angobs);
            df = mean(diff(fr));
            f0 = fr;
            for iang = 1:nang
               r_ave(iang) = sum(refl(:,iang) .* ...
                             exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
                             sum(exp(-(fr'-f0).^2/(frbw*f0)^2)*df);
            end;
         else
            r_ave = sqrt(sum(refl.*refl,2)/length(fr));
         end;
         ref(:,iband,j) = r_ave;
         if(ibl==1);
            ref(:,iband,j) = -20.*log10(abs(ref(:,iband,j)));
         end;
         ref(:,iband,j) = ref(:,iband,j).*rex(:,iband);
         res(:,iband,j) = ref(:,iband,j)-dobs(:,iband);
         stat(iband).res(:,j) = res(stat(iband).idxgood,iband,j);
         [stat(iband).axx(:,j),stat(iband).axxlags(:,j)] = ...
              xcorr(stat(iband).res(:,j),'coeff');
      end;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT DATA MISFIT
   %%
   nx = 4;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig7 = figure;hold on; box on;
   ii = 1;
   for i=1:length(bands); 
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      refmean = mean(ref,3);
      set(gca,'XLim',[angobs(1)-1 angobs(end)+1]);
      set(gca,'FontSize',14);
      set(gca,'LineWidth',1);
%      plot(rep(length(bands)+1,:),rep(i,:),'-r','LineWidth',2);
      plot(angobs,refmean(:,i).*stat(i).rnan,'-r');
%      plot(angobs,ref(:,1,1).*stat(1).rnan,'--b');
      plot(angobs,dobs(:,i).*stat(i).rnan,'xk');
      for j=1:length(angobs);
         y = ref(j,i,:);
         [nf(j,:)] = hpd(y,100,98);
      end;
      plot(angobs,nf(:,1).*stat(i).rnan,'--k');
      plot(angobs,nf(:,2).*stat(i).rnan,'--k');
      if(ibl == 0)
         text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12)
         if(max(max(dobs))< 0.45)
            ylim  = [0,0.45];
            ytick = [0,0.2,0.4];
         else
            ylim  = [0,1.01];
            ytick = [0,0.3,0.6,0.9];
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
          set(gca,'XTickLabel',[30 60],'XTick',[30 60]);
          xlabel('Angle (deg.)');
      else
          set(gca,'XTickLabel',[],'XTick',[30 60]);
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
      [n1,xout] = hist(resmean./mean(sd),x);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT ALPHA AND CP VS FREQ AT FIXED DEPTH
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear nf; 
if(NDISP > 0)

%   load Biot_halfspace_test.mat

   nx = 2;
   ny = 1;
   xim = 0.06;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig10= figure;hold on; box on;

   %% VP
   subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
   set(gca,'XLim',[frt(1)-10 frt(end)+10]);
   set(gca,'FontSize',14);
   set(gca,'LineWidth',1);
   vptmean = mean(vpt,1);
   plot(frt,vptmean,'-k','LineWidth',1);
   for j=1:length(frt);
      y = vpt(:,j);
      [nfv(j,:)] = hpd(y,100,95);
   end;
   plot(frt,nfv(:,1),'--k');
   plot(frt,nfv(:,2),'--k');
   xlabel('Frequency (Hz)');
   ylabel('Sound Velocity (m/s)');
   set(gca,'XScale','log');
%   plot(biot_f,biot_V(1,:),'-r','Linewidth',2);

   %% ALPHA
   subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
   set(gca,'XLim',[frt(1)-10 frt(end)+10]);
   set(gca,'FontSize',14);
   set(gca,'LineWidth',1);
   alfrtmean = mean(alfrt,1);
   plot(frt,alfrtmean,'-k','LineWidth',1);
   for j=1:length(frt);
      y = alfrt(:,j);
      [nfa(j,:)] = hpd(y,100,95);
   end;
   plot(frt,nfa(:,1),'--k');
   plot(frt,nfa(:,2),'--k');
   xlabel('Frequency (Hz)');
   ylabel('Attenuation (nepers/m/Hz)');
   set(gca,'XScale','log');
%   plot(biot_f,biot_A(1,:)./biot_f,'-r','Linewidth',2);

%   save(nffile,'vptmean','alfrtmean','nfv','nfa','frt');
   
end;

if(isave == 1)
   if(NDISP > 0)
%      saveas(fig1,strcat(plotfile1,plotext1),'fig');
      saveas(fig1,strcat(plotfile1,plotext2),'png');
   end;
%  saveas(fig2,strcat(plotfile2,plotext1),'fig');
  saveas(fig2,strcat(plotfile2,plotext2),'png');
%   saveas(fig3,strcat(plotfile3,plotext1),'fig');
   saveas(fig3,strcat(plotfile3,plotext2),'png');
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   saveas(fig4,strcat(plotfile4,plotext2),'png');
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   saveas(fig5,strcat(plotfile5,plotext2),'png');
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
   saveas(fig6,strcat(plotfile6,plotext2),'png');
   if(imapk == 1)
      saveas(fig66,strcat(plotfile66,plotext2),'png');
   end;
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
