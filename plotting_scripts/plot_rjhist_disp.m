%function [] = plot_rjhist_disp(filename);

filename = 'delta_sample.mat';
imarg   = 1; %% Plot depth-marginal distributions?
inorm1  = 1; %% Normalize profile marginals line by line
inorm2  = 1; %% Normalize profile marginals line by line
Tstar   = 1.;
i_vref  = 1;
isyn    = 1;
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
iratio  = 1;
isave   = 0;   % Save plots?
idatfit = 0;   % Plot data misfit maps?
iar     = 0;
iari    = 0;   % Plot ARI residuals?
NPL     = 3;   % Number parameters per layer
NMISC   = 2;   % Number parameters per layer
logLmin = -150;
logLmax = -50;

filebase    = strrep(filename,'sample.mat','')
%datafile    = 'fort.3'
%datafile    = 'seacable_p001.txt'
%datafile    = 'simulation_nocor.txt'
%datafile    = 'simulation_cor.txt'
%datafile    = 'simulation_arprocsim0.65.txt'
datafile     = strcat(filebase,'dobs.dat');
repfile     = strcat(filebase,'replica.dat');
resfile     = strcat(filebase,'residuals.dat');
resfilear     = strcat(filebase,'residualsar.dat');
resfileari    = strcat(filebase,'residualsari.dat');
vreffile    = strcat(filebase,'vel_ref.txt');
plotfile1   = strcat(filebase,'arma.');
plotfile2   = strcat(filebase,'chain_logLk.');
plotfile3   = strcat(filebase,'sigma.');
plotfile4   = strcat(filebase,'number_interfaces.');
plotfile5   = strcat(filebase,'logL.');
plotfile6   = strcat(filebase,'transdim_marg_prof.');
plotfile7   = strcat(filebase,'data.');
plotfile8   = strcat(filebase,'axx1.');
plotfile88  = strcat(filebase,'axx2.');
plotfile9   = strcat(filebase,'resmap.');
plotfile11  = strcat(filebase,'reshist.');
plotfile12  = strcat(filebase,'armap.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

corefile    = 'core.mat';

%% DELTA SITE: 
%hmax    = 325.;
%% VICTORIA SITE: 
%hmax    = 50.;
%% SIM3 SITE: 
hmax    = 100.;
%% SEACABLE SITE: 
%hmax    = 0.4;

hmax2    = hmax;
%hmax2    = 100.;
%hmax2    = 40.;
%hmax2    = 50.;

icore   = 0;
%pmin = [50     100. 0.5]';
%pmax = [2200  3000  3.0]';
pmin = [ 0.    1.7]';
pmax = [ 1200  1.9]';
if(iratio == 1)
  pmin = [  50   1.42]';
  pmax = [1200  10.00]';
%  pmin = [  50   1.42 1.2]';
%  pmax = [ 300   4.00 2.7]';
%  pmin = [0.01   1.42 1.2]';
%  pmax = [2.50   3.00 3.0]';
end;
thinstep = 100;
%% Delta
%NDATA = 30;
%NDATB = 14;
%% Victoria
NDATA = 50;
NDATB = 0;
%% SIM 3
%NDATA = 35;
%NDATB = 0;
nproc = 2;   % Number of AR processes
MMX = 1;     % Number of AR processes

if(idatfit == 1)
   tmp = dlmread(datafile);
   MMX = tmp(1,1);
   dobs   = tmp(MMX+1:end,2);
   fobs   = tmp(MMX+1:end,1);
   clear tmp;
end;
load(filename);
if(iar == 1)
if(isyn == 1)
   if(size(A,2)==55);order = 1;end;
   if(size(A,2)==56);order = 2;end;
   if(size(A,2)==57);order = 3;end;
   if(size(A,2)==58);order = 4;end;
   order = 1;
   armin = [-.6 -.6 -.6 -.6 -.6 -.6]';
   armax = [ 1.  1.  1.  1.  1.  1.]';
else;
   if(size(A,2)==57);order = 1;end;
   if(size(A,2)==59);order = 2;end;
   if(size(A,2)==61);order = 3;end;
   if(size(A,2)==63);order = 4;end;
   order = 1;
   armin = [-.6 -.6 -.6 -.6 -.6 -.6]';
   armax = [ 1.  1.  1.  1.  1.  1.]';
end;
else
   order = 0;
end;
order
narp = order * nproc

if(icore == 1)
    load(corefile);
end;
if(i_vref == 1);
  vel_ref = dlmread(vreffile);
  vel_ref(1,:) = [];
end;
if(isyn == 1)
  if(i_vref == 0);
%  mtru = [20.0,  100.0, 1200.0, 1.6,...
%          80.0,  400.0, 1500.0, 1.8,...
%         150.0,  800.0, 2000.0, 2.2,...
%                1250.0, 2500.0, 2.5]
%  mtru = [ 6., 100., 2200., 2.0, ...
%          20., 400., 2200., 2.0, ...
%          50., 800., 2200., 2.0, ...
%              1100., 2200., 2.0];
%  mtru = [ 30., 400., 2200., 2.0, ...
%              1100., 2200., 2.0];
%  mtru = [5., 100., 2200., 2.0, ...
%          5., 150., 2200., 2.0, ...
%          5., 200., 2200., 2.0, ...
%          5., 250., 2200., 2.0, ...
%         30., 400., 2200., 2.0, ...
%              800., 2200., 2.0];
  mtru = [10., 150., 2200., 2.0,...
          25., 300., 2200., 2.0,...
               800., 2200., 2.0];
  else
%    ktru = 2;
%    mtru = [20.,   0., 1.80,...
%            40.,-100., 1.80,...
%                   0., 1.80];
    ktru = 3;
    mtru = [ 10., -40., 1.80,...
             30.,  40., 1.80,...
             60.,  40., 1.80,...
                  136., 1.80];
     for ik=1:ktru+1;
      if(ik <= ktru);
        ipar = (ik-1)*NPL+1;
        z = mtru(ipar);
      else;
        ipar = (ik-1)*NPL;
        z = mtru(ipar-NPL+1);
      end;
      [vref]=disp_getref(z,vel_ref);
      mtru(ipar+1) = vref + mtru(ipar+1);
    end;
  end;
end;

NPROF = length(A);
%BURNIN = ceil(NPROF/3);
BURNIN = 1000;
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
NFP = (k*NPL)+NPL-1;
m = A(:,5:end-6-narp);
%% Read in ref velocity and convert from perturbation to Vs
if(i_vref == 1);
  for ismp=1:length(m);
    for ik=1:k(ismp)+1;
      if(ik <= k(ismp));
        ipar = (ik-1)*NPL+1;
        z = m(ismp,ipar);
      else;
        ipar = (ik-1)*NPL;
        z = m(ismp,ipar-NPL+1);
      end;
      [vref]=disp_getref(z,vel_ref);
      m(ismp,ipar+1) = vref + m(ismp,ipar+1);
    end;
  end;
end;

%% Convert ratio to Vp
if(iratio == 0)
for i = 1:length(m);
  idx = NPL*[0:A(i,4)]+NPL-1;
  idx(end) = idx(end)-1;
  m(i,idx) = m(i,idx-1).*m(i,idx);
end;
end;

sd = A(:,end-6-MMX-NMISC-narp:end-7-NMISC-narp);
size(A,2)-6-MMX-NMISC-narp
if(iar == 1)
   (size(A,2)-5-MMX-narp:size(A,2)-6-narp)
   (size(A,2)-6-narp+1:size(A,2)-6)
   alpha = A(:,end-6-narp+1:end-6);
end
logL = A(:,1);
if(imap == 1);
   [logLmap,jmap] = max(logL);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

for i = 1:size(A,1);
   idxh = [0:k(i)-1]*NPL+1;
   h(i) = sum(m(i,idxh));
   clear idxh;
end
disp('Done getting h.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig2=figure;
%nx = ceil(sqrt(max(A(:,end))));
nx = 1;
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
   set(AX(1),'YTick',[-130 -110 -90]);
   if(rem(i-1,nx)==0);
      set(get(AX(1),'Ylabel'),'String','logL')
      set(AX(1),'YTickLabel',[-130 -110 -90]);
   else;
      set(AX(1),'YTickLabel',[]);
   end;
   set(AX(2),'YTick',[2 4 6 8 10]);
   if(rem(i,nx)==0);
      set(get(AX(2),'Ylabel'),'String','No. interfaces')
      set(AX(2),'YTickLabel',[2 4 6 8 10]);
   else;
      set(AX(2),'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(AX(2),'XTickLabel',[],'YLim',[1  6]);
   set(gca,'YLim',[min(A(:,1)) max(A(:,1))]);
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
figw = 4;
figh = 6.6;
fig3=figure('visible','on');
set(fig3,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = 4;
ny = 2;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

for i=1:MMX
   subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   hold on; box off;
   set(gca,'FontSize',12);
   [n,lim]=hist(sd(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
   n = n/sum(n);
   [xx,yy]=stairs(n,lim,'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(n,lim,'k');
   clear n lim;
   xlabel('Probability');
   set(gca,'YTick',[0 2 4 6 8 10]);
   if(i == 1)
     ylabel('Standard deviation');
     set(gca,'YTickLabel',[0 2 4 6 8 10])
   else;
     set(gca,'YTickLabel',[])
   end;
   set(gca,'YLim',[0 10],'XLim',[0 0.07]);
   if(isyn==1);plot([0 1],[5. 5.],'-k');plot([0 1],[5. 5.],'--w');end;
   box on;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ALPHA PLOT
%%
if(iar == 1)
   figw = 4;
   figh = 6.6;
   fig1=figure('visible','on');
   set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = 4;
   ny = 2;
   xim = 0.05;
   yim = 0.08;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:narp
      h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
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
      [n,lim]=hist(alpha(:,i),60);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      [xx,yy]=stairs(n,lim,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(n,lim,'k');
      clear n lim;
      xlabel('Probability');
      if(i==1);ylabel('AR coefficient');end;
      if(i>1);set(gca,'YTickLabel',[]);end;
      set(gca,'YLim',[-0.6 1],'XLim',[0 0.08])
      if(isyn==1);plot([0 1],[0.8 0.8],'-k');plot([0 1],[0.8 0.8],'--w');end;
   end;
   if(order == 1)
      figw = 8;
      figh = 3;
      fig15=figure('visible','on');
      set(fig15,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
      nx = 2;
      ny = 1;
      xim = 0.05;
      yim = 0.01;
      xymarg = [0.1 0.04 0.04 0.1];
      [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
      for i=1:narp
         h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
         hold on; box off;
         set(gca,'FontSize',12);
         NAX = 800;
         axmin = -0.2;
         axmax = 1.01;
         if(i==1)
           lags = [-NDATA:NDATA];
         else
           lags = [-NDATB:NDATB];
         end;
         axlim = axmin+cumsum((axmax-axmin)/NAX*ones(1,NAX));
         C=zeros(length(lags),NPROF);
         for iprof=1:NPROF
            C(:,iprof)=1.0^2.0/(1.0-alpha(iprof,i)^2.0).* ...
            alpha(iprof,i).^(abs(lags));
            C(:,iprof)=C(:,iprof)/max(C(:,iprof));
         end
         for iz=1:size(C,1)
            Naxx(iz,:) = histc(C(iz,:),axlim);
         end;
         %
         % Normalize Histograms
         %
         if(inorm2 == 1)
            for iz=1:size(Naxx,2)
               Naxx(:,iz) = Naxx(:,iz)/max(Naxx(:,iz));
            end;
         elseif(inorm2 == 2)
            Naxx = Naxx/max(max(Naxx));
         end;
         Naxx=flipud(Naxx');
         axlim=flipud(axlim');
         imagesc(lags,axlim,Naxx);
         set(gca,'YDir','normal')
         set(gca,'XLim',[lags(1) lags(end)],'YLim',[axmin axmax]);
         set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
         if(i == 1);
            ylabel('Autocorrelation');
            set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
         end;
         set(gca,'XTickLabel',[-15 0 15],'XTick',[-15 0 15]);
         xlabel('Lag');
         set(gca,'TickDir','out');
         box on;
         if(isyn==1);
           load cov_sim.mat
           plot([-(length(tmp)-1)/2:(length(tmp)-1)/2],tmp/max(tmp),'k');
           plot([-(length(tmp)-1)/2:(length(tmp)-1)/2],tmp/max(tmp),'--w');
         end;
%         plot([-34:34],xcov(dobs,'coeff'),'-w');
         clear Naxx C axlim lags;
      end;
   end;
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
set(gca,'YLim',[0 15])

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',12);
[n,lim]=hist(k,[0:15]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim-0.5;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. interfaces');
ylabel('Probability');
set(gca,'XLim',[-0.5 15.5]);
set(gca,'YLim',[0 .6])

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
set(gca,'YLim',[logLmin logLmax])

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',12);
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
   NZ = 200;
   nsmooth = ceil(NZ/80.);
   NC = 200;
   NR = 100;
   NA = 100;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
%   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

   dz = hmax/(NZ-1);
   z = cumsum(dz*ones(1,NZ))-dz;

   h = zeros(1,NZ);
   c = zeros(NPROF,NZ);
   r = zeros(NPROF,NZ);
%   a = zeros(NPROF,NZ);
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
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+NPL-1];
        if(NPL == 3);
          idxr = (([1:k(iprof)]-1)*NPL)+3;
          %idxa = (([1:k(iprof)]-1)*NPL)+4;
          idxr = [idxr idxr(end)+NPL-1];
          %idxa = [idxa idxa(end)+NPL-1];
        end;
     else
        idxh = [];
        idxc = [1];
        idxr = [2];
        %idxa = [3];
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
     if(NPL == 3);
       prof(:,3) = m(iprof,idxr);
     elseif(NPL == 4);
       prof(:,4) = m(iprof,idxa);
     end;

     c(iprof,:) = prof(1,2);
     r(iprof,:) = prof(1,3);
     %a(iprof,:) = prof(1,4);
     for ilay=2:k(iprof)+1  %% k is # layers of current model
        if(round(prof(ilay-1,1)/dz) > 0);
           idxz = round(prof(ilay-1,1)/dz);
        else;
           idxz = 1;
        end;
        h(idxz)     = h(idxz) + 1;
        c(iprof,idxz:end) = prof(ilay,2);
        r(iprof,idxz:end) = prof(ilay,3);
        %a(iprof,idxz:end) = prof(ilay,4);
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
      %Na(iz,:) = histc(a(:,iz),alim);
   end;
   %
   % Normalize Histograms
   %
   if(inorm1 == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         %Na(iz,:) = Na(iz,:)/max(Na(iz,:));
%         Nc(iz,:) = Nc(iz,:)/NPROF;
%         Nr(iz,:) = Nr(iz,:)/NPROF;
%         Na(iz,:) = Na(iz,:)/NPROF;
      end;
   elseif(inorm2 == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      %Na = Na/max(max(Na));
   end;
   disp('Done histograms.');

c_mean = mean(c);
r_mean = mean(r);
%a_mean = mean(a);
c_mead = median(c);
r_mead = median(r);
%a_mead = median(a);

[ntmp,idxcmax] = max(Nc,[],2);
[ntmp,idxrmax] = max(Nr,[],2);
%[ntmp,idxamax] = max(Na,[],2);
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   r_max(iz) = rlim(idxrmax(iz));
%   a_max(iz) = alim(idxamax(iz));
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
%   a_max_sm(i) = sqrt(mean(a_max(NAVE1:NAVE2).^2));
end

end; % end imarg

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

figw = 12;
figh = 8;
fig6 = figure('visible','on');hold on; box on;
set(fig6,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
set(fig6, 'renderer', 'painters')
h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
hold on; box off;
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
%h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
%hold on; box off;

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
xlabel('Shear-wave Velocity (m/s)');
if(isyn == 1)
   plprof(mtru,hmax,NPL,'-k',1,0);
   plprof(mtru,hmax,NPL,'--w',1,0);
end;
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',1);
end;
if(imax == 1)
  for i=1:length(z);
     [cmx,j] = max(Nc(i,:));
     c_max(i) = clim(j);
  end;
  plot(c_max,z,'w');
end;
set(gca,'YTickLabel',[]);
if(pmax(1)>1200.);
   set(gca,'XTick',[0:300:3000]);
elseif(pmax(1)>400.);
   set(gca,'XTick',[0:200:3000]);
else
   set(gca,'XTick',[0:100:3000]);
end;
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
if(iratio == 1);
  xlabel('Vp/Vs ratio');
  set(gca,'XTickLabel',[1 2 3 4]);
  set(gca,'XTick',[1 2 3 4]);
else;
  xlabel('Compressional-wave Velocity (m/s)');
  set(gca,'XTickLabel',[500 1500 2500]);
  set(gca,'XTick',[500 1500 2500]);
end;
if(isyn == 1)
   plprof(mtru,hmax,NPL,'-k',2,0);
   plprof(mtru,hmax,NPL,'--w',2,0);
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
box on;

% subplot(h3)
% set(gca,'FontSize',12);
% if(imarg == 1)
%    pcolor(alim,z,Na);shading flat;
%    if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
%    if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
% end;
% %surf(alim,z,Na);shading flat;
% set(h3,'layer','top','CLim',[0 1])
% set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
% set(gca,'YDir','reverse');
% xlabel('Density (g/ccm)');
% if(isyn == 1)
%    plprof(mtru,hmax,'-k',3,0);
%    plprof(mtru,hmax,'--w',3,0);
% end;
% if(imap == 1)
%    plprof(mmap,hmax,'--k',3);
% end;
% if(imax == 1)
%   for i=1:length(z);
%      [amx,j] = max(Na(i,:));
%      a_max(i) = alim(j);
%   end;
%   plot(a_max,z,'w');
%   save('max_model.mat','z','c_max','r_max','a_max');
% end;
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[1.5 2.0 2.5]);
% set(gca,'XTick',[1.5 2.0 2.5]);
% %cmap = colormap(flipud(gray));
% cmap = colormap(jet);
% %cmap(1,:) = [1 1 1];
% colormap(cmap);
% %colorbar('peer',h1,'location','WestOutside');
% box on;

if(icore == 1)
   subplot(h1);
   bore = bore/1000.;
   scpt_14 = scpt_14/1000.;
   %if(hmax < 100.);
   %   plot(core(:,2),core(:,1),'-k','Linewidth',2);
   %   plot(core(:,2),core(:,1),'--w','Linewidth',2);
   %else
      plot(scpt_14(:,2),scpt_14(:,1),'-k','Linewidth',2);
      plot(scpt_14(:,2),scpt_14(:,1),'--w','Linewidth',2);
%      plot(scpt_15(:,2),scpt_15(:,1),'-k','Linewidth',2);
%      plot(scpt_15(:,2),scpt_15(:,1),'--w','Linewidth',2);
%      plot(scpt_16(:,2),scpt_16(:,1),'-k','Linewidth',2);
%      plot(scpt_16(:,2),scpt_16(:,1),'--w','Linewidth',2);
      plot(bore(:,2),bore(:,1),'-k','Linewidth',2);
      plot(bore(:,2),bore(:,1),'--w','Linewidth',2);
   %end;
end;
if(isave == 1)
   if(iar == 1);
%      saveas(fig1,strcat(plotfile1,plotext1),'fig');
      print(fig1,'-painters','-r250',strcat(plotfile1,plotext2),'-dpng');
      print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
   end;
%   saveas(fig2,strcat(plotfile2,plotext1),'fig');
   print(fig2,'-painters','-r250',strcat(plotfile2,plotext2),'-dpng');
   print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
%   saveas(fig3,strcat(plotfile3,plotext1),'fig');
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext2),'-dpng');
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext2),'-dpng');
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext2),'-dpng');
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext2),'-dpng');
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');
end;

return;
