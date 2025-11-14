function [] = plot_rjhist_auv_track(istart);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 1; %% Normalize profile marginals line by line to unit area
Tstar   = 1.;
isyn    = 1;
imap    = 1;
imead   = 0;
imean   = 0;
imax    = 0;      %%
isave   = 1;      %%
idatfit = 0;      %% plot data misfits?
ispher  = 0;      %% spherical refle.-coeff?
isd     = 1;      %% sampled over sigma?
itrackmat = 0;

files=dir('*2400_particles.mat');
%files=dir('*2400_sample.mat');

NAVEF = 3;
NPL = 4;
ipst = 1;
ipskip = 1;
kmx = 20;

IGA = 0;
frbw = 0.;
icore   = 0;
ibl     = 0;
nffile = 'nf_band3.mat';
%bands = [1000, 1250, 2400];
bands = [1000, 1200, 2000, 2400];
%bands = [1000, 1200, 2000, 2400, 2800, 3200];
NBAND = length(bands);
NSD = NBAND;
%% Sim:
pmin = [1450 1.1 0.0]';
pmax = [1750 2.2 1.0]';
sigmamin = 0.00;
sigmamax = 0.10;
logLmin = 180;
logLmax = 290;
%   NDAVE   = ceil(NPROF/10);
NDAVE   = 100;
thinstep = 1;
FBW = 37.5;
NZ = 250;
NC = 200;
NR = 200;
NA = 200;

if(isyn == 1);
  track_env = dlmread('track_environment.dat');
end;

NPING = length(files);
%NPING = 5;
ipend = ipst+(NPING*ipskip)-ipskip;
htr = zeros(1,NZ);
ctr = zeros(NPING,NZ);
rtr = zeros(NPING,NZ);
atr = zeros(NPING,NZ);

if(isyn == 1);
   mtrutmp = zeros(size(track_env));
end;
if(itrackmat == 0)
ifile2 = istart-1;
for ifile=ipst:ipskip:ipend;

   ifile2 = ifile2 + 1;
   filename = files(ifile2).name;
   disp(filename);
   filebase    = strrep(filename,'particles.mat','');
   datafile    = strrep(filename,'_particles.mat','.txt');
   repfile     = strrep(filename,'_particles.mat','_replica.mat');
%   filebase    = strrep(filename,'_sample.mat','');
%   datafile    = strrep(filename,'_sample.mat','.txt');
%   repfile     = strrep(filename,'_sample.mat','_replica.mat');
   plotext1    = 'fig';
   plotext2    = 'png';
   plotext3    = 'eps';
   plotfile1   = strcat(filebase,'afdep.');
   plotfile2   = strcat(filebase,'layer_thickness_marginal.');
   plotfile3   = strcat(filebase,'sigma.');
   plotfile4   = strcat(filebase,'number_interfaces.');
   plotfile5   = strcat(filebase,'logL.');
   plotfile6   = strcat(filebase,'transdim_marg_prof.');
   plotfile7   = strcat(filebase,'data.');
   plotfile8   = strcat(filebase,'axx.');
   plotfile9   = strcat(filebase,'reshist.');
   corefile    = 'core.mat';
%   clear A m k logL dobs angobs sd;
   clear A m k logL sd;
   load(filename);
   datafile
   tmp = dlmread(datafile);
   z_t    = tmp(1,1);
   cw     = tmp(2,1);
   rw     = tmp(3,1);
   hmax   = tmp(4,1)+tmp(4,1)/10.;
   hmax2  = hmax;
   dobs   = tmp(5:length(bands)+4,:)';
   angobs = tmp(length(bands)+5,:);
%   angobs = [0:0.25:31];
   nang   = length(angobs);
   nband  = length(bands);
%   dobs = ones(nang,nband);

   if(icore == 1)
       load(corefile);
   end;

   NPROF = length(A);
%   BURNIN = ceil(NPROF/4);
%   A = A(BURNIN:end,:);
   NPROF = length(A);

   %%
   %% Random permutation of sample is only applied to data fit computation
   %%
   idxran = randperm(length(A));
   %A = A(idxran,:);

   if(isyn == 1)
     %% P02:
     ktru = track_env(ifile2,1);
     NFPTRU = ktru*4+3;
     sigtru = track_env(ifile2,NFPTRU+2);
     mtru = track_env(ifile2,2:NFPTRU+1);
     mtrutmp(ifile2,1:NFPTRU) = mtru;
   end;

   k = A(:,4);
   NFP = (k*4)+3;
   m = A(:,5:end);
   logL = A(:,1);
   if(imap == 1);
      clear logLmap jmap kmap NFPmap mmap;
      [logLmap,jmap] = max(logL);
      kmap = k(jmap);
      NFPmap = NFP(jmap);
      mmap = m(jmap,:);
   end;

  %% Marginalize layer thickness over whole model
  if(kmx > 0);
    h = zeros(size(A,1),max(A(:,4)));
    dep = zeros(size(A,1),max(A(:,4)));
    for i = 1:size(A,1);
      if(A(:,4)>0);
        idxh = [0:k(i)-1]*NPL+1;
        dep(i,1:k(i)) = m(i,idxh);
        h(i,1:k(i)) = diff([0,dep(i,1:k(i))]);
        clear idxh;
      end;
    end;
  end;
  dep = dep(:);dep(find(dep==0))=[];
  h = h(:);h(find(h==0))=[];
  [a_enos,b_enos]=hist(h,200);a_enos=a_enos/trapz(b_enos,a_enos);
  a_enos=[0,a_enos,0];
  b_enos=[b_enos(1),b_enos,b_enos(end)];
  fig2 = figure('visible','off');
  [xx,yy]=stairs(b_enos,a_enos,'k');
  patch(xx,yy,[0.8,0.8,0.8]);

   sd = A(:,end-NSD+1:end);
%   sd = A(:,88:88+NSD);
   for i = 1:size(A,1);
      idxh = [0:k(i)-1]*4+1;
      h(i) = sum(m(i,idxh));
      clear idxh;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% SIGMA PLOT
   %%
   if(isd == 1)
      fig3=figure('visible','off');
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
         xlabel('Particle no.');
         set(gca,'XLim',[0 length(sd(:,i))])
         set(gca,'YLim',[sigmamin sigmamax],'YTick',[0.02 0.04 0.06])
         if(i > 1);set(gca,'YTickLabel',[]);end;
      
         subplot(h2);set(gca,'Layer','top');hold on;
%        [n,lim]=hist(sd(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
         [n,lim]=hist(sd(:,i),[sigmamin:(sigmamax-sigmamin)/100:sigmamax]);
         n = [0, n, 0];lim = [lim(1) lim lim(end)];n = n/sum(n);
         sigmagl(ifile2,i,:) = n;
         sigmalimgl = lim;
         [xx,yy]=stairs(n,lim,'k');
         patch(xx,yy,[0.8,0.8,0.8]);
         stairs(n,lim,'k');
         clear n lim;
         if(i == 1);ylabel('Data error standard deviation');end;
         xlabel('Probability');
         set(gca,'YLim',[sigmamin sigmamax],'YTick',[0.02 0.04 0.06])
         if(i > 1);set(gca,'YTickLabel',[]);end;
         box on;
      end;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% K PLOT
   %%
   fig4=figure('visible','off');
   subplot(2,1,1);hold on;box on;
   set(gca,'FontSize',14);
   plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
   ylabel('No. interfaces in partition');
   xlabel('Particle no.');
%   set(gca,'XLim',[0 length(k)])
   set(gca,'YLim',[0 10])

   subplot(2,1,2);hold on;box on;
   set(gca,'FontSize',14);
   [n,lim]=hist(k,[0:kmx]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
   n = n/sum(n);
   lim = lim-0.5;
   kgl(ifile2,:) = n;
   limgl = lim;
   [xx,yy]=stairs(lim,n,'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(lim,n,'k');
   clear n lim;
   xlabel('No. interfaces in partition');
   ylabel('Probability');
   set(gca,'XLim',[-0.5 kmx+0.5],'YLim',[0. 1.0]);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% logL PLOT
   %%
   fig5=figure('visible','off');
   subplot(1,2,1);hold on;box on;
   set(gca,'FontSize',14);
   plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
   ylabel('log Likelihood');
   xlabel('Particle no.');
   set(gca,'XLim',[0 length(logL)])
   set(gca,'YLim',[logLmin logLmax])
   %set(gca,'YLim',[150 400])
  
   subplot(1,2,2);hold on;box on;
   set(gca,'FontSize',14);
   [n,lim]=hist(logL,[logLmin:(logLmax-logLmin)/100:logLmax]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
   n = n/sum(n);
   logLgl(ifile2,:) = n;
   logLlimgl = lim;
   [xx,yy]=stairs(n,lim,'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(n,lim,'k');
   clear n lim;
   xlabel('Probability');
   set(gca,'YLim',[logLmin logLmax])
   %set(gca,'YLim',[150 400])

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% COMPUTE PROFILE MARGINALS
   %%
   NPARL = 4;
   if(imarg == 1)
      clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
      rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
      alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

      dz = hmax/(NZ-1);
      zlim = [0:dz:hmax];
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
        if(idxz > NZ);idxz=NZ;end;
        h(idxz)     = h(idxz) + 1;
        c(iprof,idxz:end) = prof(ilay,2);
        r(iprof,idxz:end) = prof(ilay,3);
        a(iprof,idxz:end) = prof(ilay,4);
     end;
   end;
   fprintf(1,'\n')

   %
   % Compute histograms for each depth
   %
   for iz=1:NZ
%      Nc(iz,:) = histc_tstar(c(:,iz),logL,clim,Tstar);
%      Nr(iz,:) = histc_tstar(r(:,iz),logL,rlim,Tstar);
%      Na(iz,:) = histc_tstar(a(:,iz),logL,alim,Tstar);
      [Nc(iz,:),binsc] = hist(c(:,iz),clim);
      [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
      [Na(iz,:),binsa] = hist(a(:,iz),alim);
      %% Compute 95% cred int
      y = c(:,iz);
      [nfc(iz,ifile2,:)] = hpd(y,100,95);
      y = r(:,iz);
      [nfr(iz,ifile2,:)] = hpd(y,100,95);
      y = a(:,iz);
      [nfa(iz,ifile2,:)] = hpd(y,100,95);
   end;
   [Ndep,binsdep] = hist(dep,zlim);
   %
   % Normalize Histograms
   %
   if(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
         Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
      end;
       Ndep = Ndep/trapz(binsdep,Ndep);
       Ndep = [0,Ndep,0];binsdep = [binsdep(1),binsdep,binsdep(end)];
   elseif(inorm == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      Na = Na/max(max(Na));
   end;
   Ncgl(:,:,ifile2) = Nc;
   Nrgl(:,:,ifile2) = Nr;
   Nagl(:,:,ifile2) = Na;

c_mean(:,ifile2) = mean(c);
r_mean(:,ifile2) = mean(r);
a_mean(:,ifile2) = mean(a);
c_mead(:,ifile2) = median(c);
r_mead(:,ifile2) = median(r);
a_mead(:,ifile2) = median(a);

[ntmp,idxcmax] = max(Nc,[],2);
[ntmp,idxrmax] = max(Nr,[],2);
[ntmp,idxamax] = max(Na,[],2);
for iz=1:NZ
   c_max(iz,ifile2) = clim(idxcmax(iz));
   r_max(iz,ifile2) = rlim(idxrmax(iz));
   a_max(iz,ifile2) = alim(idxamax(iz));
end;
NAVE = 8;
for i=1:length(c_max(:,ifile2))
   if(i <= NAVE)
       NAVE1 = 1;
   else
       NAVE1 = i-NAVE;
   end
   if(i >= length(c_max(:,ifile2))-NAVE)
       NAVE2 = length(c_max(:,ifile2));
   else
       NAVE2 = i+NAVE;
   end
   c_max_sm(i,ifile2) = mean(c_max(NAVE1:NAVE2,ifile2));
   r_max_sm(i,ifile2) = mean(r_max(NAVE1:NAVE2,ifile2));
   a_max_sm(i,ifile2) = mean(a_max(NAVE1:NAVE2,ifile2));
   c_mead_sm(i,ifile2) = mean(c_mead(NAVE1:NAVE2,ifile2));
   r_mead_sm(i,ifile2) = mean(r_mead(NAVE1:NAVE2,ifile2));
   a_mead_sm(i,ifile2) = mean(a_mead(NAVE1:NAVE2,ifile2));
   c_mean_sm(i,ifile2) = mean(c_mean(NAVE1:NAVE2,ifile2));
   r_mean_sm(i,ifile2) = mean(r_mean(NAVE1:NAVE2,ifile2));
   a_mean_sm(i,ifile2) = mean(a_mean(NAVE1:NAVE2,ifile2));
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

fig6 = figure('visible','off');hold on; box on;
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
   %hgl(ifile2,:) = h/sum(h);
   %h = h/sum(h);
   hgl(ifile2,:) = Ndep;
   [xx,yy]=stairs(Ndep,binsdep,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   %[xx,yy]=stairs(h,z,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax2]);
   set(gca,'YDir','reverse','XLim',[0 5.0]);
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;


subplot(h1)
set(gca,'FontSize',14);
if(imarg == 1)
   imagesc(clim,z,Nc);shading flat;
   if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
if(isyn == 1)
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',1,0);
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',1,0);
end;
%if(imap == 1)
%   plprof(mmap,hmax,'--k',1,0);
%end;
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
   imagesc(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',2,0);
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',2,0);
end
%if(imap == 1)
%   plprof(mmap,hmax,'--k',2,0);
%end
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
   imagesc(alim,z,Na);shading flat;
   if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Attenuation');
if(isyn == 1)
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'-k',3,0);
   plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'--w',3,0);
end;
%if(imap == 1)
%   plprof(mmap,hmax2,'--k',3,0);
%end;
if(imax == 1)
  for i=1:length(z);
     [amx,j] = max(Na(i,:));
     a_max(i) = alim(j);
  end;
  plot(a_max,z,'w');
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

%   for i=1:length(bands);
%      stat(i).idxnan  = find(rex(:,i) == 0);
%      stat(i).idxgood = find(rex(:,i) ~= 0);
%      stat(i).rnan = rex(:,i);
%      stat(i).rnan(stat(i).idxnan) = NaN;
%   end;

   for j=1:NDAVE;
      if(rem(j,500)==0)
         fprintf(1,'%8i',j);
      end
      l = idxran(j);
      if(k(l) == 0);l=l+1;end;
      model = m(l,:);
      clear geo vp alf1dB refl r_ave;
      geo = zeros(k(l)+2,4);
      geo(1,:) = [NaN, cw, 0., rw];
      for i=1:k(l); 
         geo(i+1,:) = [model((i-1)*4+1),1000.*model((i-1)*4+2),...
                       model((i-1)*4+4),model((i-1)*4+3)];
      end;
      geo(k(l)+2,:) = [NaN, 1000.*model(NFP(l)-2),model(NFP(l)),...
                       model(NFP(l)-1)];
%      k(l) = 2;
%      geo = [ NaN, 1511., 0., 1.029;...
%              0.633004308, 1593.9545,  0.61620242,  1.97112622;...
%              0.144615054, 1601.31862, 0.618422141, 2.01869688;...
%              NaN, 1589.17373, 0.334276375, 2.00589373];
%      afdep(l) = 0.584093125;

      for iband=1:length(bands);
         if(NAVEF > 1)
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
         else
            fr = bands(iband);
         end
         if(ispher == 0)
            [refl]=ref_nlay3(angobs,geo,fr);
            refl = abs(refl);
%            size(geo)
%            k(l)
%            model
%            size(fr)
%            size(angobs)
%            size(refl)
            if(size(refl,1)==3);
              refl = refl';
            end;
         else
            [refl, Rp] = spherical_refl(geo,z_t,fr,angobs);
            refl = abs(refl);
         end;
%         size(refl)
         %% Do freq average
         if(NAVEF > 1)
            if(IGA == 1)
               refl = abs(refl);
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
         else
            r_ave = abs(refl);
         end;
         
%         size(r_ave)
         ref(1:nang,iband,j) = r_ave;
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
   for iband=1:length(bands);
      for j=1:length(angobs);
         y = ref(j,iband,:);
         [nf(j,:,iband,ifile2)] = hpd(y,100,95);
      end;
   end;

   if(imap == 1);
      model = mmap;
      clear geo vp alf1dB refl r_ave;
      geo = zeros(kmap+2,4);
      geo(1,:) = [NaN, cw, 0., rw];
      for i=1:kmap; 
         geo(i+1,:) = [model((i-1)*4+1),model((i-1)*4+2),...
                       model((i-1)*4+4),model((i-1)*4+3)];
      end;
      geo(kmap+2,:) = [NaN, model(NFPmap-2),model(NFPmap),...
                       model(NFPmap-1)];
%      k(l) = 2;
%      geo = [ NaN, 1511., 0., 1.029;...
%              0.633004308, 1593.9545,  0.61620242,  1.97112622;...
%              0.144615054, 1601.31862, 0.618422141, 2.01869688;...
%              NaN, 1589.17373, 0.334276375, 2.00589373];
%      afdep(l) = 0.584093125;

      for iband=1:length(bands);
         if(NAVEF > 1)
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
         else
            fr = bands(iband);
         end
         if(ispher == 0)
            [refl]=ref_nlay3(angobs,geo,fr);
            refl = abs(refl);
            refl = refl';
         else
            [refl, Rp] = spherical_refl(geo,z_t,fr,angobs);
            refl = abs(refl);
         end;
         %% Do freq average
         if(NAVEF > 1)
            if(IGA == 1)
               refl = abs(refl);
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
         else
            r_ave = abs(refl);
         end;
         refmap(:,iband,ifile2) = r_ave;
      end;
   end;
   dobsmap(:,:,ifile2) = dobs;
   angobsmap(:,:,ifile2) = angobs;
   save(repfile,'ref','dobs','angobs');

fig77=figure('visible','on');
for i=1:3;subplot(1,3,i);hold off;
  plot(angobsmap(1,:,ifile2),dobsmap(:,i,ifile2),'xk');hold on;
  plot(angobsmap(1,:,ifile2),refmap(:,i,ifile2));
  plot(angobsmap(1,:,ifile2),nf(:,1,i,ifile2),'--k');
  plot(angobsmap(1,:,ifile2),nf(:,2,i,ifile2),'--k');
  set(gca,'YLim',[0 1],'XLim',[28 70]);
end;%pause(3);
close(fig77);

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
   figw = 8.5;
   figh = 2.5;
   fig7 = figure('visible','off');hold on; box on;
   set(fig7,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   ii = 1;
   for i=1:length(bands); 
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      refmean = mean(ref,3);
      set(gca,'XLim',[angobs(1)-1 angobs(end)+1]);
      set(gca,'FontSize',14);
      set(gca,'LineWidth',1);
%      plot(rep(length(bands)+1,:),rep(i,:),'-r','LineWidth',2);
      for j=1:NDAVE;
         plot(angobs,ref(:,i,j),':r');
      end;
      plot(angobs,refmean(:,i),'-r');
%      plot(angobs,ref(:,1,1).*stat(1).rnan,'--b');
      plot(angobs,dobs(:,i),'xk');
      plot(angobs,nf(:,1,iband,ifile2),'--k');
      plot(angobs,nf(:,2,iband,ifile2),'--k');
      if(ibl == 0)
%         if(max(max(dobs))< 0.64)
            ylim  = [0,0.75];
            ytick = [0,0.2,0.4,0.6];
            text(52,0.65,[num2str(bands(i)) ' Hz'],'FontSize',12)
%         else
%            ylim  = [0,1.01];
%            ytick = [0,0.3,0.6,0.9];
%            text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12)
%         end
         set(gca,'YLim',ylim,'XLim',[27 68]);
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

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT Axx
   %%
   fig8=figure('visible','off');
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
         [nfaxx(j,:)] = hpd(y,100,98);
         clear y;
      end;
      plot(stat(i).axxlags(:,1),nfaxx(:,1),'--k');
      plot(stat(i).axxlags(:,1),nfaxx(:,2),'--k');
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
      clear nfaxx;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT res hist 
   %%
   fig9=figure('visible','off');
   x = -16.375:.75:16.375;
   xx = -16.25:.01:16.25;
   ii = 1;
   for i=1:length(bands); 
      resmean = mean(stat(i).res,2);
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      resmean = resmean - mean(resmean);
      [n1,xout] = hist(resmean./mean(sd(:,i)),x);
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
%   saveas(fig2,strcat(plotfile2,plotext1),'fig');
   saveas(fig2,strcat(plotfile2,plotext2),'png');
%   saveas(fig3,strcat(plotfile3,plotext1),'fig');
   saveas(fig3,strcat(plotfile3,plotext2),'png');
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   saveas(fig4,strcat(plotfile4,plotext2),'png');
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   saveas(fig5,strcat(plotfile5,plotext2),'png');
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
   saveas(fig6,strcat(plotfile6,plotext2),'png');
   if(idatfit == 1);
%    saveas(fig7,strcat(plotfile7,plotext1),'fig');
    saveas(fig7,strcat(plotfile7,plotext2),'png');
%    saveas(fig7,strcat(plotfile7,'eps'),'epsc');
%    saveas(fig8,strcat(plotfile8,plotext1),'fig');
    saveas(fig8,strcat(plotfile8,plotext2),'png');
%    saveas(fig9,strcat(plotfile9,plotext1),'fig');
    saveas(fig9,strcat(plotfile9,plotext2),'png');
    close(fig7);
    close(fig8);
    close(fig9);
  end;
close(fig3);
close(fig4);
close(fig5);
close(fig6);
end;
   %%
   %% Plot track profile
   %%
   clear idxh idxc idxr idxa prof;
   if(isyn == 1)
      %% Find index for current model
      if(ktru > 0)
         idxh = (([1:ktru]-1)*NPARL)+1;
         idxc = (([1:ktru]-1)*NPARL)+2;
         idxr = (([1:ktru]-1)*NPARL)+3;
         idxa = (([1:ktru]-1)*NPARL)+4;
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
      if(ktru > 0)
         prof(1:ktru,1) = cumsum(mtrutmp(ifile2,idxh(1:end-1)),2);
         prof(1:ktru,1) = mtrutmp(ifile2,idxh(1:end-1));
         prof(ktru+1,1) = prof(ktru,1)+mtrutmp(ifile2,idxh(end));
      else
         prof(1,1) = zmx;
      end
      prof(:,2) = mtrutmp(ifile2,idxc);
      prof(:,3) = mtrutmp(ifile2,idxr);
      prof(:,4) = mtrutmp(ifile2,idxa);

      ctr(ifile2,:) = prof(1,2);
      rtr(ifile2,:) = prof(1,3);
      atr(ifile2,:) = prof(1,4);
      for ilay=2:ktru+1  %% ktru is # layers of current model
         idxz = round(prof(ilay-1,1)/dz);
         htr(idxz)     = htr(idxz) + 1;
         ctr(ifile2,idxz:end) = prof(ilay,2);
         rtr(ifile2,idxz:end) = prof(ilay,3);
         atr(ifile2,idxz:end) = prof(ilay,4);
      end;
      clear idxh idxc idxr idxa prof;
   end;

end; %% End loop over files

%save track_datamap.mat refmap dobsmap angobsmap nf;
%save track.mat Ncgl Nrgl Nagl ctr rtr c_max_sm r_max_sm a_max_sm c_mean_sm r_mean_sm a_mean_sm c_mead_sm r_mead_sm a_mead_sm z hgl kgl limgl logLgl logLlimgl sigmagl sigmalimgl nfc nfr nfa clim rlim alim;

else
   load track.mat
end;

figw = 12;
figh = 6.6;
%figw = 12;
%figh = 5;

ifile2 = istart - 1;
for ifile=ipst:ipskip:ipend;
   ifile2 = ifile2 + 1;
   filename = files(ifile2).name;
   disp(filename);
   filebase    = strrep(filename,'particles.mat','');
   repfile     = strrep(filename,'_particles.mat','_replica.mat');
   datafile    = strrep(filename,'_particles.mat','.txt');
   plotext1    = 'fig';
   plotext2    = 'png';
   plotfile10  = strcat(filebase,'track_uncertainty.');
   plotfile11  = strcat(filebase,'track_uncertainty_data.');
   load(repfile);
   if(isyn == 1)
     %% P02:
     ktru = track_env(ifile2,1);
     NFPTRU = ktru*4+3;
     sigtru = track_env(ifile2,NFPTRU+2);
     mtru = track_env(ifile2,2:NFPTRU+1);
     mtrutmp(ifile2,1:NFPTRU) = mtru;
   end;
   if(isyn == 1);
      NFPTRU = (track_env(ifile2,1)+1)*4-1;
   end;
   tmp = dlmread(datafile);
   z_t    = tmp(1,1);
   cw     = tmp(2,1);
   rw     = tmp(3,1);
   hmax   = tmp(4,1)+tmp(4,1)/10.;
   hmax2  = hmax;

   loc = [0.045, 0.71, 0.045,  0.14, 0.425, 0.71;...
          0.77, 0.77,   0.08,   0.08,  0.08,    0.08];
   spw1 = 0.09;
   spw2 = 0.28;
   spw3 = 0.64;
   spw4 = 0.28;
   sph1= 0.20;
   sph2= 0.60;

%   fig10=figure('visible','off');
%%   fig10=figure();
%   set(fig10,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
%   set(fig10, 'renderer', 'painters')
%   h1 = subplot('Position',[loc(1,1) loc(2,1) spw3 sph1]);
%   hold on; box off;
%   h2 = subplot('Position',[loc(1,2) loc(2,2) spw4 sph1]);
%   hold on; box off;
%   h3 = subplot('Position',[loc(1,3) loc(2,3) spw1 sph2]);
%   hold on; box off;
%   h4 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph2]);
%   hold on; box off;
%   h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph2]);
%   hold on; box off;
%   h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph2]);
%   hold on; box off;
%
%   %%
%   %% MAP ensemble velocity track profile
%   %%
%   subplot(h1);hold on;box on;
%   set(gca,'FontSize',12);
%   imagesc([ipst:ipskip:ipend],z,c_mead_sm);shading flat;
%   plot([ifile-1+ipskip/2, ifile-1+ipskip/2],[0, z(end)],'k','LineWidth',1.5);
%   plot([ifile-1+ipskip/2, ifile-1+ipskip/2],[0, z(end)],'--w','LineWidth',1.5);
%   xlabel('Ping No.');
%   ylabel('Depth (m)');
%   set(gca,'Layer','top','YDir','reverse');box on;
%   set(gca,'XLim',[ipst ipend],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
%   c1=colorbar;ylabel(c1,'Velocity (m/s)')
%
%   subplot(h2);hold on;box on;
%   set(gca,'FontSize',12);
%   [xx,yy]=stairs(limgl,kgl(ifile2,:),'k');
%   patch(xx,yy,[0.8,0.8,0.8]);
%   stairs(limgl,kgl(ifile2,:),'k');
%   xlabel('No. interfaces in partition');
%   ylabel('Probability');
%   set(gca,'XLim',[-0.5 kmx+0.5],'YLim',[0. 1.0]);
%
%   subplot(h3);hold on;box on;
%   if(imarg == 1)
%      [xx,yy]=stairs(hgl(ifile2,:),z,'-k');
%      patch(xx,yy,[0.8,0.8,0.8]);
%      [xx,yy]=stairs(hgl(ifile2,:),z,'-k');
%      set(gca,'Fontsize',12,'YLim',[0 hmax2]);
%      set(gca,'YDir','reverse','XLim',[0 0.1]);
%      set(gca,'YTickLabel',[0 1 2 3 4],'YTick',[0 1 2 3 4]);
%      xlabel('Interface probability');
%      ylabel('Depth (m)');
%      box on;
%   end;
%
%   subplot(h4)
%   set(gca,'FontSize',12);
%   if(imarg == 1)
%      imagesc(clim,z,Ncgl(:,:,ifile2));shading flat;
%      if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
%      if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
%   end;
%   %surf(clim,z,Nc);shading flat;
%   set(h4,'layer','top')
%   set(gca,'Fontsize',12,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
%   set(gca,'YDir','reverse');
%   xlabel('Velocity (m/s)');
%   if(isyn == 1)
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,'-k',1,0);
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,'--w',1,0);
%   end;
%   if(imap == 1)
%      plprof(mmap,hmax,'--k',1,0);
%   end;
%   if(imax == 1)
%     for i=1:length(z);
%        [cmx,j] = max(Nc(i,:));
%        c_max(i) = clim(j);
%     end;
%     plot(c_max,z,'w');
%   end;
%   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
%   set(gca,'XTickLabel',[1500 1600 1700]);
%   set(gca,'XTick',[1500 1600 1700]);
%   box on;
%
%   subplot(h5)
%   set(gca,'FontSize',12);
%   if(imarg == 1)
%      imagesc(rlim,z,Nrgl(:,:,ifile2));shading flat;
%      if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
%      if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
%   end;
%   %surf(rlim,z,Nr);shading flat;
%   set(h5,'layer','top')
%   set(gca,'Fontsize',12,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
%   set(gca,'YDir','reverse');
%   xlabel('Density (g/ccm)');
%   if(isyn == 1)
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,'-k',2,0);
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,'--w',2,0);
%   end
%   if(imap == 1)
%      plprof(mmap,hmax,'--k',2,0);
%   end
%   if(imax == 1)
%     for i=1:length(z);
%        [rmx,j] = max(Nr(i,:));
%        r_max(i) = rlim(j);
%     end;
%     plot(r_max,z,'w');
%   end;
%   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
%   set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
%   set(gca,'XTick',[1.4 1.6 1.8 2.0]);
%   box on;
%
%   subplot(h6)
%   set(gca,'FontSize',12);
%   if(imarg == 1)
%      imagesc(alim,z,Nagl(:,:,ifile2));shading flat;
%      if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
%      if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
%   end;
%   %surf(alim,z,Na);shading flat;
%   set(h6,'layer','top')
%   set(gca,'Fontsize',12,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
%   set(gca,'YDir','reverse');
%   xlabel('Attenuation');
%   if(isyn == 1)
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,'-k',3,0);
%      plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,'--w',3,0);
%   end;
%   if(imap == 1)
%      plprof(mmap,hmax2,'--k',3,0);
%   end;
%   if(imax == 1)
%     for i=1:length(z);
%        [amx,j] = max(Na(i,:));
%        a_max(i) = alim(j);
%     end;
%     plot(a_max,z,'w');
%   end;
%   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
%   %cmap = colormap(flipud(gray));
%   cmap = colormap(jet);
%   %cmap(1,:) = [1 1 1];
%   colormap(cmap);
%   %colorbar('peer',h1,'location','WestOutside');
%   box on;
%
%%   saveas(fig10,strcat(plotfile10,plotext2),'png');
%   print(fig10,'-painters','-r250',strcat(plotfile10,plotext2),'-dpng');
%%   print(fig10,'-painters','-r250',strcat(plotfile10,plotext3),'-depsc');
%   close(fig10);

   figw = 11;
   figh = 7.0;
          
   if(isyn == 1)
     loc2 = [0.045, 0.71, 0.045, 0.14, 0.425, 0.71, 0.0455, ...
             0.2830, 0.5205, 0.7580;...
             0.81, 0.81, 0.29, 0.29, 0.29, 0.29, 0.065, 0.065, 0.065, 0.065];
   else
   loc2 = [0.045, 0.71, 0.045, 0.14, 0.425, 0.71, 0.0455, ...
           0.40125, 0.7580;...
           0.81, 0.81, 0.29, 0.29, 0.29, 0.29, 0.065, 0.065, 0.065];
   end
   spw1 = 0.09;
   spw2 = 0.28;
   spw3 = 0.64;
   spw4 = 0.28;
   spw5 = 0.2325;
   sph1= 0.18;
   sph2= 0.45;
   sph3= 0.15;
   
   fig11=figure('visible','off');
%   fig11=figure();
   set(fig11,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   set(fig11, 'renderer', 'painters')
   h1 = subplot('Position',[loc2(1,1) loc2(2,1) spw3 sph1]);
   hold on; box off;
   h2 = subplot('Position',[loc2(1,2) loc2(2,2) spw4 sph1]);
   hold on; box off;
   h3 = subplot('Position',[loc2(1,3) loc2(2,3) spw1 sph2]);
   hold on; box off;
   h4 = subplot('Position',[loc2(1,4) loc2(2,4) spw2 sph2]);
   hold on; box off;
   h5 = subplot('Position',[loc2(1,5) loc2(2,5) spw2 sph2]);
   hold on; box off;
   h6 = subplot('Position',[loc2(1,6) loc2(2,6) spw2 sph2]);
   hold on; box off;
   h7 = subplot('Position',[loc2(1,7) loc2(2,7) spw5 sph3]);
   hold on; box off;
   h8 = subplot('Position',[loc2(1,8) loc2(2,8) spw5 sph3]);
   hold on; box off;
   h9 = subplot('Position',[loc2(1,9) loc2(2,9) spw5 sph3]);
   hold on; box off;
   if(isyn == 1)
     h10 = subplot('Position',[loc2(1,10) loc2(2,10) spw5 sph3]);
     hold on; box off;
   end;

   %%
   %% MAP ensemble velocity track profile
   %%
   subplot(h1);hold on;box on;
   set(gca,'FontSize',12);
   imagesc([ipst:ipskip:ipend],z,c_mead_sm);shading flat;
   plot([ifile, ifile],[0, z(end)],'k','LineWidth',1.5);
   plot([ifile, ifile],[0, z(end)],'--w','LineWidth',1.5);
   xlabel('Ping No.');
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[ipst ipend],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
   c1=colorbar;ylabel(c1,'Velocity (m/s)')

   subplot(h2);hold on;box on;
   set(gca,'FontSize',12);
   [xx,yy]=stairs(limgl,kgl(ifile2,:),'k');
   patch(xx,yy,[0.8,0.8,0.8]);
   stairs(limgl,kgl(ifile2,:),'k');
   if(isyn == 1)
     plot([track_env(ifile2,1), track_env(ifile2,1)],[0, 1],'k','LineWidth',1.5);
     plot([track_env(ifile2,1), track_env(ifile2,1)],[0, 1],'--w','LineWidth',1.5);
   end;
   xlabel('No. interfaces in partition');
   ylabel('Probability');
   set(gca,'YTickLabel',[0 1],'YTick',[0 1]);
   set(gca,'XLim',[-0.5 kmx+0.5],'YLim',[0. 1.0]);

   subplot(h3);hold on;box on;
   if(imarg == 1)
      [xx,yy]=stairs(hgl(ifile2,:),z,'-k');
      patch(xx,yy,[0.8,0.8,0.8]);
      [xx,yy]=stairs(hgl(ifile2,:),z,'-k');
      set(gca,'Fontsize',12,'YLim',[0 hmax2]);
      set(gca,'YDir','reverse','XLim',[0 0.1]);
      set(gca,'YTickLabel',[0 1 2 3 4],'YTick',[0 1 2 3 4]);
      xlabel('Interface probability');
      ylabel('Depth (m)');
      box on;
   end;

   subplot(h4)
   set(gca,'FontSize',12);
   if(imarg == 1)
      imagesc(clim,z,Ncgl(:,:,ifile2));shading flat;
      if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
      if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
   end;
   %surf(clim,z,Nc);shading flat;
   set(h4,'layer','top')
   set(gca,'Fontsize',12,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Velocity (m/s)');
   if(isyn == 1)
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',1,0);
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',1,0);
   end;
%   if(imap == 1)
%      plprof(mmap,hmax,'--k',1,0);
%   end;
   if(imax == 1)
     for i=1:length(z);
        [cmx,j] = max(Nc(i,:));
        c_max(i) = clim(j);
     end;
     plot(c_max,z,'w');
   end;
   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
   set(gca,'XTickLabel',[1500 1600 1700]);
   set(gca,'XTick',[1500 1600 1700]);
   box on;

   subplot(h5)
   set(gca,'FontSize',12);
   if(imarg == 1)
      imagesc(rlim,z,Nrgl(:,:,ifile2));shading flat;
      if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
      if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
   end;
   %surf(rlim,z,Nr);shading flat;
   set(h5,'layer','top')
   set(gca,'Fontsize',12,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Density (g/ccm)');
   if(isyn == 1)
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',2,0);
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',2,0);
   end
%   if(imap == 1)
%      plprof(mmap,hmax,'--k',2,0);
%   end
   if(imax == 1)
     for i=1:length(z);
        [rmx,j] = max(Nr(i,:));
        r_max(i) = rlim(j);
     end;
     plot(r_max,z,'w');
   end;
   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
   set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
   set(gca,'XTick',[1.4 1.6 1.8 2.0]);
   box on;

   subplot(h6)
   set(gca,'FontSize',12);
   if(imarg == 1)
      imagesc(alim,z,Nagl(:,:,ifile2));shading flat;
      if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
      if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
   end;
   %surf(alim,z,Na);shading flat;
   set(h6,'layer','top')
   set(gca,'Fontsize',12,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Attenuation');
   if(isyn == 1)
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'-k',3,0);
      plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'--w',3,0);
   end;
%   if(imap == 1)
%      plprof(mmap,hmax2,'--k',3,0);
%   end;
   if(imax == 1)
     for i=1:length(z);
        [amx,j] = max(Na(i,:));
        a_max(i) = alim(j);
     end;
     plot(a_max,z,'w');
   end;
   set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
   %cmap = colormap(flipud(gray));
   cmap = colormap(jet);
   %cmap(1,:) = [1 1 1];
   colormap(cmap);
   %colorbar('peer',h1,'location','WestOutside');
   box on;

   nx = 4;
   ny = ceil(length(bands)/nx);
   ii = 1;
   for i=1:length(bands); 
      if(i==1);subplot(h7);hold on;box on;end;
      if(i==2);subplot(h8);hold on;box on;end;
      if(i==3);subplot(h9);hold on;box on;end;
      if(i==4);subplot(h10);hold on;box on;end;
      set(gca,'FontSize',12);
      refmean = mean(ref(:,:,:),3);
      set(gca,'XLim',[angobs(1)-1 angobs(end)+1]);
      set(gca,'FontSize',12);
      set(gca,'LineWidth',1);
%      plot(rep(length(bands)+1,:),rep(i,:),'-r','LineWidth',2);
      for j=1:NDAVE;
         plot(angobs,ref(:,i,j),':r');
      end;
%      plot(angobs,refmean(:,i),'-r');
%      plot(angobs,ref(:,1,1).*stat(1).rnan,'--b');
      plot(angobs,dobs(:,i),'xk');
      for j=1:length(angobs);
         y = ref(j,i,:);
         [nf2(j,:)] = hpd(y,100,95);
      end;
      plot(angobs,nf2(:,1),'--k');
      plot(angobs,nf2(:,2),'--k');
      if(ibl == 0)
%         if(max(max(dobs)< 0.64)
            ylim  = [0,0.75];
            ytick = [0,0.2,0.4,0.6];
            text(52,0.65,[num2str(bands(i)) ' Hz'],'FontSize',12)
%         else
%            ylim  = [0,1.01];
%            ytick = [0,0.3,0.6,0.9];
%            text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12)
%         end
         set(gca,'YLim',ylim,'XLim',[27 68]);
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
          set(gca,'XTickLabel',[30 45 60],'XTick',[30 45 60]);
          xlabel('Angle (deg.)');
      else
          set(gca,'XTickLabel',[],'XTick',[30 45 60]);
      end
   end;
   clear nf2;

%   saveas(fig11,strcat(plotfile11,plotext2),'png');
   print(fig11,'-painters','-r250',strcat(plotfile11,plotext2),'-dpng');
%   print(fig11,'-painters','-r250',strcat(plotfile11,plotext3),'-depsc');
   close(fig11);

end; %% End loop over files


return;
