%function [] = auv_plot_rjhist_track(istart);
clear all;
close all;
istart = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 1; %% Normalize profile marginals line by line to unit area
isecondfig = 0;
isyn    = 0;
i_vref  = 1;
ibucking= 3;
ivgsnosh= 0;
ishear  = 0;
ivoro   = 0;
imap    = 0;    
imead   = 0;
imean   = 0;
imax    = 0;      %%
isave   = 0;      %%
iplot   = 0;
idatfit = 0;      %% plot data misfits?
isd     = 1;      %% sampled over sigma?
iar     = 1;
ivpvs   = 0;
icore   = 0;
itrackmat = 0;    %% 1 -> Already pre-computes track.mat
FNT = 14;

%files=dir('*_particles.mat');
files=dir('*2513_sample.mat');

kmx = 15;
if(ibucking == 0);
  NPL = 4;
  NPV = 4;
  if(ishear == 1);NPL = 6;NPV=6;end;
elseif(ibucking == 3 & ivgsnosh == 0);
  NPL = 5;
  NPV = 5;
  if(ishear == 1);NPL = 6;NPV=6;end;
else;
  NPL = 4;
end;
if(ivgsnosh == 1)
  NMISC = 5;
else;
  NMISC = 4;
end;
hmax = 7.5
hmax2 = 4.0;

pinglist = load('ping_list.txt');
NPING = pinglist(1);
%NPING = 5;
dping = 300;      %% Interval for ping axis label
pinglist(1)=[];
ipst = pinglist(1);

%% Malta
bands = [988., 1113., 1288., 1913., 2263., 2513.];
%% Simulation track spherical:
%bands = [975, 1100, 1250, 2100, 2400, 2700];
%% Simulation track plane:
%bands = [1000, 1200, 2000, 2400];

%% Bucking simulation:
if(ibucking == 0)
  %% Malta 
  %pmin = [1450 1.15 -3.]';
  %pmax = [1930 2.0 -.3]';
  %% Simulation track spherical:
  %pmin = [1450 1.2 -3]';
  %pmax = [1750 2.2 0.01]';
  %% Simulation track plane:
  %pmin = [1450 1.2 0.001]';
  %pmax = [1750 2.2 1]';
  pmin = [1400 1.3 -3.0]';
  pmax = [2000 2.4 -0.3]';
  if(isyn == 1)
    pmin = [1450 1.20 0  ]';
    pmax = [1750 2.20 1.0]';
    %pmin = [1400 1.2 0    log(sqrt(2))  0]';
    %pmax = [1700 2.0 0.1  log(100)  5]';
  end
elseif(ibucking == 3)
  if(ivgsnosh == 0)
    %% VGS AUV:
    pmin = [0.20 6.9 -1.4 -4.5]';
    pmax = [0.95 8.6 -0.4 -1.3]';
  else;
    pmin = [0.20 log(2.e7) -10]';
    pmax = [0.91 log(8.e8) -3]';
  end;
end;
if(ibucking > 0)
%  miscmin = [3.55e10, 2.3400e9,2720.,1028.]';
%  miscmax = [3.65e10, 2.3500e9,2760.,1032.]';
%  miscmin = [3.25e10, 2.0000e9,2500.,1010.]';
%  miscmax = [3.95e10, 2.5000e9,2800.,1040.]';
  miscmin = [3.00e10, 2.0000e9,2400.,1000.]';
  miscmax = [4.00e10, 2.5000e9,2800.,1050.]';
end;
%   NDAVE   = ceil(NPROF/10);
%NDAVE   = 500;
thinstep = 1;

NBAND = length(bands);
NANG = 32;
NDAVE = 1200;
NSD = NBAND;

sigmamin = 0.00;
sigmamax = 0.40;
thinstep = 1;
NZI = 400;
NZ  = 250;
NC  = 200;
NR  = 200;
NA  = 200;
NA2 = 100;

logLmin = 320;
logLmax = 480;
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NBAND);
   armax = 1.*ones(1,NBAND);
else
   order = 1;
end;

if(isyn == 1);
  track_env = dlmread('track_environment_z_logalf.dat');
  logLtru = dlmread('logL_true.txt');

  logLmin = min(logLtru) - 80;
  logLmax = max(logLtru) + 20;
end;

htr = zeros(1,NZ);
ctr = zeros(NPING,NZ);
rtr = zeros(NPING,NZ);
atr = zeros(NPING,NZ);

if(isyn == 1);
   mtrutmp = zeros(size(track_env));
end;

ifile2 = istart-1;

if(itrackmat == 0)
  ref2 = zeros(NANG,NBAND,NDAVE,NPING);
  for ifile=1:NPING;

       ifile2 = ifile2 + 1;
       if(ifile2 > pinglist(end));exit;end;
       filename = files(ifile2).name;
       disp(filename);
    %   filebase    = strrep(filename,'_particles.mat','');
       filebase    = strrep(filename,'_sample.mat','');
       datafile    = strcat(filebase,'.txt');
    %   repfile     = strcat(filebase,'_replicas.mat');
       repfile     = strcat(filebase,'_rep_ens.mat');
    %   vreffile    = strcat(filebase,'_sample.mat','_vel_ref.txt');
       vreffile    = 'p0004_pave010_0988_2513_vel_ref.txt';
       plotext1    = 'fig';
       plotext2    = 'png';
       plotext3    = 'eps';
       plotfile1   = strcat(filebase,'_afdep.');
       plotfile2   = strcat(filebase,'_layer_thickness_marginal.');
       plotfile3   = strcat(filebase,'_sigma.');
       plotfile4   = strcat(filebase,'_number_interfaces.');
       plotfile5   = strcat(filebase,'_logL.');
       plotfile6   = strcat(filebase,'_transdim_marg_prof.');
       plotfile7   = strcat(filebase,'_data.');
       plotfile8   = strcat(filebase,'_axx.');
       plotfile9   = strcat(filebase,'_reshist.');
       corefile    = 'core.mat';
    %   clear A m k logL dobs angobs sd;
       clear A m k logL sd;
       load(filename);
       load(filename);
       burnin = length(A)/2;
       A(1:burnin,:)=[];
       if(ivoro == 0);B=A;end;
       datafile
       tmp = dlmread(datafile);
       z_t    = tmp(1,1);
       cw     = tmp(2,1);
       rw     = tmp(3,1);
       hmax   = tmp(4,1)+tmp(4,1)/10.;
    %   hmax2  = hmax;
       dobs(:,:,ifile2)   = tmp(5:length(bands)+4,:)';
       angobs = tmp(length(bands)+5,:);
    %   rex(:,:,ifile2)   = tmp(6+length(bands):length(bands)+length(bands)+5,:)';
    %   angobs = [0:0.25:31];
       nang   = length(angobs);
       nband  = length(bands);
    %   dobs = ones(nang,nband);

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
      NPROF = length(A);

      %% Read in ref velocity

      NFP = (k*NPL)+NPL-1;
      m = A(:,5:end-5);
      %% Read in ref velocity
      if(i_vref == 1);
        vel_ref = dlmread(vreffile);
        vel_ref(1,:) = [];
        for ismp=1:length(m);
          for ik=1:k(ismp)+1;
            if(ik <= k(ismp));
              ipar = (ik-1)*NPL+1;
              z = m(ismp,ipar);
            else;
              ipar = (ik-1)*NPL;
              z = m(ismp,ipar-NPL+1);
            end;
            if(ibucking == 0);
              [vref,rref]=refl_getref(z,vel_ref);
              m(ismp,ipar+1) = vref + m(ismp,ipar+1);
              m(ismp,ipar+2) = rref + m(ismp,ipar+2);
            else;
              [pref]=refl_getrefpor(z,vel_ref);
              m(ismp,ipar+1) = pref + m(ismp,ipar+1);
            end;
          end;
        end;
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
      if(iar == 1)
          (size(B,2)-6-order*NBAND-NMISC+1:size(B,2)-NMISC-6)
          alpha = B(:,end-6-(order*NBAND)-NMISC+1:end-NMISC-6);
      end
      if(ibucking > 0)
        %(size(B,2)-6-NMISC-order*NBAND+1:size(B,2)-6-NMISC)
        bulkrho = B(:,end-6-(NMISC)+1:end-6);
      end

      logL = A(:,1);
      if(imap == 1);
        [logLmap,jmap] = max(logL);
        kmap = k(jmap);
        NFPmap = NFP(jmap);
        mmap = m(jmap,1:NFPmap);
      end;

      dep = dep(:);dep(find(dep==0))=[];
      h = h(:);h(find(h==0))=[];
      [a_enos,b_enos]=hist(h,200);a_enos=a_enos/trapz(b_enos,a_enos);
      a_enos=[0,a_enos,0];
      b_enos=[b_enos(1),b_enos,b_enos(end)];
      fig2 = figure('visible','off');
      [xx,yy]=stairs(b_enos,a_enos,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      if(isd == 1)
        sd(:,1:NSD) = A(:,end-6-(order*NBAND)-NSD-NMISC+1:end-6-(order*NBAND)-NMISC);
        disp('Done getting sigma.');
      end
      for i = 1:size(A,1);
        idxh = [0:k(i)-1]*4+1;
        h(i) = sum(m(i,idxh));
        clear idxh;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%
      %% SIGMA PLOT
      %%
      if(iplot == 1);
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
          h2 = subplot('Position',[loc(1,i+NSD) loc(2,i+NSD) spw sph]);
          hold on; box off;
          set(gca,'FontSize',FNT);
          subplot(h1);hold on;box on;
          plot([1:thinstep:length(sd(:,i))],sd(1:thinstep:end,i),'k');
          if(i == 1);ylabel('Data error standard deviation');end;
          xlabel('Particle no.');
          set(gca,'XLim',[0 length(sd(:,i))])
          set(gca,'YLim',[sigmamin sigmamax],'YTick',[0.0 0.1 0.4])
          if(i > 1);set(gca,'YTickLabel',[]);end;

          subplot(h2);set(gca,'Layer','top');hold on;
    %     [n,lim]=hist(sd(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
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
          set(gca,'YLim',[sigmamin sigmamax],'YTick',[0.0 0.1 0.4])
          if(i > 1);set(gca,'YTickLabel',[]);end;
          box on;
        end;
      end;
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %% ALPHA PLOT
    %%
    if(iplot == 1);
    if(iar == 1)
       figw = 12;
       figh = 8;
       fig8=figure('visible','off');
       set(fig8,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
       nx = NBAND;
       ny = 3;
       xim = 0.01;
       yim = 0.01;
       xymarg = [0.1 0.04 0.04 0.1];
       [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
       j = 1;
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
         lim = lim - (lim(3)-lim(2))/2;
         [xx,yy]=stairs(n,lim,'k');
         patch(xx,yy,[0.8,0.8,0.8]);
         stairs(n,lim,'k');
         clear n lim;
         if(j==order);xlabel('Probability');end;
         if(i==1);ylabel('AR coefficient');end;
         if(i>1);set(gca,'YTickLabel',[]);end;
         if(j<order);set(gca,'XTickLabel',[]);end;
         set(gca,'XLim',[0 0.22],'XAxisLocation','top');
         set(gca,'YLim',[-.6 1]);

         idxar = ones(size(alpha(:,(i-1)*order+j)));
         idx = find(alpha(:,(i-1)*order+j) < -0.5);
         idxar(idx) = 0;
         subplot(h2);hold on;box on;
         lim = [-1:1:2];
         [n]=hist(idxar,lim);%n = [0, n, 0];lim = [lim(1) lim lim(end)];
         n = n/sum(n);
         lim = lim - (lim(3)-lim(2))/2;
         [xx,yy]=stairs(lim,n,'k');
         patch(xx,yy,[0.8,0.8,0.8]);
         stairs(lim,n,'k');
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
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %% Bulk Moduli and Density (pore fluid/grain material) PLOT
    %%
    if(iplot == 1);
      if(ibucking > 0)
        fig33=figure('visible','off');
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
          lim = lim - (lim(3)-lim(2))/2;
          [xx,yy]=stairs(lim,n,'k');
          patch(xx,yy,[0.8,0.8,0.8]);
          stairs(lim,n,'k');
          if(isyn == 1);
            plot([mtrumisc(i) mtrumisc(i)],[0 .04],'-w');
            plot([mtrumisc(i) mtrumisc(i)],[0 .04],'--k');
          end;
          clear n lim;
          set(gca,'XLim',[miscmin(i) miscmax(i)])
          if(i == 1);xlabel('Kg (mineral grains) [Pa]');end;
          if(i == 2);xlabel('Kf (interstitial fluid) [Pa]');end;
          if(i == 3);xlabel('Density grains (kg/m^3)');end;
          if(i == 4);xlabel('Density fluid (kg/m^3)');end;
          if(i == 5);xlabel('Strain hardening');end;
    %      set(gca,'YLim',[0.0 0.06],'YTick',[0.0 0.02 0.04 0.06])
          set(gca,'YTickLabel',[]);
          box on;
        end;
      end;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %% K PLOT
    %%
    if(iplot == 1);
      fig4=figure('visible','off');
      subplot(2,1,1);hold on;box on;
      set(gca,'FontSize',FNT);
      plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
      ylabel('No. interfaces');
      xlabel('Sample No.');
    %  set(gca,'XLim',[0 length(k)])
      set(gca,'YLim',[-0.5 kmx+0.5])
      if(isyn == 1)
       plot([0,length(k)],[ktru,ktru],'k','LineWidth',1.5);
       plot([0,length(k)],[ktru,ktru],'--w','LineWidth',1.5);
      end;

      subplot(2,1,2);hold on;box on;
      set(gca,'FontSize',FNT);
      [n,lim]=hist(k,[0:kmx]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      lim = lim-0.5;
      kgl(ifile2,:) = n;
      limgl = lim;
      [xx,yy]=stairs(lim,n,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(lim,n,'k');
      clear n lim;
      xlabel('No. interfaces');
      ylabel('Probability');
      set(gca,'XLim',[-0.5 kmx+0.5],'YLim',[0. 1.0]);
      if(isyn == 1)
        plot([ktru,ktru],[0,1],'k','LineWidth',1.5);
        plot([ktru,ktru],[0,1],'--w','LineWidth',1.5);
      end;
    end;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%
     %% logL PLOT
     %%
     if(iplot == 1)
       fig5=figure('visible','off');
       subplot(1,2,1);hold on;box on;
       set(gca,'FontSize',FNT);
       plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
       ylabel('log Likelihood');
       xlabel('Sample No.');
       set(gca,'XLim',[0 length(logL)])
       set(gca,'YLim',[logLmin logLmax])
       %set(gca,'YLim',[150 400])
       if(isyn == 1)
         plot([0,length(logL)],[logLtru(ifile2),logLtru(ifile2)],'k','LineWidth',1.5);
         plot([0,length(logL)],[logLtru(ifile2),logLtru(ifile2)],'--w','LineWidth',1.5);
       end;   

       subplot(1,2,2);hold on;box on;
       set(gca,'FontSize',FNT);
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
       set(gca,'XLim',[0 1])
       %set(gca,'YLim',[150 400])
       if(isyn == 1)
         xlim = get(gca,'XLim');
    %     plot([0,xlim(2)],[logLtru(ifile2),logLtru(ifile2)],'k','LineWidth',1.5);
    %     plot([0,xlim(2)],[logLtru(ifile2),logLtru(ifile2)],'--w','LineWidth',1.5);
         plot([0,10],[logLtru(ifile2),logLtru(ifile2)],'k','LineWidth',1.5);
         plot([0,10],[logLtru(ifile2),logLtru(ifile2)],'--w','LineWidth',1.5);
       end;
     end;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%
     %% COMPUTE PROFILE MARGINALS
     %%
     if(imarg == 1);
       NZ = 500;
       NZI= 500;
       nsmooth = ceil(NZ/80.);
       NC = 300;
       NR = 300;
       NA = 300;
       NA2 = 100;
       clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
       rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
       alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         tlim = pmin(4)+cumsum((pmax(4)-pmin(4))/NA2*ones(1,NA2));
         if(ishear == 1);
         tlim2 = pmin(5)+cumsum((pmax(5)-pmin(5))/NR*ones(1,NR));end;
       end;

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
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         t = zeros(NPROF,NZ);
         if(ishear == 1);t2 = zeros(NPROF,NZ);end;
       end;
       disp('Sample size: '),disp(size(m))
       for iprof = 1:NPROF

         if(rem(iprof,10000)==0)
            fprintf(1,'%8i  ',iprof)
         end
         clear idxh idxc idxr idxa prof;
         %% Find index for current model
         if(k(iprof) > 0)
            idxh = (([1:k(iprof)]-1)*NPL)+1;
            idxc = (([1:k(iprof)]-1)*NPL)+2;
            idxr = (([1:k(iprof)]-1)*NPL)+3;
            idxa = (([1:k(iprof)]-1)*NPL)+4;
            if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
              idxt = (([1:k(iprof)]-1)*NPL)+5;
              if(ishear == 1);
                idxt2 = (([1:k(iprof)]-1)*NPL)+6;
              end;
            end;
            idxh = [idxh idxh(end)];
            idxc = [idxc idxc(end)+NPL-1];
            idxr = [idxr idxr(end)+NPL-1];
            idxa = [idxa idxa(end)+NPL-1];
            if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
              idxt = [idxt idxt(end)+NPL-1];
              if(ishear == 1);
                idxt2 = [idxt2 idxt2(end)+NPL-1];
              end;
            end;
         else
            idxh = [];
            idxc = [1];
            idxr = [2];
            idxa = [3];
            if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
              idxt = [4];
              if(ishear == 1);idxt2 = [5];end;
            end;
         end

         %% Compute the profile for current model
         if(k(iprof) > 0)
           prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
           prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
           prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
           prof(:,2) = m(iprof,idxc);
           prof(:,3) = m(iprof,idxr);
           prof(:,4) = m(iprof,idxa);
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
               prof(:,5) = m(iprof,idxt);
             if(ivpvs == 1 & ibucking == 0);
               prof(:,5) = m(iprof,idxc)./exp(m(iprof,idxt));
             end;
             if(ishear == 1);prof(:,6) = m(iprof,idxt2);end;
           end;

           for ilay=2:k(iprof)+1  %% k is # layers of current model
              idxzi = round(prof(ilay-1,1)/dzi);
              if(idxzi == 0);idxzi = 1;end;
              h(idxzi) = h(idxzi) + 1;
           end;
           c(iprof,:) = prof(1,2);
           r(iprof,:) = prof(1,3);
           a(iprof,:) = prof(1,4);
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
             t(iprof,:) = prof(1,5);
             if(ishear == 1);t2(iprof,:) = prof(1,6);end;
           end;
           for ilay=2:k(iprof)+1  %% k is # layers of current model
              idxz = round(prof(ilay-1,1)/dz);
              if(idxz == 0)idxz = 1;end;
              c(iprof,idxz:end) = prof(ilay,2);
              r(iprof,idxz:end) = prof(ilay,3);
              a(iprof,idxz:end) = prof(ilay,4);
              if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
                t(iprof,idxz:end) = prof(ilay,5);
                if(ishear == 1);t2(iprof,idxz:end) = prof(ilay,6);end;
              end;
           end;

         else
           prof(:,2) = m(iprof,idxc);
           prof(:,3) = m(iprof,idxr);
           prof(:,4) = m(iprof,idxa);
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
               prof(:,5) = m(iprof,idxt);
             if(ivpvs == 1 & ibucking == 0);
               prof(:,5) = m(iprof,idxc)./exp(m(iprof,idxt));
             end;
             if(ishear == 1);prof(:,6) = m(iprof,idxt2);end;
           end;
           c(iprof,:) = prof(1,2);
           r(iprof,:) = prof(1,3);
           a(iprof,:) = prof(1,4);
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
             t(iprof,:) = prof(1,5);
             if(ishear == 1);t2(iprof,:) = prof(1,6);end;
           end;
         end
       end;
       fprintf(1,'\n')
       disp('Done with profiles.');

       %
       % Compute histograms for each depth
       %
       disp('Starting histograms...');
       for iz=1:NZ
          [Nc(iz,:),binsc] = hist(c(:,iz),clim);
          [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
          [Na(iz,:),binsa] = hist(a(:,iz),alim);
          if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
            [Nt(iz,:),binst] = hist(t(:,iz),tlim);
            if(ishear == 1);[Nt2(iz,:),binst2] = hist(t2(:,iz),tlim2);end;
          end;
          [nfc(iz,:)] = hpd(c(:,iz),100,95);
          [nfr(iz,:)] = hpd(r(:,iz),100,95);
          [nfa(iz,:)] = hpd(a(:,iz),100,95);
          if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
            [nft(iz,:)] = hpd(t(:,iz),100,95);
            if(ishear == 1);[nft2(iz,:)] = hpd(t2(:,iz),100,95);end;
          end;
          meac(iz) = median(c(:,iz));
          mear(iz) = median(r(:,iz));
          meaa(iz) = median(a(:,iz));
          if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
            meat(iz) = median(t(:,iz));
            if(ishear == 1);meat2(iz) = median(t2(:,iz));end;
          end;
       end;
       %
       % Normalize Histograms
       %
       if(inorm == 0)
         for iz=1:NZ
           Nc(iz,:) = Nc(iz,:)/NPROF;
           Nr(iz,:) = Nr(iz,:)/NPROF;
           Na(iz,:) = Na(iz,:)/NPROF;
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
             Nt(iz,:) = Nt(iz,:)/NPROF;
             if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/NPROF;end;
           end;
         end;
       elseif(inorm == 1)
          for iz=1:NZ
             Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
             Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
             Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
            if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
              Nt(iz,:) = Nt(iz,:)/trapz(binst,Nt(iz,:));
              if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/trapz(binst2,Nt2(iz,:));end;
            end;
          end;
       elseif(inorm == 2)
         for iz=1:NZ
           Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
           Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
           Na(iz,:) = Na(iz,:)/max(Na(iz,:));
           if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
             Nt(iz,:) = Nt(iz,:)/max(Nt(iz,:));
             if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/max(Nt2(iz,:));end;
           end;
         end;
       end;
       disp('Done histograms.');
           
       Ncgl(:,:,ifile2) = Nc;
       Nrgl(:,:,ifile2) = Nr;
       Nagl(:,:,ifile2) = Na;
       Ntgl(:,:,ifile2) = Nt;

       c_mean(:,ifile2) = mean(c);
       r_mean(:,ifile2) = mean(r);
       a_mean(:,ifile2) = mean(a);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         t_mean(:,ifile2) = mean(t);
       end;
       c_nf(:,:,ifile2) = nfc;
       r_nf(:,:,ifile2) = nfr;
       a_nf(:,:,ifile2) = nfa;
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         nf_t(:,:,ifile2) = nft;
       end;

       c_mead(:,ifile2) = median(c);
       r_mead(:,ifile2) = median(r);
       a_mead(:,ifile2) = median(a);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         t_mead(:,ifile2) = median(t);
       end;

       [ntmp,idxcmax] = max(Nc,[],2);
       [ntmp,idxrmax] = max(Nr,[],2);
       [ntmp,idxamax] = max(Na,[],2);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         [ntmp,idxtmax] = max(Nt,[],2);
       end;
       for iz=1:NZ
         c_max(iz,ifile2) = clim(idxcmax(iz));
         r_max(iz,ifile2) = rlim(idxrmax(iz));
         a_max(iz,ifile2) = alim(idxamax(iz));
         if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
           t_max(iz,ifile2) = tlim(idxtmax(iz));
         end;
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
         if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));t_max_sm(i) = sqrt(mean(t_max(NAVE1:NAVE2).^2));end;
       end
     end; % end imarg
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%
     %% PLOT PROFILE MARGINALS
     %%
     if(iplot == 1);
       opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
       %[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

       if((ibucking==3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
       spw1 = 0.06;
       spw2 = 0.2;
       sph = 0.82;
       loc(:,1) = [0.06; 0.14];
       loc(:,2) = [loc(1,1)+spw1+.01; 0.14];
       loc(:,3) = [loc(1,2)+spw2+.01; 0.14];
       loc(:,4) = [loc(1,3)+spw2+.01; 0.14];
       loc(:,5) = [loc(1,4)+spw2+.01; 0.14];  
       if(ishear==1);
        loc = [0.075, 0.163, 0.331, 0.499, 0.667, 0.835;...
               0.14,   0.14,  0.14,  0.14,  0.14, 0.14];
        spw1 = 0.08;
        spw2 = 0.16;
        sph  = 0.82;
      end;

        else
          loc = [0.08,  0.18, 0.4533, 0.7266;...
                 0.14, 0.14, 0.14, 0.14];
          spw1 = 0.09;
          spw2 = 0.2633;
          sph = 0.82;
        end;

        fig6 = figure('visible','off');hold on; box on;
        figw = 18;
        figh = 10;
        fig6.PaperUnits = 'inches';
        fig6.PaperPosition = [0 0 figw figh];
        set(fig6, 'renderer', 'painters')
        h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
        hold on; box off;
        h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
        hold on; box off;
        h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
        hold on; box off;
        h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
        hold on; box off;
        if((ibucking == 3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
          h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]);
          hold on; box off;
          if(ishear == 1);
            h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph]);
            hold on; box off;
          end;
        end;

        subplot(h4);hold on;box on;
        if(imarg == 1)
           h = h/sum(h);
           hgl(ifile2,:) = h;
           [xx,yy]=stairs(h,zi,'-k');
           xx = [0;xx;0];
           yy = [yy(1);yy;yy(end)];
           patch(xx,yy,[0.8,0.8,0.8]);
           [xx,yy]=stairs(h,zi,'-k');
           set(gca,'Fontsize',14,'YLim',[0 hmax2],'XLim',[0 0.1]);
           set(gca,'YDir','reverse','TickDir','out');
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
        set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
        set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
        if(isyn == 1)
           h1=plprof(mtru,hmax,NPL,'-w',1,0);
           set(h1,'LineWidth',2);
           h1=plprof(mtru,hmax,NPL,'--k',1,0);
           set(h1,'LineWidth',2);
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
        if(ibucking > 0);
          set(gca,'XTick',[0.25 0.5 0.75]);
          xlabel('Porosity');
        else; 
          set(gca,'XTick',[1450:50 :4000]);
          xlabel('Comp. Vel. (m/s)');
        end;
        box on;

        subplot(h2)
        set(gca,'FontSize',14);
        if(imarg == 1)
           pcolor(rlim,z,Nr);shading flat;
           if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
           if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
        end;
        set(h2,'layer','top')
        set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
        set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
        if(ibucking == 0);
          xlabel('Density (g/ccm)');
          set(gca,'XTick',[1.3:.2:3]);
        elseif(ibucking == 1);
          xlabel('Tortuosity');
        elseif(ibucking == 3);
          xlabel('log10(\gamma_P)');
          set(gca,'XTick',[1:2:26]);
        end;
        if(isyn == 1)
           h1=plprof(mtru,hmax,NPL,'-w',2,0);
           set(h1,'LineWidth',2);
           h1=plprof(mtru,hmax,NPL,'--k',2,0);
           set(h1,'LineWidth',2);
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
        set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
        if(ibucking == 0);
          xlabel('log_{10} \alpha_P (dB/m/kHz)');
          %set(gca,'XTick',[0:0.3:1]);
          set(gca,'XTick',[-10:1:10]);
        elseif(ibucking == 1);
          xlabel('log Transition freq.');
        elseif(ibucking == 2);
          xlabel('strain hardening idx');
        elseif(ibucking==3 & ivgsnosh==0);
          xlabel('log10(Mat. idx)');
          set(gca,'XTick',[-10:1:2]);
        elseif(ibucking==3 & ivgsnosh==1);
          xlabel('log10(\tau_P)');
        end;
        if(isyn == 1)
           h1=plprof(mtru,hmax,NPL,'-w',3,0);
           set(h1,'LineWidth',2);
           h1=plprof(mtru,hmax,NPL,'--k',3,0);
           set(h1,'LineWidth',2);
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
        cmap = colormap(jet);
        colormap(cmap);
        box on;

        if((ibucking==3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
          subplot(h5)
          set(gca,'FontSize',14);
          if(imarg == 1)
            pcolor(tlim,z,Nt);shading flat;
            if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
            if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
          end;
          set(h5,'layer','top')
          set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax2]);
          set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
          cmap = colormap(jet);
          colormap(cmap);
          box on;
          if(isyn == 1)
            h1=plprof(mtru,hmax,NPL,'-w',4,0);
            set(h1,'LineWidth',2);
            h1=plprof(mtru,hmax,NPL,'--k',4,0);
            set(h1,'LineWidth',2);
          end;
          if(ibucking == 0);
            xlabel('V_S (m/s)');
            set(gca,'XTickLabel',[0:750:2500]);
            set(gca,'XTick',[0:750:2500]);
          elseif(ibucking == 3);
            xlabel('log10(\tau_P)');
            set(gca,'XTick',[-15:1:0]);
          end;

          if(ishear==1);
            subplot(h6)
            set(gca,'FontSize',14);
            if(imarg == 1)
              pcolor(tlim2,z,Nt2);shading flat;
              if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
              if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
            end;
            set(h6,'layer','top')
            set(gca,'Fontsize',14,'XLim',[pmin(5) pmax(5)],'YLim',[0 hmax2]);
            set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
            cmap = colormap(jet);
            colormap(cmap);
            box on;
            if(ibucking == 0);
              xlabel('\alpha_S (dB/m/kHz)');
              set(gca,'XTickLabel',[0:1:5]);
              set(gca,'XTick',[0:1:5]);
            elseif(ibucking == 3);
              xlabel('log(\gamma_S)');
              set(gca,'XTickLabel',[0:5:25]);
              set(gca,'XTick',[0:5:25]);
            end;
            if(isyn == 1)
              h1=plprof(mtru,hmax,NPL,'-w',5,0);
              set(h1,'LineWidth',2);
              h1=plprof(mtru,hmax,NPL,'--k',5,0);
              set(h1,'LineWidth',2);
            end;
          end;
        end;
       if(icore == 1)
          barstep1 = 5;
          barstep1r = 5;
          barstep2 = 10;
          barstep2r = 1;
          barstep3 = 10;
          subplot(h1);
          plot(c1(:,2),c1(:,1),'w','Linewidth',1);
       %   errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
       %              10.*ones(size(c1(1:barstep1:end,2))),...
       %              zeros(size(c1(1:barstep1:end,2))),'w','w');
          plot(c2(:,2),c2(:,1),'*g','Linewidth',1);
       %   errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
       %             10.*ones(size(c2(1:barstep2:end,2))),...
       %              zeros(size(c2(1:barstep2:end,2))),'g','g');
          plot(c3(:,2),c3(:,1),'r','Linewidth',1);
       %   errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
       %              10.*ones(size(c3(1:barstep3:end,2))),...
       %              zeros(size(c3(1:barstep3:end,2))),'r','r'); 
          plot(c4(:,2),c4(:,1),'c','Linewidth',1);
          subplot(h2);
          plot(r1(:,2),r1(:,1),'w','Linewidth',1);
       %   errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
       %              2./100.*r1(1:barstep1r:end,2),...
       %              zeros(size(r1(1:barstep1r:end,2))),'w','w');
          plot(r2(:,2),r2(:,1),'c','Linewidth',1);
       %   errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
       %              2./100.*r2(1:barstep2r:end,2),...
       %              zeros(size(r2(1:barstep2r:end,2))),'w','w');
          plot(r3(:,2),r3(:,1),'r','Linewidth',1);
       %   errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
       %              2./100.*r3(1:barstep3:end,2),...
       %              zeros(size(r3(1:barstep3:end,2))),'r','r');
       end;
     end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%  COMPUTE DATA FIT
    %%
    if(idatfit == 1)
      disp(['load ensemble files',repfile]);
    %  reptmp=dlmread(repfile);
    %  NDAVE = length(reptmp)/NBAND;
    %  for j=1:NDAVE;
    %    for iband=1:length(bands);
    %      ref(:,iband,j,ifile2) = reptmp(NBAND*(j-1)+iband,:);
    %    end;
    %  end;
    %  ref = ref(:,:,1:2:end,:);
      load(repfile);
      ntmp1 = size(ref,1);
      ntmp2 = size(ref,2);
      ntmp3 = size(ref,3);
      ref2(1:ntmp1,1:ntmp2,1:ntmp3,ifile2) = ref;
      NDAVE = size(ref,3);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%
      %% PLOT DATA MISFIT
      %%
      nx = 3;
      ny = ceil(length(bands)/nx);
      if(length(bands)>12);ny = 3;end;
      xim = 0.01;
      yim = 0.05/ny;
      xymarg = [0.07 0.04 0.04 0.14];
      [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
      fig7 = figure('visible','off');hold on; box on;
      figw = 11.6;
      figh = 7.6;
      set(fig7,'PaperUnits','inches','PaperPosition',[0 0 figw figh])

    %  if(length(bands)>12);fig77 = figure('visible','off');hold on; box on;end;
      ii = 1;
      diffang = diff(angobs)/2;
      angplt = [angobs(1)-2*diffang(1),angobs-[diffang,diffang(end)],...
                angobs(end)+2*diffang(end)];
      for i=1:length(bands);
    %    figure(fig7);
        jj = i;
    %    if(i>12);
    %      figure(fig77);
    %      jj=i-12;
    %    end;
        subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
        set(gca,'FontSize',FNT);
        refmean(:,:,ifile2) = mean(ref2(:,:,:,ifile2),3);
    %    set(gca,'XLim',[angobs(1)-2 angobs(end)+2]);
        set(gca,'XLim',[angplt(1) angplt(end)]);
        set(gca,'FontSize',FNT);
        set(gca,'LineWidth',1);
        for j=1:length(angobs);
           y = squeeze(ref2(j,i,1:ntmp3,ifile2));
           [nf(j,i,:,ifile2)] = hpd(y,100,95);
        end;
        for j=1:NDAVE;
           plot(angobs,ref2(:,i,j,ifile2),':r');
        end;
        plot(angobs,nf(:,i,1,ifile2),'-w');
        plot(angobs,nf(:,i,1,ifile2),'--k');
        plot(angobs,nf(:,i,2,ifile2),'-w');
        plot(angobs,nf(:,i,2,ifile2),'--k');
        %plot(angobs,dobs(:,i,ifile2),'xk','Linewidth',1.5);
    %    idx1 = find(rex(:,i,ifile2)==1);
    %    idx2 = find(rex(:,i,ifile2)==0);
    %    plot(angobs(idx1),dobs(idx1,i,ifile2),'xk','Linewidth',1.5);
    %    plot(angobs(idx2),dobs(idx2,i,ifile2),'x', 'color', [0.7 0.7 0.7],'Linewidth',1.5);
        plot(angobs(:),squeeze(dobs(:,i,ifile2)),'xk','Linewidth',1.5);
        clear idx1 idx2;
        ylim  = [0,1.];
        ytick = [0,0.2,0.4,0.6,0.8];
        text(55,0.9,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    %    if(max(max(dobs))> 1.)
    %      ylim  = [0,1.3];
    %      ytick = [0,0.3,0.6,0.9,1.2];
    %      text(60,1.,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    %    elseif(max(max(dobs))> 0.64)
    %      ylim  = [0,1.01];
    %      ytick = [0,0.3,0.6,0.9];
    %      text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    %    else
    %      ylim  = [0,0.66];
    %      ytick = [0,0.2,0.4,0.6];
    %      text(60,0.55,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    %    end
        set(gca,'YLim',ylim,'XLim',[angplt(1) angplt(end)]);
        if ((i == (ii-1)*nx+1))
          set(gca,'YTickLabel',ytick,'YTick',ytick);
          ylabel('Refl. Coeff.');
          ii = ii + 1;
        else
          set(gca,'YTickLabel',[],'YTick',ytick);
        end
        if (i > length(bands)-nx)
          set(gca,'XTickLabel',[10:10:90],'XTick',[10:10:90]);
          xlabel('Angle (deg.)');
        else
          set(gca,'XTickLabel',[],'XTick',[10:10:90]);
        end
      end;
    end; % idatfit

    if(isave == 1)
       saveas(fig2,strcat(plotfile2,plotext2),'png');
       saveas(fig3,strcat(plotfile3,plotext2),'png');
       saveas(fig4,strcat(plotfile4,plotext2),'png');
       saveas(fig5,strcat(plotfile5,plotext2),'png');
       saveas(fig6,strcat(plotfile6,plotext2),'png');
       if(iar == 1)
    %     saveas(fig8,strcat(plotfile2,plotext1),'fig');
         saveas(fig8,strcat(plotfile2,plotext2),'png');
       end
       if(idatfit == 1);
        saveas(fig7,strcat(plotfile7,plotext2),'png');
      end;
    end;
    %%
    %% Plot track profile
    %%
    clear idxh idxc idxr idxa prof;
    if(isyn == 1)
      %% Find index for current model
      if(ktru > 0)
         idxh = (([1:ktru]-1)*NPL)+1;
         idxc = (([1:ktru]-1)*NPL)+2;
         idxr = (([1:ktru]-1)*NPL)+3;
         idxa = (([1:ktru]-1)*NPL)+4;
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
    close all;
    %close(fig2);
    %close(fig3);
    %close(fig33);
    %close(fig4);
    %close(fig5);
    %close(fig6);
    %if(idatfit == 1); close(fig7); end;
    %if(iar == 1); close(fig8); end;
  end; %% End loop over files
  save('track.mat','-v7.3');
  ifile2 = istart - 1;
else
  ifile2tmp = istart-1;
  load track.mat
  ifile2 = ifile2tmp;
  hmax2 = 4;
end;

if isecondfig == 1;
figw = 12.6;
figh = 7.0;
%figw = 12;
%figh = 5;
ZLIM = 6;

for ifile=1:NPING;
  ifile2 = ifile2 + 1;
  if(ifile2 > pinglist(end));exit;end;
  filename = files(ifile2).name;
  disp(filename);
%  filebase    = strrep(filename,'particles.mat','');
%  repfile     = strrep(filename,'_particles.mat','_replicas.mat');
%  datafile    = strrep(filename,'_particles.mat','.txt');
  filebase    = strrep(filename,'sample.mat','');
  repfile     = strrep(filename,'_sample.mat','_rep_ens.mat');
  datafile    = strrep(filename,'_sample.mat','.txt');
  plotext1    = 'fig';
  plotext2    = 'png';
  plotfile10  = strcat(filebase,'track_uncertainty.');
  plotfile11  = strcat(filebase,'track_uncertainty_data.');
 
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
  %hmax2  = hmax;

  diffang = diff(angobs)/2;
  angplt = [angobs(1)-2*diffang(1),angobs-[diffang,diffang(end)],...
            angobs(end)+2*diffang(end)];
  
  loc2 = [0.05, 0.79, 0.05, 0.145, 0.43, 0.715, 0.05, ...
           0.3675, 0.685;...
           0.82, 0.82, 0.29, 0.29, 0.29, 0.29, 0.065, 0.065, 0.065];
   
   spw1 = 0.09;
   spw2 = 0.28;
   spw3 = 0.65;
   spw4 = 0.20;
   spw5 = 0.31;
   sph1= 0.17;
   sph2= 0.45;
   sph3= 0.15;

   figw = 13;
   figh = 8.5;         
   fig11=figure('visible','off');
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

  %%
  %% MAP ensemble velocity track profile
  %%
  subplot(h1);hold on;box on;
  set(gca,'FontSize',FNT);
  imagesc(pinglist,z,c_mead);shading flat;
  plot([pinglist(ifile2), pinglist(ifile2)],[0, z(end)],'k','LineWidth',1.5);
  plot([pinglist(ifile2), pinglist(ifile2)],[0, z(end)],'--w','LineWidth',1.5);
  xlabel('Ping No.');
  ylabel('Depth (m)');
  set(gca,'Layer','top','YDir','reverse','TickDir','out');box on;
  set(gca,'ticklength',0.5*get(gca,'ticklength'))
  set(gca,'Layer','top','YDir','reverse','XTick',[pinglist(1):dping:pinglist(end)],'XTickLabel',[pinglist(1):dping:pinglist(end)]);
  set(gca,'XLim',[pinglist(1) pinglist(end)],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
  c1=colorbar;ylabel(c1,'Velocity (m/s)')

  subplot(h2);hold on;box on;
  set(gca,'FontSize',FNT);
  [xx,yy]=stairs(limgl,kgl(ifile2,:),'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(limgl,kgl(ifile2,:),'k');
  if(isyn == 1)
    plot([track_env(ifile2,1), track_env(ifile2,1)],[0, 1],'k','LineWidth',1.5);
    plot([track_env(ifile2,1), track_env(ifile2,1)],[0, 1],'--w','LineWidth',1.5);
  end;
  xlabel('No. interfaces');
  %ylabel('Probability');
  set(gca,'YTick',[0:.2:1],'TickDir','out');
  set(gca,'XLim',[0.5 15+0.5],'YLim',[0. 0.7]);

  
  subplot(h3);hold on;box on;
  set(gca,'FontSize',FNT);
  if(imarg == 1)
     [xx,yy]=stairs(hgl(ifile2,:),zi,'-k');
     xx = [0;xx;0];
     yy = [yy(1);yy;yy(end)];
     patch(xx,yy,[0.8,0.8,0.8]);
     [xx,yy]=stairs(hgl(ifile2,:),zi,'-k');
     if(isyn == 1)
       for jk=1:ktru;
         %mtrutmp(ifile2,(jk-1)*4+1)
         %ifile2
         plot([0,1],[mtrutmp(ifile2,(jk-1)*4+1), mtrutmp(ifile2,(jk-1)*4+1)],'k','LineWidth',1);
         plot([0,1],[mtrutmp(ifile2,(jk-1)*4+1), mtrutmp(ifile2,(jk-1)*4+1)],'--w','LineWidth',1);
       end;
     end;
     set(gca,'Fontsize',12,'YLim',[0 hmax2],'TickDir','out');
     set(gca,'YDir','reverse','XLim',[0 0.1]);
     set(gca,'YTickLabel',[0:100],'YTick',[0:100]);
     xlabel('Interface prob.');
     ylabel('Depth (m)');
     box on;
  end;

  subplot(h4)
  set(gca,'FontSize',FNT);
  if(imarg == 1)
     imagesc(clim,z,Ncgl(:,:,ifile2));shading flat;
     if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
     if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
  end;
  %surf(clim,z,Nc);shading flat;
  set(h4,'layer','top')
  set(gca,'Fontsize',12,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
  set(gca,'YDir','reverse','TickDir','out');
  xlabel('Velocity (m/s)');
  if(isyn == 1)
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',1,0);
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',1,0);
  end;
%  if(imap == 1)
%     plprof(mmap,hmax,'--k',1,0);
%  end;
  if(imax == 1)
    for i=1:length(z);
       [cmx,j] = max(Nc(i,:));
       c_max(i) = clim(j);
    end;
    plot(c_max,z,'w');
  end;
  set(gca,'YTickLabel',[],'YTick',[0:1:100]);
  set(gca,'XTick',[1500:100:3000]);
  box on;

  subplot(h5)
  set(gca,'FontSize',FNT);
  if(imarg == 1)
     imagesc(rlim,z,Nrgl(:,:,ifile2));shading flat;
     if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
     if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
  end;
  %surf(rlim,z,Nr);shading flat;
  set(h5,'layer','top','TickDir','out')
  set(gca,'Fontsize',12,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
  set(gca,'YDir','reverse');
  xlabel('Density (g/cm^3)');
  if(isyn == 1)
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'-k',2,0);
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax,NPL,'--w',2,0);
  end
%  if(imap == 1)
%     plprof(mmap,hmax,'--k',2,0);
%  end
  if(imax == 1)
    for i=1:length(z);
       [rmx,j] = max(Nr(i,:));
       r_max(i) = rlim(j);
    end;
    plot(r_max,z,'w');
  end;
  set(gca,'YTickLabel',[],'YTick',[0:1:20]);
  set(gca,'XTick',[1.3:.15:3.0]);
  box on;

  subplot(h6)
  set(gca,'FontSize',FNT);
  if(imarg == 1)
     imagesc(alim,z,Nagl(:,:,ifile2));shading flat;
     if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
     if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
  end;
  %surf(alim,z,Na);shading flat;
  set(h6,'layer','top')
  set(gca,'Fontsize',12,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
  set(gca,'YDir','reverse','TickDir','out');
  %xlabel('log10(Attenuation (dB/m/kHz))');
  xlabel('Attenuation log_{10} (dB/m/kHz)');
  if(isyn == 1)
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'-k',3,0);
     plprof(mtrutmp(ifile2,1:NFPTRU),hmax2,NPL,'--w',3,0);
  end;
%  if(imap == 1)
%     plprof(mmap,hmax2,'--k',3,0);
%  end;
  if(imax == 1)
    for i=1:length(z);
       [amx,j] = max(Na(i,:));
       a_max(i) = alim(j);
    end;
    plot(a_max,z,'w');
  end;
  set(gca,'YTickLabel',[],'YTick',[0:1:100]);
  set(gca,'XTick',[-2.5:1:.5]);
  %set(gca,'XTick',[0:.2:2]);
  %cmap = colormap(flipud(gray));
  cmap = colormap(jet);
  %cmap(1,:) = [1 1 1];
  colormap(cmap);
  %colorbar('peer',h1,'location','WestOutside');
  box on;
  
  ii = 1;
  for i=1:1:length(bands);
   if(i == 1 | i == 3 | i == 6);
     if(i==1);subplot(h7);hold on;box on;end;
     if(i==3);subplot(h8);hold on;box on;end;
     if(i==6);subplot(h9);hold on;box on;end;
     set(gca,'FontSize',FNT);
     
%    set(gca,'XLim',[angobs(1)-2 angobs(end)+2]);
    set(gca,'XLim',[angplt(1) angplt(end)]);
    set(gca,'FontSize',FNT);
    set(gca,'LineWidth',1);
    for j=1:NDAVE;
       plot(angobs,ref2(:,i,j,ifile2),':r');
    end;
    if(idatfit == 1);
        plot(angobs,nf(:,i,1,ifile2),'-w');
        plot(angobs,nf(:,i,1,ifile2),'--k');
        plot(angobs,nf(:,i,2,ifile2),'-w');
        plot(angobs,nf(:,i,2,ifile2),'--k');
    %    idx1 = find(rex(:,i,ifile2)==1);
    %    idx2 = find(rex(:,i,ifile2)==0);
    %    plot(angobs(idx1),dobs(idx1,i,ifile2),'xk','Linewidth',1.5);
    %    plot(angobs(idx2),dobs(idx2,i,ifile2),'x', 'color', [0.7 0.7 0.7],'Linewidth',1.5);
    end;
    plot(angobs(:),squeeze(dobs(:,i,ifile2)),'xk','Linewidth',1.5);
    clear idx1 idx2;
    if(max(max(dobs))> 1.)
      ylim  = [0,1.3];
      ytick = [0,0.3,0.6,0.9,1.2];
      text(60,1.,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    elseif(max(max(dobs))> 0.64)
      ylim  = [0,1.01];
      ytick = [0,0.3,0.6,0.9];
      text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    else
      ylim  = [0,0.66];
      ytick = [0,0.2,0.4,0.6];
      text(60,0.55,[num2str(bands(i)) ' Hz'],'FontSize',FNT,'Color',[0,0,0])
    end
    set(gca,'YLim',ylim,'XLim',[angplt(1) angplt(end)]);
    %if ((i == (ii-1)*nx+1))
    if (i == 1)
      set(gca,'YTickLabel',ytick,'YTick',ytick);
      ylabel('Refl. coeff.','FontSize',FNT);
      ii = ii + 1;
    else
      set(gca,'YTickLabel',[],'YTick',ytick);
    end
%    if (i > length(bands)-nx)
      set(gca,'XTickLabel',[10:10:90],'XTick',[10:10:90]);
      xlabel('Angle (deg.)');
%    else
%      set(gca,'XTickLabel',[],'XTick',[10:10:90]);
%    end
   end;
  end;
%  saveas(fig11,strcat(plotfile11,plotext2),'png');
  print(fig11,'-painters','-r250',strcat(plotfile11,plotext2),'-dpng');
%  print(fig11,'-painters','-r250',strcat(plotfile11,plotext3),'-depsc');
  close(fig11);
end;
end; %% End loop over files

return;
