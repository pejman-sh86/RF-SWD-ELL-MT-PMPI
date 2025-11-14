%----------------------------------------------------------
% Plot nested sampling results
%----------------------------------------------------------
function [] = plot_nest_disp();

i_save  = 1;

i_power = 0;   % Layer with power law?
               % 0 = uniform layer
               % 1 = power law layer
               % 2 = linear grad layer
isd     = 1;   % Standard deviations were inverted for
nsd     = 1;   % Standard deviations were inverted for

ilpl    = 1;   % Layer that has power law/linear gradient

ippd1 = 0;     % Nested sampling
ippd2 = 1;     % MH sampling
ippd3 = 1;     % AIS sampling

Tstar2 = 1;
i_scaley = 1;
itrue = 0;

if(i_power == 0)
   nx = 4;
   ny = 4;
elseif(i_power == 1)
   nx = 6;
   ny = 4;
elseif(i_power == 2)
   nx = 5;
   ny = 4;
end

% INPUT:
ppd1    = 'nest_3lay_sample.mat';
ppd2    = '../mh_3lay_sample.mat';
ppd3    = 'ais_3lay_sample.mat';

% OUTPUT:
if(i_save == 1)
   mapfile       = 'mh_3lay_nestmap.txt';
   meanfile      = 'mh_3lay_mean.txt';
   nestfile      = 'mh_3lay_nest.eps';
   nestfilefig   = 'mh_3lay_nest.fig';
   aisfile       = 'mh_3lay_ais.eps';
   aisfilefig    = 'mh_3lay_ais.fig';
   plotfile      = 'mh_3lay_marg.eps';
   plotfilefig   = 'mh_3lay_marg.fig';
   sd_plotfile   = 'mh_3lay_sd_marg.eps';
   sd_plotfilefig= 'mh_3lay_sd_marg.fig';
end

nsubfig = nx*ny;
xim = 0.03;
yim = 0.3/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 160;    % 
nbin2 = 160;    % 
nbin3 = 70;    % For AIS convergence plot

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',12,'Height',5*ny,'Color','rgb',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
opts2= struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',4,'Height',3,'Color','rgb',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
opts3= struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',4,'Color','rgb',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
if(i_power == 0)
   plo = [0.,     0., 1400., 1.3; ...
          0.,     0., 1400., 1.3; ...  
          0.,     0., 1400., 1.3; ...  
          0.,     0., 1400., 1.3; ...  
          0.,     0., 1400., 1.3; ...  
          0.,     0., 1400., 1.3];
   phi = [  5.,  20., 2200., 2.2; ...
            5., 200., 2200., 2.2; ...
            5., 200., 2200., 2.2; ...
            5., 200., 2200., 2.2; ...
           10., 200., 2200., 2.2; ...
           10., 200., 2200., 2.2];
elseif(i_power == 1)
   plo = [0.,     0.,   8.5,   0., 1400., 1.3; ...
          0.,     0.,   0.,   0., 1400., 1.3; ...
          0.,     0.,   0.,   0., 1400., 1.3; ...
          0.,     0.,   0.,   0., 1400., 1.3; ...
          0.,     0.,   0.,   0., 1400., 1.3; ...
          0.,     0.,   0.,   0., 1400., 1.3];
   phi = [ 60.,   5.,  10.5, 200., 2200., 2.2; ...
           60.,  20.,  20., 400., 2200., 2.2; ...
           80., 300., 500., 400., 2200., 2.2; ...
           80., 300., 500., 400., 2200., 2.2; ...
           80., 300., 500., 400., 2200., 2.2; ...
           80., 300., 500., 400., 2200., 2.2];
elseif(i_power == 2)
   plo = [0.,     0.,   0.,  1450., 1.3; ...
          0.,     0.,   0.,  1450., 1.3; ...
          0.,     0.,   0.,  1450., 1.3; ...
          0.,     0.,   0.,  1450., 1.3; ...
          0.,     0.,   0.,  1450., 1.3; ...
          0.,     0.,   0.,  1450., 1.3];
   phi = [ 20.,   5., 100.,  2200., 2.2; ...
           40.,  10., 200.,  2200., 2.2; ...
           80., 500.,1200.,  2200., 2.2; ...
           80., 500.,1200.,  2200., 2.2; ...
           80., 500.,1200.,  2200., 2.2; ...
           80., 500.,1200.,  2200., 2.2];
end
if(isd == 1);
   plosd = [0. 0.];
   phisd = [3. 3.];
end;

if(i_power == 0)
   F(1).xticks = [0.,2.,4.;...
                  0.,2.,4. ];
   F(2).xticks = [0, 10, 20,   30,  40;...
                  0, 50, 100, 150, 200];
   F(3).xticks = [1500, 1800, 2100;...
                  1500, 1800, 2100];
   F(4).xticks = [1.3, 1.5, 1.7, 1.9, 2.1;...
                  1.3, 1.5, 1.7, 1.9, 2.1];
elseif(i_power == 1)
   F(1).xticks = [0,20,40,60];
   F(2).xticks = [0,2,4];
   F(3).xticks = [9,10];
   F(4).xticks = [0, 50, 100, 150, 200;...
                  0,100, 200, 300, 400];
   F(5).xticks = [1500, 1700, 1900, 2100;...
                  1500, 1700, 1900, 2100];
   F(6).xticks = [1.3, 1.5, 1.7, 1.9, 2.1;...
                  1.3, 1.5, 1.7, 1.9, 2.1];
elseif(i_power == 2)
   F(1).xticks = [0,10,20];
   F(2).xticks = [0,2,4];
   F(3).xticks = [0, 20,  40,  60,  80;...
                  0, 50, 100, 150, 200];
   F(4).xticks = [1500, 1700, 1900, 2100;...
                  1500, 1700, 1900, 2100];
   F(5).xticks = [1.3, 1.5, 1.7, 1.9, 2.1;...
                  1.3, 1.5, 1.7, 1.9, 2.1];
end

if(ippd1 == 1)
   load(ppd1);
   if(isd == 0)
      m1       = A(:,6:end);
   else
      m1       = A(:,6:end-nsd);
      sd1      = A(:,end-nsd+1:end);
   end
   logL     = A(:,1);
   logwidth = A(:,2);
   logwt    = A(:,3);
   H1        = A(:,4);
   logZ1     = A(:,5);
   w = exp(logwt-logZ1(end));
   fprintf(1,'Nested Sampling Evidence estimate: Znst = %12.4f \n',logZ1(end));
end;
if(ippd2 == 1)
   load(ppd2);
   if(isd == 0)
      m2  = A(:,2:end);
   else
      m2  = A(:,2:end-nsd);
      sd2 = A(:,end-nsd+1:end);
   end
   E2  = A(:,1);
end
if(ippd3 == 1)
   load(ppd3);
   if(isd == 0)
      m3     = A(:,3:end);
   else
      m3     = A(:,3:end-nsd);
      sd3    = A(:,end-nsd+1:end);
   end
   logwt3 = A(:,1);
   logL3  = A(:,2);

   logwtsum3(1) = logwt3(1);
   logZ3(1) = logwtsum3(1)-log(1);
   for i = 2:length(logwt3);
      logwtsum3(i) = log(exp(logwtsum3(i-1))+exp(logwt3(i)));
      logZ3(i) = logwtsum3(i)-log(i);
   end;
   %% Estimate Z uncertainty:
   NZ_sm = 20;
   Nsmall = floor(length(logwt3)/NZ_sm);
   for j = 1:NZ_sm;
      logwtsum_sm(1) = logwt3(((j-1)*Nsmall)+1);
      logZ_sm(j,1) = logwtsum_sm(1)-log(1);
      for i = 2:Nsmall;
         logwtsum_sm(i) = log(exp(logwtsum_sm(i-1))+exp(logwt3((j-1)*Nsmall+i)));
         logZ_sm(j,i) = logwtsum_sm(i)-log(i);
      end;
   end;
   mean(logZ_sm(:,end))
   fprintf(1,'AIS Evidence estimate: Zais = %12.4f +/- %12.4f\n',logZ3(end),std(logZ_sm(:,end)));
end

if(ippd1 == 1)
   % Plot nested sampling convergence plot
   gc1=figure(1);
   plot(exp(w));
   xlabel('No. iteration');
   ylabel('Importance weight');
   set(gca,'XLim',[1 length(w)])
end;
if(ippd3 == 1)
   % Plot AIS convergence plot
   gc3=figure(3);
   subplot(1,2,1);hold on;box on;
%   [nZ,limZ] = hist(logwt3,30);
   [nZ,limZ] = hist_ais(logwt3,ones(size(logwt3)),80);
   limZ_2 = [min(logwt3);limZ'; max(logwt3)];
   nZ_2 = [0; nZ'; 0];
   fill(limZ_2,nZ_2,[0.8,0.8,0.8]);
   set(gca,'YTick',[],'FontSize',14);
   xlabel('log Weight');
   subplot(1,2,2);hold on;box on;
   plot(logZ3);
   set(gca,'XLim',[0 11000],'XScale','log')
   set(gca,'XTick',[10 100 1000 10000])
   xlabel('No. Trajectory');
   ylabel('log(Z)');

   % Compute weight variance
   var_wt = var(exp(logwt3));
   var_logwt = var(logwt3);
   wt_star = exp(logwt3)/length(logwt3)*sum(exp(logwt3));
   var_wt_star = var(wt_star);
   fprintf(1,'AIS weight variance: var(wt) = %12.4f \n',var_wt);
   fprintf(1,'AIS logweight variance: var(logwt) = %12.4f \n',log(var_logwt));
   fprintf(1,'AIS weight star variance: var(wt*) = %12.4f \n',log(var_wt_star));
   fprintf(1,'AIS effective sample size: N/(1+var(wt*)) = %12.4f \n',length(logwt3)/(1+var_wt_star));
   % Plot Neal style plot:
   gc5=figure(5);hold on;box on;
   for ipar=1:8
      subplot(3,3,ipar);hold on;box on;
      plot(m3(:,ipar),logwt3,'.k')
   end;
end

%xparname = [{'a1'},{'a2'},{'a3'},{'a4'},...
%            {'a5'},{'a6'},{'a7'},{'a8'},...
%            {'a9'},{'a10'},{'a11'},{'a12'}];
if(i_power == 0)
  xparname = [{'h1'},{'vs1'},{'vp1'},{'r1'};...
              {'h2'},{'vs2'},{'vp2'},{'r2'};...
              {'h3'},{'vs3'},{'vp3'},{'r3'};...
              {'h4'},{'vs4'},{'vp4'},{'r4'};...
              {'h5'},{'vs5'},{'vp5'},{'r5'};...
              {'h6'},{'vs6'},{'vp6'},{'r6'};...
              {'h7'},{'vs7'},{'vp7'},{'r7'};...
              {'h8'},{'vs8'},{'vp8'},{'r8'};...
              {'h9'},{'vs9'},{'vp9'},{'r9'};...
              {'h10'},{'vs10'},{'vp10'},{'r10'}];
elseif(i_power == 1)
  xparname = [{'h1'},{'vs1'},{'vs1'},{'vs1'},{'vp1'},{'r1'};...
              {'h2'},{'vs2'},{'vs2'},{'vs2'},{'vp2'},{'r2'};...
              {'h3'},{'vs3'},{'vs3'},{'vs3'},{'vp3'},{'r3'};...
              {'h4'},{'vs4'},{'vs4'},{'vs4'},{'vp4'},{'r4'};...
              {'h5'},{'vs5'},{'vs5'},{'vs5'},{'vp5'},{'r5'};...
              {'h6'},{'vs6'},{'vs6'},{'vs6'},{'vp6'},{'r6'};...
              {'h7'},{'vs7'},{'vs7'},{'vs7'},{'vp7'},{'r7'};...
              {'h8'},{'vs8'},{'vs8'},{'vs8'},{'vp8'},{'r8'};...
              {'h9'},{'vs9'},{'vs9'},{'vs9'},{'vp9'},{'r9'};...
              {'h10'},{'vs10'},{'vs10'},{'vs10'},{'vp10'},{'r10'}];
elseif(i_power == 2)
  xparname = [{'h1'},{'vs1'},{'vs1'},{'vp1'},{'r1'};...
              {'h2'},{'vs2'},{'vs2'},{'vp2'},{'r2'};...
              {'h3'},{'vs3'},{'vs3'},{'vp3'},{'r3'};...
              {'h4'},{'vs4'},{'vs4'},{'vp4'},{'r4'};...
              {'h5'},{'vs5'},{'vs5'},{'vp5'},{'r5'};...
              {'h6'},{'vs6'},{'vs6'},{'vp6'},{'r6'};...
              {'h7'},{'vs7'},{'vs7'},{'vp7'},{'r7'};...
              {'h8'},{'vs8'},{'vs8'},{'vp8'},{'r8'};...
              {'h9'},{'vs9'},{'vs9'},{'vp9'},{'r9'};...
              {'h10'},{'vs10'},{'vs10'},{'vp10'},{'r10'}];
end

if(ippd1 == 1)
   [logL_max i_max] = max(logL);
   xmap = m1(i_max,:);
   if(i_save == 1)
       save(mapfile,'xmap','-ASCII');
   end;
end;

if(ippd1 == 1)
   npar = size(m1,2);
elseif(ippd2 == 1)
   npar = size(m2,2);
elseif(ippd3 == 1)
   npar = size(m3,2);
end

if(i_power == 0)
   nlay = (npar-3)/4;
elseif(i_power == 1)
   nlay = (npar-5)/4;
elseif(i_power == 2)
   nlay = (npar-4)/4;
end

npar
nlay

%pplo = [plo(1:nlay*4), plo(end-2:end)];
%pphi = [phi(1:nlay*4), phi(end-2:end)];

if(ippd1 == 1)
   nmods = size(m1,1);
elseif(ippd2 == 1)
   nmods = size(m2,1);
elseif(ippd3 == 1)
   nmods = size(m3,1);
end
nfig = ceil(npar/nsubfig);

%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

%
% Make subplot mask:
%
nlay+1
nx
idxsub = zeros(nlay+1,nx);
idxsub(1:end-1,[1,end-2:end]) = 1;
if(i_power > 0)
   idxsub(ilpl,:) = 1;
end
idxsub(end,end-2:end) = 1;

ipar = 1;
for iy = 1:nlay+1
for ix = 1:nx

   if(idxsub(iy,ix) == 1)
   if(ippd1 ==1)
      clear idx;
      idx  = find(m1(:,ipar)>plo(iy,ix) & m1(:,ipar)<phi(iy,ix));
      [n1(:,ipar),lim1(:,ipar)] = hist_wt(m1(idx,ipar),w(idx),nbin1);
      clear idx;
   end
   if(ippd2 ==1)
      idx = find(m2(:,ipar)>plo(iy,ix) & m2(:,ipar)<phi(iy,ix));
      [n2(:,ipar),lim2(:,ipar)] = hist_tstar(m2(idx,ipar),E2(idx),nbin2,Tstar2);
      clear idx;
   end
   if(ippd3 ==1)
      idx = find(m3(:,ipar)>plo(iy,ix) & m3(:,ipar)<phi(iy,ix));
      [n3(:,ipar),lim3(:,ipar)] = hist_ais(m3(idx,ipar),logwt3(idx),nbin3);
      clear idx;
   end

   ipar = ipar +1;
   end
end
end
for ipar = 1:nsd
   if(isd == 1)
      if(ippd1 ==1)
         [n1sd(:,ipar),lim1sd(:,ipar)] = hist_wt(sd1(:,ipar),w,nbin1);
      end
      if(ippd2 ==1)
         [n2sd(:,ipar),lim2sd(:,ipar)] = hist_tstar(sd2(:,ipar),E2,nbin2,Tstar2);
      end
      if(ippd3 ==1)
         [n3sd(:,ipar),lim3sd(:,ipar)] = hist_ais(sd3(:,ipar),logwt3,nbin2);
      end
   end;
end;

%----------------------------------------------------------
%  Plot histograms:
%----------------------------------------------------------
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc2=figure(2);
isubfig = 1;
ipar = 1;
for iy = 1:nlay+1
for ix = 1:nx
   if(idxsub(iy,ix) == 1)
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      hold on;box on;

      if(ippd1 ==1)
         lim1_2 = [plo(iy,ix);lim1(1,ipar);
                  lim1(:,ipar);lim1(end,ipar); phi(iy,ix)];
         n1_2 = [0; 0; n1(:,ipar); 0; 0];
         n1_2 = n1_2*(phi(iy,ix)-plo(iy,ix));
         fill(lim1_2,n1_2,[0.8,0.8,0.8]);
         ylims_tmp1(iy,ix) = max(n1_2);
      end

      if(ippd2 ==1)
         lim2_2 = [plo(iy,ix);lim2(1,ipar);
                  lim2(:,ipar);lim2(end,ipar); phi(iy,ix)];
         n2_2 = [0; 0; n2(:,ipar); 0; 0];
         n2_2 = n2_2*(phi(iy,ix)-plo(iy,ix));
%         stairs(lim2_2,n2_2,'b');
         fill(lim2_2,n2_2,[0.8,0.8,0.8]);
         ylims_tmp2(iy,ix) = max(n2_2);
      end
      if(ippd3 ==1)
         lim3_2 = [plo(iy,ix);lim3(1,ipar);
                  lim3(:,ipar); lim3(end,ipar); phi(iy,ix)];
         n3_2 = [0;0; n3(:,ipar);0; 0];
         n3_2 = n3_2*(phi(iy,ix)-plo(iy,ix));
%         stairs(lim3_2,n3_2,'r');
         plot(lim3_2,n3_2,'r','LineWidth',1);
         ylims_tmp3(iy,ix) = max(n3_2);
      end

      set(gca,'layer','top')
      ipar = ipar +1;

   end
   isubfig = isubfig + 1;
%   end
end
end

if(ippd1 == 1)
   ylims1 = max(max(ylims_tmp1));
else
   ylims1 = 0;
end
if(ippd2 == 1)
   ylims2 = max(max(ylims_tmp2));
else
   ylims2 = 0;
end
if(ippd3 == 1)
   ylims3 = max(max(ylims_tmp3));
else
   ylims3 = 0;
end
ylims = max([ylims1,ylims2,ylims3]);

ipar = 1;
isubfig = 1;
for iy = 1:nlay+1
for ix = 1:nx
   if(idxsub(iy,ix) == 1)
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      set(gca,'YTick',[],'FontSize',14);
      set(gca,'XTick',[F(ix).xticks(iy,:)],'FontSize',14);
      if(iy == nlay+1 | (iy == nlay & ix < nx-2) | (idxsub(iy+1,ix) == 0)) 
         xlabel(xparname(iy,ix),'FontSize',14);
      else
%         set(gca,'XTickLabel',[]);
      end 
      set(gca,'XLim',[plo(iy,ix) phi(iy,ix)],'FontSize',14);

      if(i_scaley == 1)
         set(gca,'YLim',[0 ylims+.05*ylims]);
      end
      ipar = ipar +1;
   end
   isubfig = isubfig + 1;
end
end

if(isd == 1)
   gc4=figure(4);
   for ipar = 1:nsd;
      subplot(1,nsd,ipar);hold on;box on;
      set(gca,'XLim',[plosd(ipar) phisd(ipar)],'FontSize',14);
      if(ippd1 ==1)
         lim1sd_2 = [min(sd1(:,ipar));min(sd1(:,ipar));
                  lim1sd(:,ipar); max(sd1(:,ipar)); max(sd1(:,ipar))];
         n1sd_2 = [0; n1sd(1,ipar); n1sd(:,ipar); n1sd(end,ipar); 0];
         fill(lim1sd_2,n1sd_2,[0.8,0.8,0.8]);
      end

      if(ippd2 ==1)
         lim2sd_2 = [min(sd2(:,ipar));min(sd2(:,ipar));
                  lim2sd(:,ipar); max(sd2(:,ipar)); max(sd2(:,ipar))];
         n2sd_2 = [0; n2sd(1,ipar); n2sd(:,ipar); n2sd(end,ipar); 0];
%         stairs(lim2sd_2,n2sd_2,'b');
         fill(lim2sd_2,n2sd_2,[0.8,0.8,0.8]);
      end
      if(ippd3 ==1)
         lim3sd_2 = [min(sd3(:,ipar));min(sd3(:,ipar));
                  lim3sd(:,ipar); max(sd3(:,ipar)); max(sd3(:,ipar))];
         n3sd_2 = [0; n3sd(1,ipar); n3sd(:,ipar); n3sd(end,ipar); 0];
%         stairs(lim3sd_2,n3sd_2,'r');
         plot(lim3sd_2,n3sd_2,'r','LineWidth',1);
      end
      set(gca,'YTick',[],'FontSize',14);
      xlabel('Standard deviation (m/s)');

      set(gca,'layer','top')
   end;
end;

if(i_save == 1)
   if(ippd1 == 1)
      exportfig(gc1,nestfile,opts);
      hgsave(gc1,nestfilefig);
   end
   if(ippd3 == 1)
      exportfig(gc3,aisfile,opts3);
      hgsave(gc3,aisfilefig)
   end
   exportfig(gc2,plotfile,opts);
   hgsave(gc2,plotfilefig)
   if(isd == 1)
      exportfig(gc4,sd_plotfile,opts2);
      hgsave(gc4,sd_plotfilefig)
   end
end

return;
