%----------------------------------------------------------
% Compute Prior from PPD
%----------------------------------------------------------
function [] = plot_histrg();

i_save = 1;

Tstar = 1;
Tstarc = 1;
nx = 5;
ny = 8;

% INPUT:
ppd     = 'x_s21_1_5_25_4layrg_sample.mat';
data    = 'x_s21_1_5_25_4layrg.mat';
%ppd     = 'x_s16_1_10_40_4layrg_sample.mat';
%data    = 'x_s16_1_10_40_4layrg.mat';
%ppd     = 'x_s13_2_10_50_5layrgb_sample.mat';
%data    = 'x_s13_2_10_50_5layrgb.mat';
%ppd     = 'x_s01_1_3_20_2layrg_sample.mat';
%data    = 'x_s01_1_3_20_2layrg.mat';
%ppd     = 'x_s02_1_3_16_2layrg_sample.mat';
%data    = 'x_s02_1_3_16_2layrg.mat';
%ppd     = 'sim_C_1_7_40_6layrg_sample.mat';
%data    = 'sim_C_1_7_40_6layrg.mat';

% OUTPUT:
if(i_save == 1)
%   mapfile  = 'x_s01_1_3_20_2layrg_fgsmap.txt';
%   meanfile = 'x_s01_1_3_20_2layrg_mean.txt';
%   hpdfile  = 'x_s01_1_3_20_2layrg_hpds.txt';
%   plotfile = 'x_s01_1_3_20_2layrg_marg_a.eps';
%   mapfile  = 'x_s02_1_3_16_2layrg_fgsmap.txt';
%   meanfile = 'x_s02_1_3_16_2layrg_mean.txt';
%   hpdfile  = 'x_s02_1_3_16_2layrg_hpds.txt';
%   plotfile = 'x_s02_1_3_16_2layrg_marg_a.eps';
%   mapfile  = 'x_s13_2_10_50_5layrgb_fgsmap.txt';
%   meanfile = 'x_s13_2_10_50_5layrgb_mean.txt';
%   hpdfile  = 'x_s13_2_10_50_5layrgb_hpds.txt';
%   plotfile = 'x_s13_2_10_50_5layrgb_marg_a.eps';
%   mapfile  = 'x_s16_1_10_40_4layrg_fgsmap.txt';
%   meanfile = 'x_s16_1_10_40_4layrg_mean.txt';
%   hpdfile  = 'x_s16_1_10_40_4layrg_hpds.txt';
%   plotfile = 'x_s16_1_10_40_4layrg_marg_a.eps';
   mapfile  = 'x_s21_1_5_25_4layrg_fgsmap.txt';
   meanfile = 'x_s21_1_5_25_4layrg_mean.txt';
   hpdfile  = 'x_s21_1_5_25_4layrg_hpds.txt';
   plotfile = 'x_s21_1_5_25_4layrg_marg_a.eps';
%   mapfile  = 'sim_C_1_7_40_6layrg_fgsmap.txt';
%   meanfile = 'sim_C_1_7_40_6layrg_mean.txt';
%   hpdfile  = 'sim_C_1_7_40_6layrg_hpds.txt';
%   plotfile = 'sim_C_1_7_40_6layrg_marg_a.eps';
end

nsubfig = nx*ny;
xim = 0.03;
yim = 0.1/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 30;
nbin2 = 10;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',12,'Height',1*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
mtrue = [0.30, 1480, 1.30, 1.6, 0.01, ...
         0.35, 1510, 1.65, 0.30, ...
         1.20, 1500, 1.55, 0.30, ...
         0.40, 1700, 1.90, 0.30, ...
               1560, 1.70, 0.01 ];

load(data);
load(ppd);
m = A(:,2:end);
%m = A(:,2:end-1);
A2 = A;
E = A2(:,1);

lo = F(1).minlim(1:end-3);
hi = F(1).maxlim(1:end-3);

lob = F(1).minlim(end-2:end);
hib = F(1).maxlim(end-2:end);

F(1).xticks = [0.5, 1, 1.5];
F(2).xticks = [1500, 1600, 1700];
F(3).xticks = [1.4, 1.6, 1.8, 2.0];
F(4).xticks = [1.4, 1.6, 1.8, 2.0];
F(5).xticks = [0, 0.5, 1.0 ];
xparname = [{'h1'},{'c1'},{'r1t'},{'r1b'},{'a1'},...
            {'h2'},{'c2'},{'r2t'},{'r2'},{'a2'},...
            {'h3'},{'c3'},{'r3t'},{'r3'},{'a3'},...
            {'h4'},{'c4'},{'r4t'},{'r4'},{'a4'},...
            {'h5'},{'c5'},{'r5t'},{'r5'},{'a5'},...
            {'h6'},{'c6'},{'r6t'},{'r6'},{'a6'},...
            {'h7'},{'c7'},{'r7t'},{'r7'},{'a7'},...
            {'h8'},{'c8'},{'r8t'},{'r8'},{'a8'}];

[cost_min i_min] = min(E);
xmap = m(i_min,:);
if(i_save == 1)
    save(mapfile,'xmap','-ASCII');
end

for ipar = 1:size(m,2)

    xmean(ipar) = mean(m(:,ipar));

end;
if(i_save == 1)
     save(meanfile,'xmean','-ASCII');
end;

[EE idxE] = min(E);
mstart = m(idxE,:);

npar = size(m,2);
nlay = (npar-4)/4;
%plo = [0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0,...
%       0.0, 1450, 1.2, 1.2, 0.0];
%phi = [2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0,...
%       2.0, 1750, 2.1, 2.1, 1.0];
%
% Site 01
%
plo = [0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 1.2, 0.0];
phi = [2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0,...
       2.0, 1850, 2.2, 2.2, 1.0];


nmods = size(m,1);

%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

for ipar = 1:npar

   [n1(:,ipar),lim(:,ipar)] = hist_tstar(m(:,ipar),E,nbin1,Tstar);
%   [n1(:,ipar),lim(:,ipar)] = hist(m(:,ipar),nbin1);
   nf_int(ipar,:) = hpd(m(:,ipar),100,95);

end

%----------------------------------------------------------
%  Plot histograms:
%----------------------------------------------------------

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc1=figure(1);
isubfig = 1;
ipar = 1;
ipanel = 1;
for iy = 1:ny
for ix = 1:nx

   if(ipar <= npar)
      if(ix == 3 && iy >= 2) 
         isubfig = isubfig + 1;
      elseif(ipanel == nx*(nlay+1)-4) 
         isubfig = isubfig + 1;
      else
         if(ipar > npar-2) 
            ylims_tmp(iy,ix+1) = max(n1(:,ipar));
         else
            ylims_tmp(iy,ix) = max(n1(:,ipar));
         end
         subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
         hold on;box on;

         lim_2 = [min(m(:,ipar));min(m(:,ipar));
                  lim(:,ipar); max(m(:,ipar)); max(m(:,ipar))];
         n1_2 = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
         fill(lim_2,n1_2,[0.8,0.8,0.8]);

         set(gca,'layer','top')
         ipar = ipar +1;
         isubfig = isubfig + 1;
      end;
      ipanel = ipanel + 1;
   end
end
end

%ylims = max(max(ylims_tmp));
ylims = max(ylims_tmp);
maxtmp= max(ylims(3:4));
ylims(3:4) = maxtmp;


ipar = 1;
isubfig = 1;
ipanel = 1;
for iy = 1:ny
for ix = 1:nx
   if(ipar <= npar)
      if(ix == 3 && iy >= 2) 
         isubfig = isubfig + 1;
      elseif(ipanel == nx*(nlay+1)-4) 
         isubfig = isubfig + 1;
      else
         subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
         set(gca,'YTick',[],'FontSize',14);
         set(gca,'XTick',[F(ix).xticks],'FontSize',14);
         set(gca,'XTickLabel',[],'FontSize',14);
         if(iy == 1 && ix == 3) 
            set(gca,'XTick',[F(ix).xticks],'FontSize',14);
            set(gca,'XTickLabel',[F(ix).xticks],'FontSize',14);
            xlabel(xparname(ipanel));
         elseif(iy == nlay && ix == 1)
            set(gca,'XTick',[F(ix).xticks],'FontSize',14);
            set(gca,'XTickLabel',[F(ix).xticks],'FontSize',14);
            xlabel(xparname(ipanel));
         elseif(iy > nlay)
            set(gca,'XTick',[F(ix).xticks],'FontSize',14);
            set(gca,'XTickLabel',[F(ix).xticks],'FontSize',14);
            xlabel(xparname(ipanel));
         end 
         set(gca,'XLim',[plo(ipanel) phi(ipanel)],'FontSize',14);
         set(gca,'YLim',[0 ylims(ix)+.05*ylims(ix)]);
         ipar = ipar +1;
         isubfig = isubfig + 1;
      end;
      ipanel = ipanel + 1;
   end
end
end

if(i_save == 1)
   exportfig(gc1,plotfile,opts);
end

if(i_save == 1)
    save(hpdfile,'nf_int','-ASCII');
end

return;
