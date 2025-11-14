%----------------------------------------------------------
% Compute Prior from PPD
%----------------------------------------------------------
function [] = compute_prior();

i_rot = 1;
i_save = 0;

Tstar = 1;
Tstarc = 1;
NLFIX = 0
if(NLFIX>0)
    NPFIX = NLFIX*4;
else
    NPFIX = 7;
end
nx = 4
ny = 7
nsubfig = nx*ny;
xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 30;
nbin2 = 10;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
% INPUT:
%ppd     = 'sim_C_1_7_40_1lay_sample.mat';
%ppd2    = 'sim_C_1_7_40_1lay_sample.mat';
%data    = 'sim_C_1_7_40_1lay.mat';

ppd     = 'x_jan_1_sph_ci2_sample.mat';
ppd2    = 'x_jan_1_sph_ci2_sample.mat';
data    = 'x_jan_1.mat';

%ppd     = 'x_s13_2_10_50_3lay_sample.mat';
%ppd2    = 'x_s13_2_10_50_3lay_sample.mat';
%data    = 'x_s13_2_10_50_3lay.mat';

mtrue = [0.30, 1480, 1.30, 0.01, ...
               1560, 1.70, 0.01 ]
%        0.35, 1510, 1.65, 0.30, ...
%        1.20, 1500, 1.55, 0.30, ...
%        0.40, 1700, 1.90, 0.30, ...
%         0.60, 1500, 1.55, 0.30, ...
%         0.15, 1480, 1.50, 0.30, ...
plo = [0.0, 1450, 1.2, 0.0,...
            1450, 1.2, 0.0];
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
phi = [2.0, 1750, 2.1, 1.0, ...
            1750, 2.1, 1.0];
%       2.0, 1750, 2.2, 1.0,...
%       2.0, 1750, 2.2, 1.0,...
%       2.0, 1750, 2.2, 1.0,...
%       2.0, 1750, 2.0, 1.0,...
%       2.0, 1750, 2.0, 1.0,...
%       6.0, 1750, 2.2, 1.0,...
%       6.0, 1750, 2.2, 1.0,...
%       2.0, 1750, 2.0, 1.0,...

if(i_rot == 1)
%    data2   = 'sim_A_1.mat';
%    data2   = 'sim_A_1.mat';
    data2   = 'x_jan_1.mat';
end

% OUTPUT:
if(i_save == 1)
   prior         = 'sim_C_1_7_40_1lay_prior.txt';
   priormat      = 'sim_C_1_7_40_1lay_prior.mat';
   mapfile       = 'sim_C_1_7_40_1lay_map.txt';
   meanfile      = 'sim_C_1_7_40_1lay_mean.txt';
   hpdfile       = 'sim_C_1_7_40_1lay_hpds.txt';
   plotfile(1,:) = 'sim_C_1_7_40_1lay_marg_a.eps';
   plotfile(2,:) = 'sim_C_1_7_40_1lay_marg_b.eps';
   plotfile(3,:) = 'sim_C_1_7_40_1lay_marg_c.eps';

%    prior         = 'x_s21_1_10-50_4lay_ci_prior.txt';
%    priormat      = 'x_s21_1_10-50_4lay_ci_prior.mat';
%    mapfile       = 'x_s21_1_10-50_4lay_ci_map.txt';
%    meanfile      = 'x_s21_1_10-50_4lay_ci_mean.txt';
%    hpdfile       = 'x_s21_1_10-50_4lay_ci_hpds.txt';
%    plotfile(1,:) = 'x_s21_1_10-50_4lay_ci_marg_a.eps';
%    plotfile(2,:) = 'x_s21_1_10-50_4lay_ci_marg_b.eps';
%    plotfile(3,:) = 'x_s21_1_10-50_4lay_ci_marg_c.eps';

%    prior         = 'x_s16_1_7_40_3layrg_prior.txt';
%    priormat      = 'x_s16_1_7_40_3layrg_prior.mat';
%    mapfile       = 'x_s16_1_7_40_3layrg_map.txt';
%    meanfile      = 'x_s16_1_7_40_3layrg_mean.txt';
%    hpdfile       = 'x_s16_1_7_40_3layrg_hpds.txt';
%    plotfile(1,:) = 'x_s16_1_7_40_3layrg_marg_a.eps';
%    plotfile(2,:) = 'x_s16_1_7_40_3layrg_marg_b.eps';
%    plotfile(3,:) = 'x_s16_1_7_40_3layrg_marg_c.eps';

%    prior         = 'x_s13_2_10_50_3lay_prior.txt';
%    priormat      = 'x_s13_2_10_50_3lay_prior.mat';
%    mapfile       = 'x_s13_2_10_50_3lay_map.txt';
%    meanfile      = 'x_s13_2_10_50_3lay_mean.txt';
%    hpdfile       = 'x_s13_2_10_50_3lay_hpds.txt';
%    plotfile(1,:) = 'x_s13_2_10_50_3lay_marg_a.eps';
%    plotfile(2,:) = 'x_s13_2_10_50_3lay_marg_b.eps';
%    plotfile(3,:) = 'x_s13_2_10_50_3lay_marg_c.eps';

%    prior         = 'sim_A_2_ci_prior.txt';
%    priormat      = 'sim_A_2_ci_prior.mat';
%    mapfile       = 'sim_A_2_ci_map.txt';
%    meanfile      = 'sim_A_2_ci_mean.txt';
%    hpdfile       = 'sim_A_2_ci_hpds.txt';
%    plotfile(1,:) = 'sim_A_2_ci_marg_a.eps';
%    plotfile(2,:) = 'sim_A_2_ci_marg_b.eps';
%    plotfile(3,:) = 'sim_A_2_ci_marg_c.eps';

%    prior         = 'tmp_ci_prior.txt';
%    priormat      = 'tmp_ci_prior.mat';
%    mapfile       = 'tmp_ci_map.txt';
%    meanfile      = 'tmp_ci_mean.txt';
%    hpdfile       = 'tmp_ci_hpds.txt';
%    plotfile(1,:) = 'tmp_ci_marg_a.eps';
%    plotfile(2,:) = 'tmp_ci_marg_b.eps';
%    plotfile(3,:) = 'tmp_ci_marg_c.eps';
end

load(data);
load(ppd);
m = A(:,2:NPFIX+1);
%m = A(:,2:end-1);
A2 = A;
E = A2(:,1);

load(ppd2);
mc = A(:,2:NPFIX+1);
%mc = A(:,2:end-1);
Ec = A(:,1);

size(m)
size(mc)

lo = F(1).minlim(1:end-3);
hi = F(1).maxlim(1:end-3);

lob = F(1).minlim(end-2:end);
hib = F(1).maxlim(end-2:end);

F(1).xticks = [0.5, 1, 1.5];
F(2).xticks = [1500, 1600, 1700];
F(3).xticks = [1.4, 1.6, 1.8, 2.0];
F(4).xticks = [0, 0.5, 1.0 ];
xparname = [{'h1'},{'c1'},{'r1'},{'a1'},...
            {'h2'},{'c2'},{'r2'},{'a2'},...
            {'h3'},{'c3'},{'r3'},{'a3'},...
            {'h4'},{'c4'},{'r4'},{'a4'},...
            {'h5'},{'c5'},{'r5'},{'a5'},...
            {'h6'},{'c6'},{'r6'},{'a6'},...
            {'h7'},{'c7'},{'r7'},{'a7'}];
xparname2 = [{'m1'},{'m2'},{'m3'},{'m4'}];

if(i_rot == 1)
    load(data2);
    minlim = F(1).minlim(1:4);
    maxlim = F(1).maxlim(1:4);
    load(data);
    minlim = [minlim F(1).minlim];
    maxlim = [maxlim F(1).maxlim];
    minlim
    maxlim
end

[cost_min i_min] = min(E)
xmap = m(i_min,:);
if(i_save == 1)
    save(mapfile,'xmap','-ASCII');
end
xmap

for ipar = 1:size(m,2)

    xmean(ipar) = mean(m(:,ipar));

end;
if(i_save == 1)
     save(meanfile,'xmean','-ASCII');
end;

[EE idxE] = min(E);
mstart = m(idxE,:);

npar = size(m,2)
nmods = size(m,1)
nmodsc = size(mc,1)
nfig = ceil(npar/nsubfig);

%----------------------------------------------------------
%  Calc correlation matrix:
%----------------------------------------------------------
if(i_rot == 1)
%m2 = m;
%mm = zeros(size(m2));
%for i = 1:npar
%  m2(:,i) = (m2(:,i) - minlim(i))/(maxlim(i)-minlim(i));
%  mm(:,i) = m2(:,i) - sum(m2(:,i))/nmods;
%end
mm = m;

mcov = zeros(npar,npar);
mcor = zeros(size(mcov));
mcorp = zeros(size(mcov));
for i = 1:npar
  for j = 1:npar
    mcov(i,j) = sum(mm(:,i).*mm(:,j));
  end
end
mcov = mcov/nmods;

for i = 1:npar
  for j = 1:npar
    mcor(i,j) = mcov(i,j)/sqrt(mcov(i,i)*mcov(j,j));
  end
end
[U,S,V] = svd(mcov);
%U = U(1:NPFIX,1:NPFIX);
%mcor;
end


%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

for ipar = 1:npar

    [n1(:,ipar),lim(:,ipar)] = hist_tstar(m(:,ipar),E,nbin1,Tstar);
%    [n1(:,ipar),lim(:,ipar)] = hist(m(:,ipar),nbin1);
    n1(:,ipar) = n1(:,ipar) * (phi(ipar) - plo(ipar));
    nf_int(ipar,:) = hpd(m(:,ipar),100,95);
    [lim(1,ipar) lim(end,ipar)];

    [nc(:,ipar),limc(:,ipar)] = hist_tstar(mc(:,ipar),Ec,nbin1,Tstarc);
%    [nc(:,ipar),limc(:,ipar)] = hist(mc(:,ipar),nbin1);
    nc(:,ipar) = nc(:,ipar) * (phi(ipar) - plo(ipar));
    nf_intc(ipar,:) = hpd(mc(:,ipar),100,95);

end

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ipar = 1;
for ifig = 1:nfig
    gc1=figure(ifig);
    isubfig = 1;
    while (isubfig <= nsubfig & ipar <= npar)

        if(ipar == npar-2) 
            isubfig = isubfig + 1;
        end;
        subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
        hold on;

        lim_2 = [min(m(:,ipar));min(m(:,ipar));
                 lim(:,ipar); max(m(:,ipar)); max(m(:,ipar))];
        n1_2 = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
        fill(lim_2,n1_2,[0.8,0.8,0.8]);

        lim_2 = [min(mc(:,ipar)); min(mc(:,ipar)); 
                 limc(:,ipar)+mean(diff(limc(:,ipar)))/2];
        n1_2 = [0; nc(:,ipar); 0];
%        stairs(lim_2,n1_2,'-r');
        set(gca,'layer','top')
        box on;
        set(gca,'YTick',[],'FontSize',14);
        set(gca,'XLim',[plo(ipar) phi(ipar)],'FontSize',14);
%        plot([mtrue(ipar) mtrue(ipar)],[0 ylims(2)],'-k')
%        plot([xmap(ipar) xmap(ipar)],[0 ylims(2)],'--k')

        xlabel(xparname(ipar));
        if(ipar > npar-3) 
          xlabel(xparname(ipar+1));
        end
        set(gca,'YLim',[0 20]);
        ipar = ipar +1;
        isubfig = isubfig + 1;

    end
    if(i_save == 1)
        exportfig(gc1,plotfile(ifig,:),opts);
    end
end

if(i_rot == 1)

fprintf(1,'rotating...\n')
for imod = 1:nmods

    m_rot(imod,1:NPFIX) = U'*m(imod,1:NPFIX)';

end

save tmp.mat m_rot

fprintf(1,'done rotating!\n')

for ipar = 1:NPFIX

    [n2(:,ipar),lim2(:,ipar)] = hist(m_rot(:,ipar),nbin2);
    width(ipar) = lim2(end,ipar)-lim2(1,ipar);
    n2(:,ipar) = n2(:,ipar)/(nmods);
    lim2(end,ipar) = max(m_rot(:,ipar));
    lim2(1,ipar) = min(m_rot(:,ipar));
    nf_int(ipar,:) = hpd(m_rot(:,ipar),100,99.9);

%    tmp1=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2-...
%         abs(lim2(1,ipar)-lim2(end,ipar));
%    tmp2=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2+...
%         abs(lim2(1,ipar)-lim2(end,ipar));
%    tmp3(:,ipar)=[tmp1:abs(tmp1-tmp2)/9:tmp2];

end
%lim3 = lim2;
%lim2 = tmp3;

for i = 1:NPFIX
    m3(:,i) = (rand(10000,1)-0.5)*(abs(lim2(1,i)-lim2(end,i)))+...
              min(lim2(:,i))+abs(max(lim2(:,i))-min(lim2(:,i)))/2;
end

for i = 1:10000
    m4(i,:) = U*m3(i,:)';
end;

for i = 1:NPFIX
    lim(1,i) = min(m4(:,i));
    lim(end,i) = max(m4(:,i));
end;

if(i_save == 1)
    save(prior,'n1','lim','mstart','n2','lim2','U','-ASCII');
    save(priormat,'n1','lim','mstart','n2','lim2','U');
end
n2
nfig2 = 2 * nfig;
ipar = 1;
for ifig = nfig+1:nfig2
%for ifig = 1:nfig
    gc1=figure(ifig);
    isubfig = 1;
    while (isubfig <= nsubfig & ipar <= NPFIX)

        subplot(4,4,isubfig);hold on;box on;
%        subplot('Position',[llw(ipar) llh(2) spw sph]);
%        stairs(lim3(:,ipar),n2(:,ipar),'--k');
        stairs(lim2(:,ipar),n2(:,ipar),'-k');
        ylims = get(gca,'Ylim');
%        xlabel(xparname2(ipar));
        set(gca,'YTick',[],'FontSize',14);
        ipar = ipar +1;
        isubfig = isubfig + 1;

    end
end
else % else if rotate if

for ipar = 1:NPFIX

    [n2(:,ipar),lim2(:,ipar)] = hist(m(:,ipar),nbin2);
    width(ipar) = lim2(end,ipar)-lim2(1,ipar);
    n2(:,ipar) = n2(:,ipar)/(nmods);
    lim2(end,ipar) = max(m(:,ipar));
    lim2(1,ipar) = min(m(:,ipar));
    nf_int(ipar,:) = hpd(m(:,ipar),100,99.9);
    
    tmp1=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2-...
         abs(lim2(1,ipar)-lim2(end,ipar));
    tmp2=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2+...
         abs(lim2(1,ipar)-lim2(end,ipar));
    tmp3=[tmp1:abs(tmp1-tmp2)/9:tmp2];

end
lim2 = tmp3;

if(i_save == 1 && NLFIX >1)
    save(prior,'mstart','n2','lim2','-ASCII');
    save(priormat,'mstart','n2','lim2');
end


end % end if rotate if

if(i_save == 1)
    save(hpdfile,'nf_int','-ASCII');
end

return;
