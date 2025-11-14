%----------------------------------------------------------
% Compute Prior from PPD
%----------------------------------------------------------
function [] = compute_prior();

i_rot = 0;
i_save = 1;

Tstar = 1;
Tstarc = 1;
NLFIX = 0
if(NLFIX>0)
    NPFIX = NLFIX*4+1;
else
    NPFIX = 20;
end
nx = 5
ny = 5
nsubfig = nx*ny;
xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 20;
nbin2 = 10;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
% INPUT:
ppd     = 'sim_C_1_7_40_4layrg_sample.mat';
ppd2    = 'sim_C_1_7_40_4layrg_sample.mat';
data    = 'sim_C_1_7_40_4layrg.mat';
mtrue = [0.30, 1480, 1.30, 1.6, 0.01, ...
         0.35, 1510, 1.65, 0.30, ...
         0.60, 1500, 1.55, 0.30, ...
         0.40, 1700, 1.90, 0.30, ...
               1560, 1.70, 0.01 ]
%         0.40, 1700, 1.90, 0.30, ...
%         0.60, 1500, 1.55, 0.30, ...

%ppd     = 'x_s13_2_10_50_4layrg_sample.mat';
%ppd2    = 'x_s13_2_10_50_4layrg_sample.mat';
%data    = 'x_s13_2_10_50_4layrg.mat';
%ppd     = 'x_s16_1_15_40_3layrg_sample.mat';
%ppd2    = 'x_s16_1_15_40_3layrg_sample.mat';
%data    = 'x_s16_1_15_40_3layrg.mat';
%ppd     = 'x_s21_1_10-50_4layrg_b_ci_samplec.mat';
%ppd2    = 'x_s21_1_10-50_4layrg_b_ci_samplec.mat';
%data    = 'x_s21_1_10-50_4layrg_b.mat';
if(i_rot == 1)
    data2   = 'x_s16_1_5-40_4layrg2.mat';
%    data    = 'x_jan_2.mat';
%    data2   = 'x_jan_1.mat';
end

% OUTPUT:
if(i_save == 1)

   prior         = 'sim_C_1_7_40_4layrg_prior.txt';
   priormat      = 'sim_C_1_7_40_4layrg_prior.mat';
   mapfile       = 'sim_C_1_7_40_4layrg_map.txt';
   meanfile      = 'sim_C_1_7_40_4layrg_mean.txt';
   hpdfile       = 'sim_C_1_7_40_4layrg_hpds.txt';
   plotfile(1,:) = 'sim_C_1_7_40_4layrg_marg_a.eps';
   plotfile(2,:) = 'sim_C_1_7_40_4layrg_marg_b.eps';
   plotfile(3,:) = 'sim_C_1_7_40_4layrg_marg_c.eps';

%    prior         = 'x_s21_1_10-50_4layrg_b_ci_prior.txt';
%    priormat      = 'x_s21_1_10-50_4layrg_b_ci_prior.mat';
%    mapfile       = 'x_s21_1_10-50_4layrg_b_ci_map.txt';
%    meanfile      = 'x_s21_1_10-50_4layrg_b_ci_mean.txt';
%    hpdfile       = 'x_s21_1_10-50_4layrg_b_ci_hpds.txt';
%    plotfile(1,:) = 'x_s21_1_10-50_4layrg_b_ci_marg_a.eps';
%    plotfile(2,:) = 'x_s21_1_10-50_4layrg_b_ci_marg_b.eps';
%    plotfile(3,:) = 'x_s21_1_10-50_4layrg_b_ci_marg_c.eps';

%   prior         = 'x_s16_1_15_40_3layrg_prior.txt';
%   priormat      = 'x_s16_1_15_40_3layrg_prior.mat';
%   mapfile       = 'x_s16_1_15_40_3layrg_map.txt';
%   meanfile      = 'x_s16_1_15_40_3layrg_mean.txt';
%   hpdfile       = 'x_s16_1_15_40_3layrg_hpds.txt';
%   plotfile(1,:) = 'x_s16_1_15_40_3layrg_marg_a.eps';
%   plotfile(2,:) = 'x_s16_1_15_40_3layrg_marg_b.eps';
%   plotfile(3,:) = 'x_s16_1_15_40_3layrg_marg_c.eps';

%    prior         = 'x_s13_2_10_50_4layrg_prior.txt';
%    priormat      = 'x_s13_2_10_50_4layrg_prior.mat';
%    mapfile       = 'x_s13_2_10_50_4layrg_map.txt';
%    meanfile      = 'x_s13_2_10_50_4layrg_mean.txt';
%    hpdfile       = 'x_s13_2_10_50_4layrg_hpds.txt';
%    plotfile(1,:) = 'x_s13_2_10_50_4layrg_marg_a.eps';
%    plotfile(2,:) = 'x_s13_2_10_50_4layrg_marg_b.eps';
%    plotfile(3,:) = 'x_s13_2_10_50_4layrg_marg_c.eps';

%    prior         = 'x_jan_2_sph_ciw_prior.txt';
%    priormat      = 'x_jan_2_sph_ciw_prior.mat';
%    mapfile       = 'x_jan_2_sph_ciw_map.txt';
%    meanfile      = 'x_jan_2_sph_ciw_mean.txt';
%    hpdfile       = 'x_jan_2_sph_ciw_hpds.txt';
%    plotfile(1,:) = 'x_jan_2_sph_ciw_marg_a.eps';
%    plotfile(2,:) = 'x_jan_2_sph_ciw_marg_b.eps';
%    plotfile(3,:) = 'x_jan_2_sph_ciw_marg_c.eps';

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

plo = [0.0, 1450, 1.2, 1.2, 0.0,...
       0.0, 1450, 1.2, 0.0,...
       0.0, 1450, 1.2, 0.0,...
       0.0, 1450, 1.2, 0.0,...
            1450, 1.2, 0.0];
%       0.1, 1450, 1.2, 0.0,...
%       0.1, 1450, 1.2, 0.0,...
phi = [2.5, 1750, 2.2, 2.4, 1.0, ...
       2.5, 1750, 2.2, 1.0,...
       2.5, 1750, 2.2, 1.0,...
       2.5, 1750, 2.2, 1.0,...
            1750, 2.2, 1.0];
%       2.5, 1750, 2.4, 1.0,...
%       2.5, 1750, 2.4, 1.0,...
xparname = [{'h1'},{'c1'},{'r1t'},{'r1b'},{'a1'},...
            {'h2'},{'c2'},{'r2'},{'a2'},...
            {'h3'},{'c3'},{'r3'},{'a3'},...
            {'h4'},{'c4'},{'r4'},{'a4'},...
            {'h5'},{'c5'},{'r5'},{'a5'},...
            {'h6'},{'c6'},{'r6'},{'a6'},...
            {'h7'},{'c7'},{'r7'},{'a7'}];
xparname2 = [{'m1'},{'m2'},{'m3'},{'m4'}];

if(i_rot == 1)
    load(data2);
    minlim = F(1).minlim(1:5);
    maxlim = F(1).maxlim(1:5);
    load(data);
    minlim = [minlim F(1).minlim];
    maxlim = [maxlim F(1).maxlim];
    minlim
    maxlim
end

[cost_min i_min] = min(A2(:,1))
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

conv = min(A2(:,end))
[EE idxE] = min(A2(:,1));
mstart = m(idxE,:);

npar = size(m,2)
nmods = size(m,1);
nmodsc = size(mc,1);
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
%    width(ipar) = lim(end,ipar)-lim(1,ipar);
%    n1(:,ipar) = n1(:,ipar)/(width(ipar)/nbin1*nmods);
    nf_int(ipar,:) = hpd(m(:,ipar),100,95);
    [lim(1,ipar) lim(end,ipar)];

    [nc(:,ipar),limc(:,ipar)] = hist_tstar(mc(:,ipar),Ec,nbin1,Tstarc);
%    [nc(:,ipar),limc(:,ipar)] = hist(mc(:,ipar),nbin1);
%    widthc(ipar) = limc(end,ipar)-limc(1,ipar);
%    nc(:,ipar) = nc(:,ipar)/(widthc(ipar)/nbin1*nmodsc);
    nf_intc(ipar,:) = hpd(mc(:,ipar),100,95);

end

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ipar = 1;
for ifig = 1:nfig
    gc1=figure(ifig);
    isubfig = 1;
    while (isubfig <= nsubfig & ipar <= npar)

        if(ipar == 8) 
            isubfig = isubfig + 1;
        elseif(ipar == 12) 
            isubfig = isubfig + 1;
        elseif(ipar == 16) 
            isubfig = isubfig + 1;
        elseif(ipar == 20) 
            isubfig = isubfig + 1;
        elseif(ipar == 24) 
            isubfig = isubfig + 1;
        elseif(ipar == 26) 
            isubfig = isubfig + 1;
        elseif(ipar == 27) 
            isubfig = isubfig + 1;
        end;
        subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
        hold on;

        lim_2 = [min(m(:,ipar));min(m(:,ipar));
                 lim(:,ipar); max(m(:,ipar)); max(m(:,ipar))];
        n1_2 = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
        fill(lim_2,n1_2,[0.8,0.8,0.8]);
%        stairs(lim_2,n1_2,'-k');

        lim_2 = [min(mc(:,ipar)); min(mc(:,ipar)); 
                 limc(:,ipar)+mean(diff(limc(:,ipar)))/2];
        n1_2 = [0; nc(:,ipar); 0];
%        stairs(lim_2,n1_2,'-r');
        set(gca,'layer','top')
        box on;
        set(gca,'YTick',[],'FontSize',14);
        set(gca,'XLim',[plo(ipar) phi(ipar)],'FontSize',14);
        ylims(1) = 0;
        ylims(2) = max([max(n1(:,ipar)) max(nc(:,ipar))]);
        ylims(2) = ylims(2)+ylims(2)/10;
%        plot([nf_int(ipar,1) nf_int(ipar,1)],[0 ylims(2)],':b')
%        plot([nf_int(ipar,2) nf_int(ipar,2)],[0 ylims(2)],':b')
%        plot([nf_intc(ipar,1) nf_intc(ipar,1)],[0 ylims(2)],':r')
%        plot([nf_intc(ipar,2) nf_intc(ipar,2)],[0 ylims(2)],':r')
        plot([mtrue(ipar) mtrue(ipar)],[0 ylims(2)],'-k')
%        plot([xmap(ipar) xmap(ipar)],[0 ylims(2)],'--k')

%        plot([xmap(ipar) xmap(ipar)],[0 ylims(2)],'--k')
%        if(isubfig ~= 9)
%        if(isubfig <= 12)
%            set(gca,'XTick',[]);
%        end
%        end
        xlabel(xparname(ipar));
        if(ipar > npar-3) 
          xlabel(xparname(ipar+1));
        end
        set(gca,'YLim',ylims);
%        set(gca,'XLim',[xlo(ipar) xhi(ipar)]);
        ipar = ipar +1;
        isubfig = isubfig + 1;

    end
    if(i_save == 1)
%        saveas(gca,plotfile(ifig,:),'epsc2');
        exportfig(gc1,plotfile(ifig,:),opts);
    end
end

if(i_rot == 1)

fprintf(1,'rotating...\n')
for imod = 1:nmods

    m_rot(imod,1:NPFIX) = U'*m(imod,1:NPFIX)';

end
fprintf(1,'done rotating!\n')

clear n1 lim;
for ipar = 1:NPFIX

    [n1(:,ipar),lim(:,ipar)] = hist(m(:,ipar),nbin2);
    n1(:,ipar) = n1(:,ipar)/(nmods);
    lim(end,ipar) = max(m(:,ipar));
    lim(1,ipar) = min(m(:,ipar));

    [n2(:,ipar),lim2(:,ipar)] = hist(m_rot(:,ipar),nbin2);
    n2(:,ipar) = n2(:,ipar)/(nmods);
    lim2(end,ipar) = max(m_rot(:,ipar));
    lim2(1,ipar) = min(m_rot(:,ipar));

%    tmp1=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2-...
%         abs(lim2(1,ipar)-lim2(end,ipar));
%    tmp2=min(lim2(:,ipar))+abs(lim2(1,ipar)-lim2(end,ipar))/2+...
%         abs(lim2(1,ipar)-lim2(end,ipar));
%    tmp3(:,ipar)=[tmp1:abs(tmp1-tmp2)/9:tmp2];

end
%lim3 = lim2;
%lim2 = tmp3;

%for i = 1:NPFIX
%    m3(:,i) = (rand(10000,1)-0.5)*(abs(lim2(1,i)-lim2(end,i)))+...
%              min(lim2(:,i))+abs(max(lim2(:,i))-min(lim2(:,i)))/2;
%end

%save bla

%for i = 1:10000
%    m4(i,:) = U*m3(i,:)';
%end;

%for i = 1:NPFIX
%    lim(1,i) = min(m4(:,i));
%    lim(end,i) = max(m4(:,i));
%end;

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
