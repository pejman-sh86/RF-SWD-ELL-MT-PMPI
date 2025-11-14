%----------------------------------------------------------
% Compute Prior from PPD
%----------------------------------------------------------
function [] = compute_prior_tt();

i_save = 1;
i_saveplot = 1;
ndat = 1;

nx = 3;
ny = 4;
xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

% INPUT:
ppd    = 'sample_tt_int.mat';
ppd2     = 'sample_real.mat';
ppd3     = 'sample_tt_int_old.mat';

% OUTPUT:
if(i_save == 1)
    prior   = 'picked_tt_int_covstat_prior.txt';
    mapfile = 'picked_tt_int_covstat_map.txt';
    meanfile = 'picked_tt_int_covstat_mean.txt';
    stdfile = 'picked_tt_int_covstat_std.txt';
    hpdfile = 'picked_tt_int_covstat_hpds.txt';
end
if(i_saveplot == 1)
    plotfile(1,:) = 'picked_tt_covstat_marg_a.eps';
    plotfile(2,:) = 'picked_tt_covstat_marg_b.eps';
    plotfile(3,:) = 'picked_tt_covstat_marg_c.eps';
%    plotfile(1,:) = 'picked_tt_int_covstat_marg_a.eps';
%    plotfile(2,:) = 'picked_tt_int_covstat_marg_b.eps';
%    plotfile(3,:) = 'picked_tt_int_covstat_marg_c.eps';
%    plotfile(1,:) = 'picked_tt2_compare_marg_a.eps';
%    plotfile(2,:) = 'picked_tt2_compare_marg_b.eps';
%    plotfile(3,:) = 'picked_tt2_compare_marg_c.eps';
end

load(ppd);
m = A(:,2:end);
%m = [A(:,2:7) A(:,14:17)];
nlay = 3;
npar = size(m,2)
A1 = A;
min_A1 = min(A1);
max_A1 = max(A1);
if(ndat > 1)
    load(ppd2);
%    mc = A(:,2:end);
    mc = [A(:,2:7) A(:,14:17)];
    A2 = A;
    min_A2 = min(A2)
    max_A2 = max(A2)
end
if(ndat>2)
   load(ppd3);
   mcc = A(:,2:end);
   A3 = A;
   min_A3 = min(A3)
   max_A3 = max(A3)
end

mtrue=[3.7958381        1540.46     0.93985244      1508.3639       9.428378      1506.6589      0.24088261     0.48829459     0.45668536      0.4973452];
xlo = [0.5, 1450.,...
       0.5, 1450.,...
       0.5, 1450.,...
       0.0, 0.0, 0.0, 0.0];
xhi = [4., 1750.,... 
       4., 1750.,... 
       4., 1750.,... 
       1.0, 1.5, 1.5, 1.5]; 
if(length(m(1,:)) == 10)
  m(:,(nlay*2)+1) = m(:,(nlay*2)+1)-xhi((nlay*2)+1)/2;
  m(:,(nlay*2)+2) = m(:,(nlay*2)+2)-xhi((nlay*2)+2)/2;
  m(:,(nlay*2)+3) = m(:,(nlay*2)+3)-xhi((nlay*2)+3)/2;
  m(:,(nlay*2)+4) = m(:,(nlay*2)+4)-xhi((nlay*2)+4)/2;
  if(ndat > 1)
    mc(:,(nlay*2)+1) = mc(:,(nlay*2)+1)-xhi((nlay*2)+1)/2;
    mc(:,(nlay*2)+2) = mc(:,(nlay*2)+2)-xhi((nlay*2)+2)/2;
    mc(:,(nlay*2)+3) = mc(:,(nlay*2)+3)-xhi((nlay*2)+3)/2;
    mc(:,(nlay*2)+4) = mc(:,(nlay*2)+4)-xhi((nlay*2)+4)/2;
  end
  if(ndat > 2)
    mcc(:,(nlay*2)+1) = mcc(:,(nlay*2)+1)-xhi((nlay*2)+1)/2;
    mcc(:,(nlay*2)+2) = mcc(:,(nlay*2)+2)-xhi((nlay*2)+2)/2;
    mcc(:,(nlay*2)+3) = mcc(:,(nlay*2)+3)-xhi((nlay*2)+3)/2;
    mcc(:,(nlay*2)+4) = mcc(:,(nlay*2)+4)-xhi((nlay*2)+4)/2;
  end

end
xhi(:,(nlay*2)+1) = xhi(:,(nlay*2)+1)-xhi((nlay*2)+1)/2;
xhi(:,(nlay*2)+2) = xhi(:,(nlay*2)+2)-xhi((nlay*2)+2)/2;
xhi(:,(nlay*2)+3) = xhi(:,(nlay*2)+3)-xhi((nlay*2)+3)/2;
xhi(:,(nlay*2)+4) = xhi(:,(nlay*2)+4)-xhi((nlay*2)+4)/2;

xlo(:,(nlay*2)+1) = xlo(:,(nlay*2)+1)-xhi((nlay*2)+1);
xlo(:,(nlay*2)+2) = xlo(:,(nlay*2)+2)-xhi((nlay*2)+2);
xlo(:,(nlay*2)+3) = xlo(:,(nlay*2)+3)-xhi((nlay*2)+3);
xlo(:,(nlay*2)+4) = xlo(:,(nlay*2)+4)-xhi((nlay*2)+4);
%lo = [1.6 1450 1.60 0.00];
%hi = [1.8 1490 1.75 0.15];

%lob = [1560 1.65 0.00];
%hib = [1590 1.80 0.15];

%xparname = [{'h_1 (m)'},{'c_1 (m/s)'},...
%            {'h_2 (m)'},{'c_2 (m/s)'},...
%            {'h_3 (m)'},{'c_3 (m/s)'},...
%            {'h_4 (m)'},{'c_4 (m/s)'},...
%            {'h_5 (m)'},{'c_5 (m/s)'},...
%            {'h_6 (m)'},{'c_6 (m/s)'},...
%            {'a_1 (ms)'},{'a_2 (ms)'},{'a_3 (ms)'},{'a_4 (ms)'},...
%                   {'a_5 (ms)'},{'a_6 (ms)'},{'a_7 (ms)'}];
xparname = [{'h1'},{'c1'},...
            {'h2'},{'c2'},...
            {'h3'},{'c3'},...
            {'a1'},{'a2'},{'a3'},{'a4'}];
%            {'h4'},{'c4'},...
%           {'h5'},{'c5'},...
%           {'h6'},{'c6'},...
%            {'a1'},{'a2'},{'a3'},{'a4'},...
%                   {'a5'},{'a6'},{'a7'}];
%xparname2 = [{'m1'},{'m2'},{'m3'},{'m4'}];

nbin = 30;
nsubfig = nx*ny;
nfig = ceil(npar/nsubfig);
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);


[cost_min i_min] = min(A1(:,1));
xmap = m(i_min,:);
if(i_save == 1)
    save(mapfile,'xmap','-ASCII');
end

for ipar = 1:size(m,2)

    xmean(ipar) = mean(m(:,ipar));
    xstd(ipar) = std(m(:,ipar));

end;
if(i_save == 1)
     save(meanfile,'xmean','-ASCII');
     save(stdfile,'xstd','-ASCII');
end;

conv = min(A1(:,end));
[E idxE] = min(A1(:,1));
mstart = m(idxE,:);

nmods = size(m,1);
if(ndat>1)
    nmodsc = size(mc,1);
end
if(ndat>2)
    nmodscc = size(mcc,1);
end

%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

p_offset1 = zeros(1,npar);p_offset2 = zeros(1,npar);p_offset3 = zeros(1,npar);

for ipar = 1:npar

    [n1(:,ipar),lim(:,ipar)] = hist(m(:,ipar),nbin);
    width(ipar) = lim(end,ipar)-lim(1,ipar);
    n1(:,ipar) = n1(:,ipar)/(width(ipar)/nbin*nmods);
    nf_int(ipar,:) = hpd(m(:,ipar),100,95);
    lim_2(:,ipar) = [min(m(:,ipar));min(m(:,ipar));lim(:,ipar);max(m(:,ipar));max(m(:,ipar))];
    n1_2(:,ipar) = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
    p_offset1(ipar) = 1.2*max(n1_2(:,ipar));

  if(ndat > 1)
    [nc(:,ipar),limc(:,ipar)] = hist(mc(:,ipar),nbin);
    widthc(ipar) = limc(end,ipar)-limc(1,ipar);
    nc(:,ipar) = nc(:,ipar)/(widthc(ipar)/nbin*nmodsc);
    nf_intc(ipar,:) = hpd(mc(:,ipar),100,95);
    lim_2c(:,ipar) = [min(mc(:,ipar));min(mc(:,ipar));limc(:,ipar); max(mc(:,ipar)); max(mc(:,ipar))];
    n1_2c(:,ipar) = [0; nc(1,ipar); nc(:,ipar); nc(end,ipar); 0];
    p_offset2(ipar) = 1.2*max(n1_2c(:,ipar));
%    p_offset2 = p_offset1;
  end

  if(ndat > 2)
    [ncc(:,ipar),limcc(:,ipar)] = hist(mcc(:,ipar),nbin);
    widthcc(ipar) = limcc(end,ipar)-limcc(1,ipar);
    ncc(:,ipar) = ncc(:,ipar)/(widthcc(ipar)/nbin*nmodscc);
    nf_intcc(ipar,:) = hpd(mcc(:,ipar),100,95);
    lim_2cc(:,ipar) = [min(mcc(:,ipar)); min(mcc(:,ipar)); limcc(:,ipar);max(mcc(:,ipar));max(mcc(:,ipar))];
    n1_2cc(:,ipar) = [0; ncc(1,ipar); ncc(:,ipar); ncc(end,ipar); 0];
    p_offset3(ipar) = 1.2*max(n1_2cc(:,ipar));
%    p_offset3 = p_offset1;
  end

end

idx = [1 2 7 3 4 8 5 6 9 10];
p_offset = max([p_offset1; p_offset2; p_offset3]);
ipar = idx(1);
for ifig = 1:nfig
    gc1=figure(ifig);
    set(gc1,'PaperUnits','inches','Units','inches')
    isubfig = 1;
    while (isubfig <= nsubfig & ipar < npar)
        ipar = idx(isubfig);

        if(ipar == npar)
            subplot('Position',[loc(1,isubfig+nx-1) loc(2,isubfig+nx-1) spw sph]);
        else
            subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
        end
        hold on;
%
%       DATA SET 1
%
        fill(lim_2(:,ipar),n1_2(:,ipar),[0.8,0.8,0.8]);
        
        
%
%       DATA SET 2
%
        if(ndat > 1)
            stairs(lim_2c(:,ipar),n1_2c(:,ipar)+p_offset(ipar),'-k');
            plot([xlo(ipar) xhi(ipar)],[p_offset(ipar) ... 
                  p_offset(ipar)],'-k','LineWidth',0.5)
        end

%
%       DATA SET 3
%
        if(ndat > 2)
            stairs(lim_2cc(:,ipar),n1_2cc(:,ipar)+p_offset(ipar)+ ...
                   p_offset(ipar),'-k');
            plot([xlo(ipar) xhi(ipar)],[p_offset(ipar)+p_offset(ipar) ... 
                  p_offset(ipar)+ p_offset(ipar)],'-k','LineWidth',0.5)
        end

        set(gca,'layer','top')
        box on;
        set(gca,'YTick',[],'FontSize',14);
        ylims(1) = 0;
%        ylims(2) = max([max(n1(:,ipar)) max(nc(:,ipar))]);
        ylims(2) = p_offset(ipar);
        if(ndat > 1)
            ylims(2) = ylims(2) + p_offset(ipar);
        end
        if(ndat > 2)
            ylims(2) = ylims(2) + p_offset(ipar);
        end
        xlabel(xparname(ipar));
        set(gca,'YLim',ylims);
        set(gca,'XLim',[xlo(ipar) xhi(ipar)]);
%        set(gca,'XLim',[-0.75 0.75]);
%        set(gca,'XLim',[-0.5 0.5]);
        isubfig = isubfig + 1;

    end
    if(i_saveplot == 1)
%        saveas(gca,plotfile(ifig,:),'epsc2');
        exportfig(gc1,plotfile(ifig,:),opts);
    end
end

for ipar = 1:npar

    [n2(:,ipar),lim2(:,ipar)] = hist(m(:,ipar),nbin);
    width(ipar) = lim2(end,ipar)-lim2(1,ipar);
    n2(:,ipar) = n2(:,ipar)/(nmods);
    lim2(end,ipar) = max(m(:,ipar));
    lim2(1,ipar) = min(m(:,ipar));
    nf_int(ipar,:) = hpd(m(:,ipar),100,99.99);

end

if(i_save == 1)
    save(prior,'n1','lim','mstart','n2','lim2','-ASCII');
end
lim2(1,:)
lim2(end,:)

if(i_save == 1)
    save(hpdfile,'nf_int','-ASCII');
end

return;
