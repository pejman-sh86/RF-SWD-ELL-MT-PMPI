function [] = plot_rjmh_profs(filename);

corefile    = 'core.mat';
%% DELTA SITE: 
hmax    = 320.;
%hmax    = 45.;
order = 1;
icore   = 1;
pmin = [   0  100 1.2]';
pmax = [1200 3000 2.7]';

load(filename);
if(icore == 1)
    load(corefile);
end;

NPROF = length(A);
BURNIN = ceil(NPROF/3);
A = A(BURNIN:end,:);


k = A(:,4);

NFP = (k*4)+3;
m = A(:,5:end-5-order);

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
set(gca,'Fontsize',12,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
set(gca,'Fontsize',12,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;
set(gca,'Fontsize',12,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse');

for i=1:1000
   mm = m(i,1:k(i)*4+3);
   subplot(h1);
   plprof(mm,hmax,'--k',1,0);
   subplot(h2);
   plprof(mm,hmax,'--k',2,0);
   subplot(h3);
   plprof(mm,hmax,'--k',3,0);
   mm
end;

return;
