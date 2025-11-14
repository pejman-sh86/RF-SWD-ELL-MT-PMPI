function [] = plot_model_hol_lay();

%set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);
set(0, 'DefaultFigurePaperPosition', [0 0 14 8]);

sourcefile = 'source_13_low.dat';
corefile  = 'core13.mat';
source = load(sourcefile);
source(:,1) = source(:,1) - source(1,1);
plotfile1 = 'sim_B_1_core.eps';
icore = 1;
iplot = 1;

if(icore == 1)
    load(corefile);
end;

hmax = 6;

%M(1).m = load(mapfile);
M(1).m = [0.10, 1490, 1.45, 0.30, ...
          0.30, 1550, 1.75, 0.25, ...
          0.35, 1580, 1.95, 0.30, ...
          3.50, 1540, 1.70, 0.10, ...
                1600, 1.75, 0.50]; 

pmin = [1450 1.2 0.0];
pmax = [1750 2.1 1.0];
         
npl = 4;

col = {'-k' '-r' '-r' '--k' '--r' ':k'};
col = char(col);

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',6,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf1 = figure(1);hold on; box on;
set(hf1, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box off;


for i = 1:size(M,2)

subplot(h1)
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
%plprof(M(i).m,hmax,'k',1);
plot((source(1:100,2)*100)+1650,(1540*source(1:100,1)-.1),'k','Linewidth',1);
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
box on;

subplot(h2)
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
%plprof(M(i).m,hmax,'k',2);
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
box on;

subplot(h3)
set(h3,'layer','top')
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)');
%plprof(M(i).m,hmax,'k',3);
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;

end;

barstep1 = 5;
barstep1r = 5;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'g','Linewidth',1);
    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
               10.*ones(size(c1(1:barstep1:end,2))),...
               zeros(size(c1(1:barstep1:end,2))),'g','g');

    subplot(h2);
    plot(r1(:,2),r1(:,1),'g','Linewidth',1);
    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
               2./100.*r1(1:barstep1r:end,2),...
               zeros(size(r1(1:barstep1r:end,2))),'g','g');

end;


if(iplot == 1)
%    saveas(hf1,plotfile1,'epsc2');
    exportfig(hf1,plotfile1,opts);
end;


return;
