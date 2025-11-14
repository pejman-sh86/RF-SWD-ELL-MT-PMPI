function [] = plot_sinrg_core();

icore = 1;
iplot = 1;
NLGRAD = 100;
hmax = 2.0;
pmin = [1450 1.2 0.0];
pmax = [1550 1.8 1.0];

sample    = 'x_s19_1_4_25_0lay_sample.mat';
mapfile   = 'x_s19_1_4_25_0lay_map.dat';
plotfile1 = 'x_s19_1_4_25_0lay_mapcore.eps';
corefile  = 'core.mat';

map = load(mapfile);
if(icore == 1)
   load(corefile);
   sc_5 = 1.003551;
   sc_6 = 1.010277;
   sr_5 = 1.006025;
   sr_6 = 0.996599;
end;

for i = 1:NLGRAD
   znorm(i) = ((1./NLGRAD)-(1./NLGRAD/2.))+((1./NLGRAD)*(i-1));
   m_rg(((i-1)*4)+1) = map(1)/NLGRAD;
   m_rg(((i-1)*4)+2) = map(2)+(map(3)-map(2))/(NLGRAD-1)*(i-1);
   m_rg(((i-1)*4)+3) = map(4)+sin(znorm(i)*pi/2.)^map(6)...
                       *(map(5)-map(4));
   m_rg(((i-1)*4)+4) = map(7);
end
m_rg(end-2:end) = [map(end-4),map(end-2),map(end)];

NLAY  = NLGRAD;

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure(1);hold on; box on;
set(hf2, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box on;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box on;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box on;

subplot(h1);
plprof(m_rg,hmax,'k',1);
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[1500 1600 1700 1800]);
set(gca,'XTick',[1500 1600 1700 1800]);

subplot(h2);
plprof(m_rg,hmax,'k',2);
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0 2.2]);
set(gca,'XTick',[1.4 1.6 1.8 2.0 2.2]);

subplot(h3);
plprof(m_rg,hmax,'k',3);
set(h3,'layer','top')
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)');
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);

barstep1 = 10;
barstep1r = 10;
barstep2 = 2;
barstep2r = 1;
barstep3 = 10;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'g','Linewidth',1);
    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
               10.*ones(size(c1(1:barstep1:end,2))),...
               zeros(size(c1(1:barstep1:end,2))),'g','g');
%    plot(c2(:,2),c2(:,1),'b','Linewidth',1);
%    errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%               2*10.*ones(size(c2(1:barstep2:end,2))),...
%               zeros(size(c2(1:barstep2:end,2))),'b','b');
%    plot(c3(:,2),c3(:,1),'r','Linewidth',1);
%    errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%               10.*ones(size(c3(1:barstep3:end,2))),...
%               zeros(size(c3(1:barstep3:end,2))),'r','r');


    subplot(h2);
    plot(r1(:,2),r1(:,1),'g','Linewidth',1);
    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
               2./100.*r1(1:barstep1r:end,2),...
               zeros(size(r1(1:barstep1r:end,2))),'g','g');
%    plot(r2(:,2),r2(:,1),'b','Linewidth',1);
%    errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%               4./100.*r2(1:barstep2r:end,2),...
%               zeros(size(r2(1:barstep2r:end,2))),'b','b');
%    plot(r3(:,2),r3(:,1),'r','Linewidth',1);
%    errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%               2./100.*r3(1:barstep3:end,2),...
%               zeros(size(r3(1:barstep3:end,2))),'r','r');
%

end;
if(iplot == 1)
    exportfig(hf2,plotfile1,opts);
end;

return;
