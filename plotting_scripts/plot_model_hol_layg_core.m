function [] = plot_model_hol_layg();

set(0, 'DefaultFigurePaperPosition', [0 0 11 5]);
%set(0, 'DefaultFigurePaperPosition', [0 0 14 8]);

corefile  = 'core15.mat';
plotfile1 = 'sim_C_trueprofile.eps';

icore = 0;

if(icore == 1)
    load(corefile);
end;

M(1).m = [0.30,  1480, 1.30, 1.6, 0.01, ...
          0.35,  1510, 1.65,      0.30, ...
          1.20,  1500, 1.55,      0.30, ...
          0.40,  1700, 1.90,      0.30, ...
                 1560, 1.70,      0.01 ];
%M(4).m = [];
%M(5).m = [];
         
npl = 4;

col = {'-k' '-b' '-r' '--k' '--r' ':k'};
col = char(col);

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.04 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',6,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf1=figure(1);hold on; box on;
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
set(gca,'Fontsize',14,'XLim',[1450 1750]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
hold on; box on;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
set(gca,'Fontsize',14,'XLim',[1.2 2.1]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.2 1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.2 1.4 1.6 1.8 2.0]);
hold on; box on;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
set(gca,'Fontsize',14,'XLim',[0.0 1.0]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[0.0 0.25 0.5 0.75 1.0]);
set(gca,'XTick',[0.0 0.25 0.5 0.75 1.0]);
hold on; box on;

hmax = 3;

for i = 1:size(M,2)

  nl = (length(M(i).m)-4)/npl
  mtmp = [M(i).m(1:3) M(i).m(5:end)];
% Velocity
  subplot(h1);box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);
  plprof(mtmp,hmax,col(i,:),1);
  xlabel('Velocity (m/s)');
  ylabel('Depth (m)');
  ylims = get(gca,'Ylim');

% Density
  subplot(h2);box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);
  plprofg(M(i).m,hmax,col(i,:),2);
  xlabel('Density (g/ccm)');

% Attenuation
  subplot(h3);hold on; box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);
  plprof(mtmp,hmax,col(i,:),3);
  xlabel('Attenuation (dB/L)');

end;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'r','Linewidth',2);
%    plot(c2(:,2),c2(:,1),'r','Linewidth',2);

    subplot(h2);
    plot(r1(:,2),r1(:,1),'r','Linewidth',2);
%    plot(r2(:,2),r2(:,1),'r','Linewidth',2);
end;

subplot(h1); box on;
subplot(h2); box on;
subplot(h3); box on;
%saveas(hf1,plotfile1,'epsc2');
exportfig(hf1,plotfile1,opts);

%subplot(h1);hold on; box on;
%legend([handle(1) handle(2) handle(3)],'true model','sph. inv.','pl. inv.')

return;
