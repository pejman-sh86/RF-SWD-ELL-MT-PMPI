%
% Plot 2-D Histograms from FGS samples
%
function plot_2dhist;

set(0, 'DefaultFigurePaperPosition', [0 0 12 7]);
set(gcf, 'renderer', 'painters')
%set(gcf, 'renderer', 'OpenGL')
%set(gcf, 'renderer', 'zbuffer')

sample = 'picked_tt_int_sample.mat';
plotfile1 = strrep(sample,'.mat','joints.eps');
load(sample);

ip1 = [ 1     4     4    11    15];
ip2 = [14    15    16    18    16];
%ip1 = [ 1     3     5     6     9    10];
%ip2 = [ 3     4     7     8    11    12];

minlim =  [1.0, 1450.,1.5, 1500.,...
           0.5, 1500.,3.0, 1500.,...
           1.0, 1500.,1.0, 1500.,...
           0.0, 0.0, 0.0, 0.0,...
           0.0, 0.0, 0.0];
maxlim =  [3., 1600.,4., 1800.,...
           4., 1800.,10., 1800.,...
           5., 1800.,10., 1800.,...
           0.5, 1.5, 1.5, 1.5,...
           1.5, 1.5, 0.5];

xlabels = {'h1','c1','h2','c2','h3','c3','h4','c4',...
           'h5','c5','h6','c6','a1','a2','a3','a4','a5','a6','a7'};

nplot = length(ip1);
nrow = ceil(nplot/3)

m = A(:,2:end);
for i = 1:7
  m(:,12+i) = m(:,12+i)-maxlim(12+i)/2;
end

npar = length(m(1,:));
ndat = length(m(:,1));

nx = 3
ny = 2
xim = 0.08;
yim = 0.08;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

figure(1);
for i = 1:nplot
  [n,x,nbins] = histmulti5([m(:,ip2(i)) m(:,ip1(i))],[100 100]);

%set(gca,'ydir','reverse');
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);
  hold on;box on;
%  colormap(gray);
%  pcolor(x(:,2),x(:,1),n),shading interp;
  pcolor(x(:,2),x(:,1),n);
  shading flat;
  xlabel(xlabels(ip1(i)),'FontSize',14);
  ylabel(xlabels(ip2(i)),'FontSize',14);
  set(gca,'FontSize',14)
  set(gca,'XLim',[min(x(:,2)) max(x(:,2))]);
  set(gca,'YLim',[min(x(:,1)) max(x(:,1))]);
  if(max(x(:,2))>1000)
      set(gca,'XTick',[1520 1570 1620 1670]);
  else
%      set(gca,'XTick',[1 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.0 5.4 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0]);
%      set(gca,'XTick',[1 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.0 5.4 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0]);
  end
  if(max(x(:,1))>1000)
      set(gca,'YTick',[1520 1570 1620 1670]);
  else
%      set(gca,'YTick',[1 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.0 5.4 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0]);
%      set(gca,'YTick',[1 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 4.6 5.0 5.4 5.8 6.2 6.6 7.0 7.4 7.8 8.2 8.6 9.0]);
  end
end

saveas(gca,plotfile1,'epsc2');

return;
