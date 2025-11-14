
dex = 0;
set(0, 'DefaultFigurePaperPosition', [0 0 6 3]);

set(gcf, 'renderer', 'painters')
%set(gcf, 'renderer', 'OpenGL')
%set(gcf, 'renderer', 'zbuffer')

cov_mat_file = 'cov_offstat.mat';
res_file = 'residuals_map.mat';
plotfile = 'cov_offstat.eps';

load(cov_mat_file);
load(res_file);
nfreq = 7;
freq = [315 400 500 630 800 1000 1250 1600];

width = 0.38
height = 0.76
llw = [0.1 0.6];
llh = [0.2];

gca1 = figure(1);

for ifreq = 4:4
%for ifreq = nfreq:nfreq
  ang = F(ifreq).ang;
  subplot('Position',[llw(1) llh(1) width height]);
  hold on; box on;
  pcolor(F(ifreq).csave);
  shading flat;
%  colormap('gray');
  set(gca,'YDir','reverse','FontSize',14);
  set(gca,'XTick',[40  80 120],'XTickLabel',{'40'; '80'; '120'})
  set(gca,'YTick',[40  80 120],'YTickLabel',{'40'; '80'; '120'})
  xlabel('Data points','FontSize',14);
  ylabel('Data points','FontSize',14);
  set(gca,'XLim',[1 length(F(ifreq).csave)],'YLim',[1 length(F(ifreq).csave)]);
  text(4.5*length(F(ifreq).csave)/6,length(F(ifreq).csave)/6,...
      [num2str(freq(ifreq)) ' Hz'],'FontSize',10,'BackgroundColor',[1 1 1]);



  subplot('Position',[llw(2) llh(1) width height]);
  hold on; box on;
  pcolor(ang,ang,F(ifreq).csave);
  shading flat;
%  colormap('gray');
  set(gca,'YDir','reverse','FontSize',14);
  set(gca,'XTick',[20 40 60],'XTickLabel',{'20'; '40'; '60'})
  set(gca,'YTick',[20 40 60],'YTickLabel',{'20'; '40'; '60'})

  xlabel('Angle [deg.]','FontSize',14);
  ylabel('Angle [deg.]','FontSize',14);
  set(gca,'XLim',[ang(1) ang(end)],'YLim',[ang(1) ang(end)]);
  text(ang(end)-20,ang(1)+12,...
      [num2str(freq(ifreq)) ' Hz'],'FontSize',10,'BackgroundColor',[1 1 1]);

end

%saveas(gca1,plotfile,'epsc2')
%saveas(gca1,'CAA_plot1.tif','tiffn')
return;
