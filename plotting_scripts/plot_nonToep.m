function []=plot_nonToep();
set(0, 'DefaultFigurePaperPosition', [0 0 8 6]);

cov_mat_file = ('x_jan_1_sph_ci2_cov.mat');
plotfile     = ('x_jan_1_sph_ci2_cov.eps');

load(cov_mat_file);
nfreq = length(F(1).freq);
ndat = length(F(1).csave);

nx = 3;
ny = 2;
xim = 0.01;
yim = 0.06;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

gca1 = figure(1);
set(gca1, 'renderer', 'painters')

k = 1;
for ifreq = 1:nfreq
  subplot('Position',[loc(1,ifreq) loc(2,ifreq) spw sph]);
  hold on; box off;

  pcolor(F(ifreq).ang(1:ndat),F(ifreq).ang(1:ndat),F(ifreq).csave);
  shading flat;
  cmap = colormap(flipud(gray));
  colormap(cmap);
  set(gca,'layer','top')
  box on;
  set(gca,'YDir','reverse','FontSize',12);
  set(gca,'YTick',[30 50 70])
  set(gca,'YTickLabel',[]);
  set(gca,'XTick',[30 50 70])
  set(gca,'XTickLabel',[]);
  if ifreq >= 3
    set(gca,'XTick',[30 50 70])
    set(gca,'XTickLabel',[30 50 70]);
  end
  if ifreq >= 3; xlabel('Angle (deg.)','FontSize',14);end;
  if(ifreq == 1 | ifreq == 4)
    set(gca,'YTick',[30 50 70])
    set(gca,'YTickLabel',[30 50 70]);
    ylabel('Angle (deg.)','FontSize',14);
  end

  set(gca,'XLim',[F(ifreq).ang(1) F(ifreq).ang(ndat)],'YLim',[F(ifreq).ang(1) F(ifreq).ang(ndat)],'FontSize',14);
  title(num2str(F(1).freq(ifreq)),'Position',[F(ifreq).ang(1)+5 F(ifreq).ang(1)]);

end

saveas(gca1,plotfile,'eps2')
return;
