function [] = mat_plotc(A,r);
set(0, 'DefaultFigurePaperPosition', [0 0 5 5]);

%pcolor(newang,newang,F(6).csave)
pcolor(r,r,A)
shading flat
%colorbar
set(gca,'YDir','reverse')
set(gca,'FontSize',14)
%xlabel('Angle (deg.)')
%ylabel('Aangle (deg.)')
%    set(gca,'XTick',[20 40 60 80])
%    set(gca,'XTickLabel',[20 40 60 80])
%    set(gca,'YTick',[20 40 60 80])
%    set(gca,'YTickLabel',[20 40 60 80])
    xlabel('Range (m)')
    ylabel('Rrange (m)')
    set(gca,'XLim',[r(1) r(end)])
    set(gca,'XTick',[100 200 300 400])
    set(gca,'XTickLabel',[100 200 300 400])
    set(gca,'YLim',[r(1) r(end)])
    set(gca,'YTick',[100 200 300 400])
    set(gca,'YTickLabel',[100 200 300 400])

%saveas(gca,'matrix.eps','epsc2')
%saveas(gca,'matrix.png','png')

