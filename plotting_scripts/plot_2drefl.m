function []=plot_2drefl();
set(gcf, 'renderer', 'painters')

load('x_s13_1_broad.mat');

%r = sqrt(150^2./(sin(F(1).ang/90*pi/2)).^2-150^2);
r = F(1).r;
ang = F(1).ang;
freq = F(1).freq;
r = fliplr(r);

for i=1:length(F(1).freq)

    R(i,:) = F(i).dat;

end

size(ang)
size(freq)
%R = abs(R)
hf1=figure(1);
pcolor(ang,freq,R);shading flat;
set(gca,'layer','top')
box on;
xlabel('Angle (deg.)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
text(92,1100,'|R|','Rotation',90,'FontSize',12)
caxis([0 1.])

h = colorbar;
set(gca,'FontSize',14);
set(h,'FontSize',14);

return;

hf2 = figure(2);hold on; box on;
pcolor(r,freq,R),shading flat;
set(gca,'layer','top')
box on;
xlabel('Range (m)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
set(gca,'XLim',[20 560]);
set(gca,'YLim',[380 1680]);
caxis([0 0.8])
plot([300 300],[800 1200],'k','LineWidth',2);
plot([240 360],[1000 1000],'k','LineWidth',2);

h = colorbar;
set(gca,'FontSize',14);
set(h,'FontSize',14);
text(690,1100,'R','Rotation',90,'FontSize',12)
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',7,'Height',5,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
exportfig(hf2,'plot.eps',opts);
return;
