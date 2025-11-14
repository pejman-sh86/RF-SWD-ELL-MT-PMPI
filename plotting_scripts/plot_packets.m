function [] = plot_packets();
set(0, 'DefaultFigurePaperPosition', [0 0 8 6]);

data1    = 'x_jan_1.mat';
load(data1);
F1 = F;
data2    = 'x_jan_2.mat';
load(data2);
F2 = F;
%data3    = 'x_jan_4.mat';
%load(data3);
%F3 = F;

N = length(F1(1).ang);
freq = [400 600 900 1250 1600];

nx = 3
ny = 3
xim = 0.02;
yim = 0.06;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

ipar = 1;
for j = 1:length(freq)

    fig1 = figure(1);
    hold on;box on;
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    clear idx;
    idx = find(F1(1).Rex(1:N,j) == 1);
    N2 = length(idx);
    plot(F1(j).ang(idx),F1(j).dat(idx),'-k','MarkerSize',5);

    clear idx;
    idx = find(F2(1).Rex(1:N,j) == 1);
    N2 = length(idx);
    plot(F2(j).ang(idx),F2(j).dat(idx),'--k','MarkerSize',5);

%    clear idx;
%    idx = find(F3(1).Rex(1:N,j) == 1);
%    N2 = length(idx);
%    plot(F3(j).ang(idx),F3(j).dat(idx),'.-k','MarkerSize',5);

    set(gca,'FontSize',12);
    text(48,0.9,[num2str(freq(j)) ' Hz'],'FontSize',12)
    set(gca,'XLim',[15 80],'YLim',[0 1]);
    if ((j == 1) | (j == 4))
        set(gca,'YTickLabel',[0 0.3 0.6 0.9],'YTick',[0 0.3 0.6 0.9]);
        ylabel('|V|');
    else
        set(gca,'YTickLabel',[],'YTick',[0 0.3 0.6 0.9]);
    end
    if ((j == 3) | (j == 4) | (j == 5))
        set(gca,'XTickLabel',[30 60],'XTick',[30 60]);
        xlabel('Angle (deg.)');
    else
        set(gca,'XTickLabel',[],'XTick',[30 60]);
    end
%    errorbar(ang(1:step:N),mes(j).dat(1:step:N),L,U,'xk');

end;
saveas(gca,'packets.eps','epsc2');

return;
