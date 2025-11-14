function [] = compare_data_rough();

%data1    = 'x_jan_1.mat';
%data1    = 'sim_A_4.mat';
%load(data1);
%F1 = F;
%F2 = load('spher_ref.dat');
F2 = load('spher_ch.dat');
%F3 = load('plane_ref.dat');
F3 = load('plane_ch.dat');

%N = length(F1(1).ang);
%freq = [400 600 900 1250 1600];
freq = [700, 850, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250];

nx = 4
ny = 4
nsubfig = nx*ny;
xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ipar = 1;
for j = 1:length(freq)

    fig1 = figure(1);
    hold on;box on;
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    clear idx;
%    idx = find(F1(1).Rex(1:N,j) == 1);
%    N2 = length(idx);
%    plot(F1(j).ang(idx),F1(j).dat(idx),':k','MarkerSize',5);

    plot(F2(end,:),F2(j,:),'-k','MarkerSize',5);
    plot(F3(end,:),F3(j,:),'-r','MarkerSize',5);
%    plot(F4(end,:),F4(j,:),'-b','MarkerSize',5);
%    plot(F5(end,:),F5(j,:),':b','MarkerSize',5);

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
%legend('1','0.05','0.1','0.2');
legend('spher','plane');


return;
