function [] = compute_std_tt();

set(0, 'DefaultFigurePaperPosition', [0 0 12 8]);
i_save = 1;
factor = 2;
trace = 'x.mat'
%trace = 'pro_sim_A.mat'
sdfile = 'sigma_map.txt'
residuals = 'residuals_map.mat';
plotfile1 = 'sigma_map.eps';
plotfile2 = 'residuals_map.eps';
replica = 'replica_map.mat';
load(replica);
rep = rep2;
picked = 'picked_tt_int.mat';
%picked = 'picked_tt2.mat';
load(picked);
%res2 = (rep-tt(:,1:70))*1000;
%averes2 = mean(res2')
load(trace);

offsets = [0.70195118      1.2069232      1.0393571     0.51960475];
offsets = offsets - [1.0, 1.5, 1.5, 1.5]./2;
%offsets = [0.1888662,     0.34376936,     0.27528262,...
%           0.38473811, 0.79711843,     0.32147155,...      
%           0.2211022];
%offsets = offsets - [0.3, 1.0, 1.0, 1.0,1.0, 1.0, 0.5]./2;

NAVE = 10;
NRAN = length(rep)
%r = x.r(1:NRAN) * -1;
NLAY = size(rep,1);
fact = [1, 1, 1, 1, 1, 1, 1];
res = rep - tt(:,1:NRAN);

nx = 2
ny = 4
xim = 0.08;
yim = 0.04;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

sd = zeros(size(rep));
for j=1:NLAY
%    for i=1+NAVE:NRAN-NAVE
    for i=1:NRAN
        
        if(i <= NAVE)
            NAVE1 = 1;
        else
            NAVE1 = i-NAVE;
        end
        if(i >= NRAN-NAVE)
            NAVE2 = NRAN;
        else
            NAVE2 = i+NAVE;
        end
        sd(j,i) = factor*sqrt(mean(res(j,NAVE1:NAVE2).^2));
%        sd(j,i) = sqrt(mean(res(j,i-NAVE:i+NAVE).^2));

    end
%    sd(j,1:NAVE) = sd(j,NAVE+1);
%    sd(j,end-NAVE+1:end) = sd(j,end-NAVE);
    sd(j,:) = sd(j,:)*fact(j);
%    subplot(4,2,j)
    fig1 = figure(1);
%    subplot('Position',[loc(1,2) loc(2,2) spw sph]);
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    plot(r,sd(j,:)*1000,'k');
    plot([0 r(end)],[0 0],'--k');
    set(gca,'XLim',[20 r(end)])
    ylims = get(gca,'YLim');
    set(gca,'FontSize',16)
    set(gca,'XTickLabel',[])
%    set(gca,'YTick',[0 0.1 0.2])
%    set(gca,'YTickLabel',[0.0 0.1 0.2]);
    set(gca,'XTick',[100 200 300 400])
    set(gca,'XTickLabel',[]);
    set(gca,'LineWidth',1);
%    xlabel('Range (m)');
%    ylabel('Std. dev. (ms)');
    if(j == 7 | j == 6) set(gca,'XTickLabel',[100 200 300 400]); end;
    if(j == 7 | j == 6) xlabel('Range (m)'); end;
    if(j==3) ylabel('Std. dev. (ms)'); end;

    fig2 = figure(2);
%    subplot('Position',[loc(1,1) loc(2,1) spw sph]);
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    plot(r,res(j,:)*1000,'k');
    plot([0 r(end)],[0 0],'--k');
    plot([0 r(end)],[offsets(j) offsets(j)],':k','LineWidth',1);
%    plot([0 70],[averes2(j) averes2(j)],'--k');
    set(gca,'XLim',[20 r(end)])
    ylims = get(gca,'YLim');
    if(min(res(j,:))>0) set(gca,'YLim',[0 ylims(2)]);end;
    if(max(res(j,:))<0) set(gca,'YLim',[ylims(1) 0]);end;
    set(gca,'FontSize',16)
    set(gca,'YLim',[-.2 .2])
%    set(gca,'YTick',[-0.2 0 0.2])
%    set(gca,'YTickLabel',[-0.2 0.0 0.2]);
    set(gca,'XTickLabel',[])
    set(gca,'XTick',[100 200 300 400])
    set(gca,'LineWidth',1);
%    ylabel('Error (ms)');
    if(j == 7 | j == 6) set(gca,'XTickLabel',[100 200 300 400]); end;
    if(j == 7 | j == 6) xlabel('Range (m)'); end;
    if(j == 3) ylabel('Error (ms)'); end;
%    if(j == 3) ylabel('Residuals (ms)'); end;
    mis1(j) = sum(((res(j,:))./sd(j,:)).^2);
    mis2(j) = sum(((res(j,:)-(offsets(j)/1000.))./sd(j,:)).^2);
    mis3(j) = mean(res(j,:));
    mis4(j) = mean(res(j,:)-(offsets(j)/1000.));
end

mis1
mis2
mis3*1000
mis4*1000
sum(mis1)/7
sum(mis2)/7
(sum(abs(mis3))/7)*1000
(sum(abs(mis4))/7)*1000

if(i_save == 1)
    saveas(fig1,plotfile1,'epsc2');
    saveas(fig2,plotfile2,'epsc2');
    save(sdfile,'sd','-ascii');
    save(residuals,'res');
end

return;
