%----------------------------------------------------------
% Compute Standard Deviation from PPD
%----------------------------------------------------------
function [] = compute_std();
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);
isave = 1;

%data      = 'x_s21_1_5_40_6lay.mat';
%plotfile1 = 'x_s21_1_5_40_6laycov_data.png';
%plotfile2 = 'x_s21_1_5_40_6laycov_res.png';
%plotfile3 = 'x_s21_1_5_40_6laycov_std.png';
%replica   = 'x_s21_1_5_40_6laycov_rep.dat';
%replica2  = 'x_s21_1_5_40_6laycov_rep.dat';
%sdevs     = 'x_s21_1_5_40_6laycov_sdev.txt';
%residuals = 'x_s21_1_5_40_6laycov_res.txt';
%ibl = 0; % Are we using bottom loss?

%data      = 'x_s20_1_8_40_0lay.mat';
%plotfile1 = 'x_s20_1_8_40_0lay_data.png';
%plotfile2 = 'x_s20_1_8_40_0lay_res.png';
%plotfile3 = 'x_s20_1_8_40_0lay_std.png';
%replica   = 'x_s20_1_8_40_0lay_rep.dat';
%replica2  = 'x_s20_1_8_40_0lay_rep.dat';
%sdevs     = 'x_s20_1_8_40_0lay_sdev.txt';
%residuals = 'x_s20_1_8_40_0lay_res.txt';
%ibl = 1; % Are we using bottom loss?

%data      = 'x_s19_1_8_25_0lay.mat';
%plotfile1 = 'x_s19_1_8_25_0lay_data.png';
%plotfile2 = 'x_s19_1_8_25_0lay_res.png';
%plotfile3 = 'x_s19_1_8_25_0lay_std.png';
%replica   = 'x_s19_1_8_25_0lay_rep.dat';
%replica2  = 'x_s19_1_8_25_0lay_rep.dat';
%sdevs     = 'x_s19_1_8_25_0lay_sdev.txt';
%residuals = 'x_s19_1_8_25_0lay_res.txt';
%ibl = 1; % Are we using bottom loss?

%data      = 'x_s16_1_5_25_7lay.mat';
%plotfile1 = 'x_s16_1_5_25_7laycov_data.png';
%plotfile2 = 'x_s16_1_5_25_7laycov_res.png';
%plotfile3 = 'x_s16_1_5_25_7laycov_std.png';
%replica   = 'x_s16_1_5_25_7lay_rep.dat';
%replica2  = 'x_s16_1_5_25_7lay_rep.dat';
%sdevs     = 'x_s16_1_5_25_7laycov_sdev.txt';
%residuals = 'x_s16_1_5_25_7laycov_res.txt';
%ibl = 0; % Are we using bottom loss?

%data      = 'x_s13_2_10_50_9lay.mat';
%plotfile1 = 'x_s13_2_10_50_9lay_data.png';
%plotfile2 = 'x_s13_2_10_50_9lay_res.png';
%plotfile3 = 'x_s13_2_10_50_9lay_std.png';
%replica   = 'x_s13_2_10_50_9lay_rep.dat';
%replica2  = 'x_s13_2_10_50_9lay_rep.dat';
%sdevs     = 'x_s13_2_10_50_9lay_sdev.txt';
%residuals = 'x_s13_2_10_50_9lay_res.txt';
%ibl = 0; % Are we using bottom loss?

data      = 'x_s07_1_1_30_nodisp.mat';
plotfile1 = 'x_s07_1_1_30_nodisp_data.png';
plotfile2 = 'x_s07_1_1_30_nodisp_res.png';
plotfile3 = 'x_s07_1_1_30_nodisp_std.png';
replica   = 'x_s07_1_1_30_nodisp_rep.dat';
replica2  = 'x_s07_1_1_30_nodisp_rep.dat';
sdevs     = 'x_s07_1_1_30_nodisp_sdev.txt';
residuals = 'x_s07_1_1_30_nodisp_res.txt';
ibl = 0; % Are we using bottom loss?

%data      = 'x_s05_1_8_32_0lay.mat';
%plotfile1 = 'x_s05_1_8_32_0lay_data.png';
%plotfile2 = 'x_s05_1_8_32_0lay_res.png';
%plotfile3 = 'x_s05_1_8_32_0lay_std.png';
%replica   = 'x_s05_1_8_32_0lay_rep.dat';
%replica2  = 'x_s05_1_8_32_0lay_rep.dat';
%sdevs     = 'x_s05_1_8_32_0lay_sdev.txt';
%residuals = 'x_s05_1_8_32_0lay_res.txt';
%ibl = 1; % Are we using bottom loss?

%data      = 'x_s05_2_8_25_4lay.mat';
%plotfile1 = 'x_s05_2_8_25_4laycov_data.png';
%plotfile2 = 'x_s05_2_8_25_4laycov_res.png';
%plotfile3 = 'x_s05_2_8_25_4laycov_std.png';
%replica   = 'x_s05_2_8_25_4laycov_rep.dat';
%replica2  = 'x_s05_2_8_25_4laycov_rep.dat';
%sdevs     = 'x_s05_2_8_25_4laycov_sdev.txt';
%residuals = 'x_s05_2_8_25_4laycov_res.txt';
%ibl = 0; % Are we using bottom loss?

%data      = 'x_s04_1_3_32_0layb.mat';
%plotfile1 = 'x_s04_1_3_32_0layb_data.png';
%plotfile2 = 'x_s04_1_3_32_0layb_res.png';
%plotfile3 = 'x_s04_1_3_32_0layb_std.png';
%replica   = 'x_s04_1_3_32_0layb_rep.dat';
%replica2  = 'x_s04_1_3_32_0layb_rep.dat';
%sdevs     = 'x_s04_1_3_32_0layb_sdev.txt';
%residuals = 'x_s04_1_3_32_0layb_res.txt';
%ibl = 1; % Are we using bottom loss?

%data      = 'x_s02_1_3_25_7lay.mat';
%plotfile1 = 'x_s02_1_3_25_7laycov_data.png';
%plotfile2 = 'x_s02_1_3_25_7laycov_res.png';
%plotfile3 = 'x_s02_1_3_25_7laycov_std.png';
%replica   = 'x_s02_1_3_25_7laycov_rep.dat';
%replica2  = 'x_s02_1_3_25_7laycov_rep.dat';
%sdevs     = 'x_s02_1_3_25_7laycov_sdev.txt';
%residuals = 'x_s02_1_3_25_7laycov_res.txt';
%ibl = 0; % Are we using bottom loss?

%data      = 'x_s01_1_3_20_6lay.mat';
%plotfile1 = 'x_s01_1_3_20_6lay_data.png';
%plotfile2 = 'x_s01_1_3_20_6lay_res.png';
%plotfile3 = 'x_s01_1_3_20_6lay_std.png';
%replica   = 'x_s01_1_3_20_6lay_rep.dat';
%replica2  = 'x_s01_1_3_20_6lay_rep.dat';
%sdevs     = 'x_s01_1_3_20_6lay_sdev.txt';
%residuals = 'x_s01_1_3_20_6lay_res.txt';
%ibl = 0; % Are we using bottom loss?

%data = 'sim_A_4.mat';
%plotfile1 = 'sim_A_4_data.png';
%plotfile2 = 'sim_A_4_res.png';
%plotfile3 = 'sim_A_4_std.png';
%replica = 'sim_A_4_rep.txt';
%replica2 = 'sim_A_4_rep.txt';
%sdevs = 'sim_A_4_sdev.txt';
%residuals = 'sim_A_4_res.txt';
%ibl = 0; % Are we using bottom loss?

NAVE = 30;

load(data);
mes = F;
rep = load(replica);
rep2 = load(replica2);
%tru = load(true);
 
freq = F(1).freq(:);
M = length(F(1).minlim);
N = size(rep,2);
Nf = length(freq);
disp('    M,   N,   Nf,ceil(M/Nf): '),fprintf(1,'%5i',M,N,Nf,ceil(M/Nf));

sd2 = zeros(1,length(freq));
step = 1;
sd = zeros(length(freq),N);

nx = 5;
ny = 5;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

for j = 1:length(freq)
    ang = mes(1).ang(1:N);
    res(j,:) = mes(j).dat(1:N)'-rep(j,1:N);
    res(j,:) = res(j,:).*F(1).Rex(:,j)';
 
    clear idx;
    clear idx2;
    idx  = find(F(1).Rex(1:N,j) == 1);
    idx2 = find(F(1).Rex(1:N,j) == 0);
    N2 = length(idx);
    sd2(j) = sqrt(1/(N2-((M+Nf)/Nf))*sum((mes(j).dat(idx)'-rep(j,idx)).^2));

    resabs = abs(res);
    resabs(j,idx2) = mean(resabs(j,:));
    for i=1:N

        if(i <= NAVE)
            NAVE1 = 1;
        else
            NAVE1 = i-NAVE;
        end
        if(i >= N-NAVE)
            NAVE2 = N;
        else
            NAVE2 = i+NAVE;
        end
        sd(j,i) = sqrt(mean(resabs(j,NAVE1:NAVE2).^2));

    end

    fig1 = figure(1);
%    subplot('Position',[loc(1,2) loc(2,2) spw sph]);
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    plot(ang,sd(j,:),'k');
    set(gca,'XLim',[ang(1) ang(end)])
    set(gca,'FontSize',14)
    set(gca,'LineWidth',1);
    if(ibl == 0)
       text(48,0.3,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[0 0.3]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[0.0 0.1 0.2],'YTick',[0.0 0.1 0.2]);
           ylabel('Residuals');
       else
           set(gca,'YTickLabel',[],'YTick',[0.0 0.1 0.2]);
       end
    else
       text(48,3.5,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[0. 4.]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[0 1 2 3 4],'YTick',[0 1 2 3 4]);
           ylabel('Std. dev.');
       else
           set(gca,'YTickLabel',[],'YTick',[0 1 2 3 4]);
       end
    end
    if (j > length(freq)-nx) 
        set(gca,'XTickLabel',[30 60],'XTick',[30 60]);
        xlabel('Angle (deg.)');
    else
        set(gca,'XTickLabel',[],'XTick',[30 60]);
    end

    fig2 = figure(2);
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    plot(ang(idx),res(j,idx),'k');
    plot([0 90],[0 0],'--k');
    set(gca,'XLim',[ang(1) ang(end)])
    ylims = get(gca,'YLim');
    if(min(res(j,:))>0) set(gca,'YLim',[0 ylims(2)]);end;
    if(max(res(j,:))<0) set(gca,'YLim',[ylims(1) 0]);end;
    set(gca,'FontSize',14)
    set(gca,'LineWidth',1);
    if(ibl == 0)
       text(48,0.2,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[-0.3 0.3]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[-0.2 -0.1 0.0 0.1 0.2],'YTick',[-0.2 -0.1 0.0 0.1 0.2]);
           ylabel('Residuals');
       else
           set(gca,'YTickLabel',[],'YTick',[-0.2 -0.1 0.0 0.1 0.2]);
       end
    else
       text(48,3.,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[-4. 4.]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[-3 -2 -1 0 1 2 3],'YTick',[-3 -2 -1 0 1 2 3]);
           ylabel('Residuals (dB)');
       else
           set(gca,'YTickLabel',[],'YTick',[-3 -2 -1 0 1 2 3]);
       end
    end
    if (j > length(freq)-nx) 
        set(gca,'XTickLabel',[30 60],'XTick',[30 60]);
        xlabel('Angle (deg.)');
    else
        set(gca,'XTickLabel',[],'XTick',[30 60]);
    end

    fig3 = figure(3);
    subplot('Position',[loc(1,j) loc(2,j) spw sph]);
    hold on;box on;
    set(gca,'FontSize',12);
    plot(ang(idx),mes(j).dat(idx),'.:k','MarkerSize',5);
    plot(ang(:),rep(j,:),'-k','LineWidth',1);
    if(ibl == 0)
       text(48,0.9,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[0 1]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[0 0.3 0.6 0.9],'YTick',[0 0.3 0.6 0.9]);
           ylabel('Refl. Coeff.');
       else
           set(gca,'YTickLabel',[],'YTick',[0 0.3 0.6 0.9]);
       end
    else
       text(48,15,[num2str(freq(j)) ' Hz'],'FontSize',12)
       set(gca,'XLim',[5 80],'YLim',[10 40]);
       if ((j == 1) | (j == nx+1) | (j == (2*nx)+1)) 
           set(gca,'YTickLabel',[15 25 35],'YTick',[15 25 35]);
           ylabel('Bottom Loss (dB)');
       else
           set(gca,'YTickLabel',[],'YTick',[15 25 35]);
       end
    end
    if (j > length(freq)-nx) 
        set(gca,'XTickLabel',[30 60],'XTick',[30 60]);
        xlabel('Angle (deg.)');
    else
        set(gca,'XTickLabel',[],'XTick',[30 60]);
    end
    L = sd2(j);
    U = sd2(j);
    disp(sd2(j)),disp(length(idx));
%    errorbar(ang(1:step:N),mes(j).dat(1:step:N),L,U,'xk');

end
res(j+1,:) = rep(j,1:N)';
disp('mean sd2:'),disp(mean(sd2));
disp('mean sd:'),disp(mean(sd'));

if(isave == 1)
%  saveas(fig3,plotfile1,'epsc2');
%  saveas(fig2,plotfile2,'epsc2');
%  saveas(fig1,plotfile3,'epsc2');
  saveas(fig3,plotfile1,'png');
  saveas(fig2,plotfile2,'png');
  saveas(fig1,plotfile3,'png');

  save(sdevs,'sd','-ascii');
  save(residuals,'res','-ascii');
end

return;
