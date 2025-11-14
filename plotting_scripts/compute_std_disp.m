%----------------------------------------------------------
% Compute Standard Deviation from PPD
%----------------------------------------------------------
function [] = compute_std_disp();
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);
isave = 1;

data      = '2_b_80m_y_m0_2lay.txt';
plotfile1 = '2_b_80m_y_m0_2lay_data.png';
plotfile2 = '2_b_80m_y_m0_2lay_res.png';
plotfile3 = '2_b_80m_y_m0_2lay_std.png';
replica   = '2_b_80m_y_m0_2lay_rep.dat';
replica2  = '2_b_80m_y_m0_2lay_rep.dat';
sdevs     = '2_b_80m_y_m0_2lay_sdev.txt';
residuals = '2_b_80m_y_m0_2lay_res.txt';
ibl = 1; % Are we using bottom loss?

NAVE = 60;

mes  = load(data);
rep  = load(replica);
rep2 = load(replica2);
 
freq = mes(2:end,1);
N = size(rep,2);

step = 1;
sd = zeros(1,N);

res(j,:) = mes(j).dat(1:N)'-rep(j,1:N);
res(j,:) = res(j,:).*F(1).Rex(:,j)';
 
sd2(j) = sdev()
		
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
       set(gca,'XLim',[5 80],'YLim',[0 20]);
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
