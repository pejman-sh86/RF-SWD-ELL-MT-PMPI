%----------------------------------------------------------
% Compute Standard Deviation from PPD
%----------------------------------------------------------
function [] = compute_std(filebase);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);
isave = 1;

%filebase  = 'elba_0630_2000';
data      = strcat(filebase,'.mat');
plotfile1 = strcat(filebase,'_data.png');
plotfile2 = strcat(filebase,'_res.png');
plotfile3 = strcat(filebase,'_std.png');
replica   = strcat(filebase,'_rep.dat');
replica2  = strcat(filebase,'_rep.dat');
sdevs     = strcat(filebase,'_sdev.txt');
residuals = strcat(filebase,'_res.txt');
ibl = 1; % Are we using bottom loss?

NAVE = 20;

load(data);
mes = F;
rep = dlmread(replica);
rep2 = dlmread(replica2);
 
freq = mes(1).freq;
NBAND = length(freq);
M = 7;

step = 1;

nx = 5;
ny = 3;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

for j = 1:NBAND

   N = sum(mes(1).Rex(j,:));
   F(j).ang = mes(j).ang(find(mes(1).Rex(j,:)==1));
   F(j).dat = mes(j).dat(find(mes(1).Rex(j,:)==1));
   F(j).rep = rep(j,find(mes(1).Rex(j,:)==1));
   fprintf(1,'    M,   N,   NBAND,ceil(M/NBAND): \n');
   fprintf(1,'%6i',M,N,NBAND,ceil(M/NBAND));
   fprintf(1,'\n');

   F(j).res = F(j).dat-F(j).rep;
 
   F(j).sd2 = sqrt(1./(N-((M+NBAND)/NBAND))*sum((F(j).res).^2.));

   F(j).resabs = abs(F(j).res);
   for i=1:N;

       if(i <= NAVE/2)
           NAVE1 = 1;
       else
           NAVE1 = i-NAVE/2;
       end
       if(i >= N-NAVE/2)
           NAVE2 = N;
       else
           NAVE2 = i+NAVE/2;
       end
       %disp([NAVE1,NAVE2])
       F(j).sd(i) = sqrt(mean(F(j).resabs(NAVE1:NAVE2).^2.));

   end

   fig1 = figure(1);
   subplot('Position',[loc(1,j) loc(2,j) spw sph]);
   hold on;box on;
   plot(F(j).ang,F(j).sd,'k');
   set(gca,'XLim',[F(j).ang(1) F(j).ang(end)])
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
    plot(F(j).ang,F(j).res,'k');
    plot([0 90],[0 0],'--k');
    set(gca,'XLim',[F(j).ang(1) F(j).ang(end)])
    ylims = get(gca,'YLim');
    if(min(F(j).res)>0) set(gca,'YLim',[0 ylims(2)]);end;
    if(max(F(j).res)<0) set(gca,'YLim',[ylims(1) 0]);end;
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
    plot(F(j).ang,F(j).dat,'.:k','MarkerSize',5);
    plot(F(j).ang,F(j).rep,'-k','LineWidth',1);
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
       text(48,35,[num2str(freq(j)) ' Hz'],'FontSize',12)
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
    L = F(j).sd2;
    U = F(j).sd2;
%    errorbar(ang(1:step:N),mes(j).dat(1:step:N),L,U,'xk');

end

%res(j+1,:) = rep(j,1:N)';
%disp('mean sd2:'),disp(mean(sd2));
%disp('mean sd:'),disp(mean(sd'));

if(isave == 1)
%  saveas(fig3,plotfile1,'epsc2');
%  saveas(fig2,plotfile2,'epsc2');
%  saveas(fig1,plotfile3,'epsc2');
  saveas(fig3,plotfile1,'png');
  saveas(fig2,plotfile2,'png');
  saveas(fig1,plotfile3,'png');

  for i = 1:NBAND
    tmp(i) = length(F(i).res);
  end
  save(sdevs,'tmp','-ascii');
  %tmp=F(1).minlim;save(sdevs,'tmp','-ascii','-append');
  %tmp=F(1).maxlim;save(sdevs,'tmp','-ascii','-append');
  %tmp=F(1).mtrue;save(sdevs,'tmp','-ascii','-append');
  for i = 1:NBAND
    tmp = F(i).sd;
    save(sdevs,'tmp','-ascii','-append');
  end
  for i = 1:NBAND
    tmp = F(i).ang;
    save(sdevs,'tmp','-ascii','-append');
  end
  clear tmp;

  for i = 1:NBAND
    tmp(i) = length(F(i).res);
  end
  save(residuals,'tmp','-ascii');
  %tmp=F(1).minlim;save(residuals,'tmp','-ascii','-append');
  %tmp=F(1).maxlim;save(residuals,'tmp','-ascii','-append');
  %tmp=F(1).mtrue;save(residuals,'tmp','-ascii','-append');
  for i = 1:NBAND
    tmp = F(i).res;
    save(residuals,'tmp','-ascii','-append');
  end
  for i = 1:NBAND
    tmp = F(i).ang;
    save(residuals,'tmp','-ascii','-append');
  end

end

return;
