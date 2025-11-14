function [] = rep_vs_data_hol();
%
%PLOT reflection model and data 
%
paper = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 9 6]);

ml_file = 'ml_std_1.mat';
plotfile(1,:) = 'rep_vs_data_0a.eps';
plotfile(2,:) = 'rep_vs_data_0b.eps';

load('sim_A_1dB_4_mapdat.mat');
Frep = F;
load('sim_A_1dB_4_truedat.mat');
Frep2 = F;
load('sim_A_sph_4.mat');
Fsph = F;
load('sim_A_1dB_4.mat');
nfreq = length(F(1).freq);

nsubpfig = 9;
nfig = ceil(nfreq/nsubpfig);

if(paper == 1)
  llw_tmp = [0.1 0.4 0.7];
  llw = [llw_tmp llw_tmp llw_tmp];
  llh = [0.7 0.4 0.1];
else
  llw_tmp = [0.15 0.43 0.71];
  llw = [llw_tmp llw_tmp llw_tmp];
  llh = [0.69 0.42 0.15];
end
spw = 0.28
sph = 0.28


kdat = 1;
for ifig = 1:nfig

hand(ifig) = figure(ifig);hold on;
i = 1;

for k=1:nsubpfig
   if k <= 3
     subplot('Position',[llw(i) llh(1) spw sph]);
     set(gca,'XTickLabel',[]);
     i = i + 1;
   elseif k <= 6
     subplot('Position',[llw(i) llh(2) spw sph]);
     set(gca,'XTickLabel',[]);
     i = i + 1;
   else
     subplot('Position',[llw(i) llh(3) spw sph]);
     set(gca,'XTickLabel',[20 40 60 80]);
     xlabel('Angle (deg.)','FontSize',14);
     i = i + 1;
   end
   if ((k == 1)|(k == 4)|(k == 7))
       set(gca,'YTickLabel',[0 5 10 15]);
       ylabel('BL (dB)','FontSize',14);
   else
       set(gca,'YTickLabel',[]);
   end
   set(gca,'XTick',[20 40 60 80],'FontSize',14);
   set(gca,'YTick',[0 5 10 15]);
   hold on;box on
   Err = ones(length(F(k).ang(1:3:end)),1);
   plot(F(kdat).ang, F(kdat).dat,'-k','LineWidth',1);
%   dashline(F(kdat).ang, F(kdat).dat,'.k',60,'.k',60);
   errorbar(F(k).ang(1:3:end),F(k).dat(1:3:end),Err,'k.');
   hold on;
   
   plot(Frep(kdat).ang,Frep(kdat).dat,'b');
   plot(Frep2(kdat).ang,Frep2(kdat).dat,'r');
%     plot(Fsph(kdat).ang,Fsph(kdat).dat,':k');
%   title([num2str(freq(k)) ' Hz'] )
   text(65,18,[num2str(F(1).freq(kdat)) ' Hz'])
   axis([10 85 0 20]);
%   legend('simA','pl map','pl true','sp true','Location','SouthEast'); 

  kdat = kdat+1;
  if kdat>nfreq;break;end;
end
%saveas(hand(ifig),plotfile(ifig,:),'epsc2');
end

return
% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
