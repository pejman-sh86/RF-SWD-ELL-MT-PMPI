function [] = rep_vs_data_hol();
%
%PLOT reflection model and data 
%
paper = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 9 7]);

ml_file = 'site20b_ml_std_1.mat';
plotfile(1,:) = 'site20b_rep_vs_data_0a.eps';
plotfile(2,:) = 'site20b_rep_vs_data_0b.eps';

load('site20b_rep0.mat');
%load('sim_A_1dB_4_mapdat.mat');
Frep = F;
load('hol_site20b.mat');
%load('sim_A_1dB_4.mat');
nfreq = length(F(1).freq);

nsubpfig = 9;
nfig = ceil(nfreq/nsubpfig);

if(paper == 1)
  llw_tmp = [0.15 0.43 0.71];
  llw = [llw_tmp llw_tmp llw_tmp];
%  llw = [0.10 0.335 0.57 0.805 0.10 0.335 0.57 0.805 0.10 0.335 0.57 0.805];
  llh = [0.69 0.42 0.15];
else
  llw_tmp = [0.15 0.43 0.71];
  llw = [llw_tmp llw_tmp llw_tmp];
  llh = [0.69 0.42 0.15];
end
spw = 0.195
sph = 0.26


kdat = 1;
for ifig = 1:nfig

hand(ifig) = figure(ifig);hold on;
i = 1;
for k=1:nsubpfig
   if k <= 3
     subplot('Position',[llw(i) llh(1) spw sph]);
     i = i + 1;
   elseif k <= 6
     subplot('Position',[llw(i) llh(2) spw sph]);
     i = i + 1;
   else
     subplot('Position',[llw(i) llh(3) spw sph]);
     i = i + 1;
   end
   hold on;box on
%   Err = ml_sd(k)*ones(nang(k),1);
   plot(F(kdat).ang, F(kdat).dat,'--k.');
%   errorbar(F(k).ang,F(k).dat,Err,'k.');
   hold on;
   
   if(paper == 1)
     plot(Frep(kdat).ang,Frep(kdat).dat,'b');
   else
     plot(Frep(kdat).ang,Frep(kdat).dat,'b');
   end
%   title([num2str(freq(k)) ' Hz'] )
   text(55,40,[num2str(F(1).freq(kdat)) ' Hz'])
   axis([0 90 10 45]);
%   legend('hol','rep'); 
%   if(paper == 1)
%     set(gca,'XTick',[0 30 60 90],'YTick',[0 1],...
%     'YLim',[0 1.1],'FontSize',12);
%   else
%     set(gca,'XTick',[0 40 80],'YTick',[0 1],...
%     'YLim',[0 1.1],'FontSize',16);
%   end
%   if k >=  6;xlabel('Angle [deg]');end
%   if k == 1;ylabel('BL [dB]');end
%   if(k < 6) 
%     set(gca,'XTickLabel',[]);
%   else 
%     if(paper == 1)
%       set(gca,'XTickLabel',[0 30 60 90]);
%     else
%       set(gca,'XTickLabel',[0 40 80]);
%     end
%   end

  kdat = kdat+1;
  if kdat>nfreq;break;end;
end
saveas(hand(ifig),plotfile(ifig,:),'epsc2');
end

return
% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
