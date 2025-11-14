function [] = rep_vs_data_har();
%
%PLOT reflection model and data 
%
paper = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);

nlay = 2;

tab = readtextfile('measured_parameters_3_750-1265.dat');
[outfile,F,converg_cov,converg_marg] =...
convert_parfile(tab);

repfile = strrep(outfile,'sample.dat','rep_data_3_750-1265.mat');
ml_file = strrep(outfile,'sample.dat','ml_std_3.mat');
harfile = strrep(outfile,'sample.dat','har_data_750-1265.mat');
plotfile1 = strrep(outfile,'sample.dat','rep_vs_data_3.eps');

load(repfile);
Frep = B;
load(harfile);
for i = 1:length(F(1).freq)
  F(i).dat = B(i).dat;
end
F(1).ang = [0:2:90];
nfreq = length(fr);
F(1).freq

figure(1);
if(paper == 1)
  llw = [0.10 0.335 0.57 0.805 0.10 0.335 0.57 0.805 0.10 0.335 0.57 0.805];
  llh = [0.69 0.42 0.15];
else
  llw = [0.15 0.43 0.71 0.15 0.43 0.71 0.15 0.43 0.71];
  llh = [0.69 0.42 0.15];
end
spw = 0.195
sph = 0.26

i = 1;
j = 1;
for k=1:nfreq
   if k <= 4
     subplot('Position',[llw(i) llh(1) spw sph]);
     i = i + 1;
   elseif k <= 8
     subplot('Position',[llw(i) llh(2) spw sph]);
     i = i + 1;
   else
     subplot('Position',[llw(i) llh(3) spw sph]);
     i = i + 1;
   end
   hold on;box on
%   Err = ml_sd(k)*ones(nang(k),1);
   plot(F(1).ang(F(1).astart:F(1).aend), F(k).dat(F(1).astart:F(1).aend),'k.');
%   errorbar(F(k).ang,F(k).dat,Err,'k.');
   hold on;
   
   if(paper == 1)
     plot(F(1).ang(F(1).astart:F(1).aend),Frep(k).dat(F(1).astart:F(1).aend),'b');
   else
     plot(F(1).ang(F(1).astart:F(1).aend),Frep(k).dat(F(1).astart:F(1).aend),'b');
   end
%   title([num2str(freq(k)) ' Hz'] )
   text(55,0.9,[num2str(fr(k)) ' Hz'])
   axis([0 90 0 1]);
   legend('har','rep'); 
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
end
saveas(gca,plotfile1,'epsc2');

return
% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
