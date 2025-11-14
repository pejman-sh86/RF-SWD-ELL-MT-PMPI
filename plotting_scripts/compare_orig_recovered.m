function [] = compare_orig_recovered();
%
% PLOT reflection model and data 
%
paper = 0;
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz

sample = 'sample.dat';
plotfile1 = strrep(sample,'sample.dat','compare_cov.eps');

%load six_cov_mat_est3;
%GG = F;
load real2_cov_mat_est2;
H = F;
load cov_mat_est1c;
%G = F;
%load cov_mat_est2;
%load ../full_cov_spher/cov_mat_est0;

freq = [315 400 500 630 800 1000 1250 1600];
%freq = [315 400 500 630 800 1000 1250 1600 2000 2500];
%freq = [500 630 800 1000 1250 1600 2000 2500 3160];
nfreq = length(freq);

figure(1);
if(paper == 1)
  width = 0.26;
  height = 0.24;
  llw = [0.10 0.40 0.7  0.10 0.40 0.7  0.10 0.40 0.7 ];
  llh = [0.75 0.44 0.13];
else
  width = 0.26;
  height = 0.24;
  llw = [0.15 0.43 0.71 0.15 0.43 0.71 0.15 0.43 0.71];
  llh = [0.75 0.44 0.13];
end
i = 1;
j = 1;
for k=1:nfreq
   if k <= 3
     subplot('Position',[llw(i) llh(1) width height]);
     i = i + 1;
   elseif k <= 6
     subplot('Position',[llw(i) llh(2) width height]);
     i = i + 1;
   else
     subplot('Position',[llw(i) llh(3) width height]);
     i = i + 1;
   end
   hold on;box on;
   
   if(paper == 0)
     plot(F(k).csave(1,:),'r');
%     plot(G(k).csave(1,:),'-b');
     plot(H(k).csave(1,:),'--b');
%     plot(GG(k).csave(1,:),':b');
   else
     plot(F(k).csave(1,:),'r');
     plot(G(k).csave(1,:),'-b');
%     plot(H(k).csave(1,:),'--b');
%     plot(GG(k).csave(1,:),':b');
   end
   legend('0','1','2','3');
   text(2*length(F(k).csave(1,:))/3,4.8,[num2str(freq(k)) ' Hz'],'FontSize',10);
%   text(55,40,[num2str(freq(k)) ' Hz'])
   if(paper == 1)
     set(gca,'XTick',[50 100],'XLim',[0 length(F(k).csave(1,:))],...
     'YTick',[-1 0 1 2 3 4 5],'YLim',[-1.5 6],'FontSize',12);
   else
     set(gca,'XTick',[50 100],'XLim',[0 length(F(k).csave(1,:))],...
     'YTick',[-1 0 1 2 3 4 5],'YLim',[-1.8 6],'FontSize',12);
   end
   if k >=  6;xlabel('Lag');end
   if k == 1;ylabel('A_{nn}');end
   if k == 4;ylabel('A_{nn}');end
   if k == 7;ylabel('A_{nn}');end
   if(paper == 1)
     set(gca,'XTickLabel',{'50';'100'});
   else
     set(gca,'XTickLabel',{'50';'100'});
   end
   if k == 1; set(gca,'YTickLabel',{'';'0';'';'2';'';'4';''});end
   if k == 2; set(gca,'YTickLabel',{'';'';'';'';'';'';''});end
   if k == 3; set(gca,'YTickLabel',{'';'';'';'';'';'';''});end
   if k == 4; set(gca,'YTickLabel',{'';'0';'';'2';'';'4';''});end
   if k == 5; set(gca,'YTickLabel',{'';'';'';'';'';'';''});end
   if k == 6; set(gca,'YTickLabel',{'';'';'';'';'';'';''});end
   if k == 7; set(gca,'YTickLabel',{'';'0';'';'2';'';'4';''});end
   if k == 8; set(gca,'YTickLabel',{'';'';'';'';'';'';''});end
end
saveas(gca,plotfile1,'epsc2');
% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
