function [] = compare_data_hol_lay_R();
%
%PLOT reflection model and data 
%
paper = 1;
isave = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 9 6]);

ml_file = 'ml_std_1.mat';
plotfile(1,:) = 'rep_vs_data_0a.eps';
plotfile(2,:) = 'rep_vs_data_0b.eps';

map = load('spher_refldBa.dat');
true = load('spher_refldBb.dat');
%sdev = load('sim_A_4_sdev.txt');
%load('sim_A_sph_4.mat');
%Fsph = F;
%load('sim_A_4.mat');
load('x_jan_1_15-30.mat');
nfreq = length(F(1).freq);

ndat = length(map)-1;


%for ifreq = 1:nfreq
%
%    Fsph(ifreq).R = 10.^(-Fsph(ifreq).dat/20);
%
%end

nx = 2
ny = 3
nsubpfig = 5;
nfig = 1;
xim = 0.06;
yim = 0.06;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

kdat = 1;
for ifig = 1:nfig

hand(ifig) = figure(ifig);hold on;
i = 1;

for k=1:nfreq
%for k= [1 4]

   
   subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   hold on;box on;

   if k <= 2
     set(gca,'XTickLabel',[]);
   elseif k <= 2
     set(gca,'XTickLabel',[]);
   else
     set(gca,'XTickLabel',[20 40 60 80]);
     xlabel('Angle (deg.)','FontSize',14);
   end
   if ((k == 1)|(k == 4)|(k == 7))
       set(gca,'YTickLabel',[0 0.3 0.6 0.9]);
       ylabel('|V|','FontSize',14);
   else
       set(gca,'YTickLabel',[]);
   end
   set(gca,'XTick',[20 40 60 80],'FontSize',14);
   set(gca,'YTick',[0 0.3 0.6 0.9]);
   hold on;box on
%   Err = mean(sdev(k,:));
%   Err = sdev(k);
   plot(F(k).ang(1:ndat), F(k).dat(1:ndat).*F(1).Rex(1:ndat,k),'xk');
%   plot(F(k).ang(1:5:ndat), F(k).dat(1:5:ndat),'.k');
%   dashline(F(kdat).ang, F(kdat).dat,'.k',60,'.k',60);
%   errorbar(F(k).ang(1:5:ndat),F(k).dat(1:5:ndat),Err,'k.');
   hold on;
   
%   plot(map(6,1:ndat),map(kdat,1:ndat),'-k','LineWidth',0.5);
   plot(F(k).ang(1:ndat),map(k,1:ndat),'--k','LineWidth',0.5);
   plot(F(k).ang(1:ndat),true(k,1:ndat),':k');
%   plot(Fsph(kdat).ang,Fsph(kdat).R,':k');
%   plot(F(kdat).ang,F(kdat).dat,':k');
%   title([num2str(freq(k)) ' Hz'] )
   text(65,.9,[num2str(F(1).freq(k)) ' Hz'])
   axis([10 85 0 1.3]);
%   legend('simA','pl map','pl true','sp true','Location','SouthEast'); 

  kdat = kdat+1;
  i = i+1;
  if kdat>nfreq;break;end;
end
  if(isave == 1)
    saveas(hand(ifig),plotfile(ifig,:),'epsc2');
  end
end

return
% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
