
%load R_20090522_p175_2752_psCMP_fringe.mat;
load R_20090522_p175_2752_ALL_psCMP;

R=RR_psCMP;

idx = [8, 16, 64];
figure(3);box on;hold on;
subplot(3,1,1);box on;
pcolor([1:3433],360/2/pi*ANG_all_tilt(:,1),squeeze(R(idx(1),:,:)));
shading flat;
subplot(3,1,2);box on;
pcolor([1:3433],360/2/pi*ANG_all_tilt(:,1),squeeze(R(idx(2),:,:)));
shading flat;
subplot(3,1,3);box on;
pcolor([1:3433],360/2/pi*ANG_all_tilt(:,1),squeeze(R(idx(3),:,:)));
shading flat;

ipskip = 10;
for k=1:ipskip:3433;
   figure(1);box on;hold on;
   pcolor(360/2/pi*ANG_all_tilt(:,k),freq,R(:,:,k));
   shading flat;colorbar;
   set(gca,'XLim',[30 70],'YLim',[800 3200],'CLim',[0 .6],'YTick',[900:100:3200]);
   plot([30 70],[950 950],'-k');plot([30 70],[1300 1300],'-k');
   plot([30 70],[1850 1850],'-k');plot([30 70],[2400 2400],'-k');
   hold off;

%   figure(2);box on;
%   for j=1:3;
%      subplot(1,3,j);box on;
%      plot(360/2/pi*ANG_all_tilt(:,k),R(idx(j),:,k));
%      shading flat;
%      set(gca,'XLim',[30 70],'YLim',[0 .6]);
%   end;
   pause(.2);
end;


