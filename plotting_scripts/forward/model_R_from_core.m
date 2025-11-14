function [] = model_R();
set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);

load core.mat;

fstep = 75.;
freq = [100.:fstep:10000.];
nfreq = length(freq);
ang = [1.:.5:90.];
nang = length(ang);

%----------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

  rhos = r1(:,2);
  cs = c1(:,2);
  hs = [c1(1,1); diff(c1(:,1))];
  as = 0.02*ones(size(hs));

  geo_sin=[    NaN     1511.         0.      1.029;  ...
               hs      cs           as      rhos; ...
               NaN     cs(end)       as(end)       rhos(end)];
%geo_sin

gc1 = figure(1);hold on;
%subplot('Position',[left bottom width height])
subplot('Position',[0.05 .61 .3 .3]);hold on;
plot(c1(:,2),c1(:,1),'r');
%plot(c2(:,2),c2(:,1),'r');
set(gca,'YDir','reverse','XLim',[1450 1550],'YLim',[0 2],'layer','top');
box on;
text(1450,-.2,['alpha = ',num2str(as(end)),'dB/L']) 


subplot('Position',[0.05 .05 .3 .3]);hold on;
plot(r1(:,2),r1(:,1),'r');
%plot(r2(:,2),r2(:,1),'r');
set(gca,'YDir','reverse','Xlim',[1.2 1.8],'YLim',[0 2],'layer','top');
box on;

for ifreq= 1:nfreq

  fr = [freq(ifreq)-fstep/2.:fstep/8.:freq(ifreq)+fstep/2.];
  for i = 1:8   
     [ref(i,:)] = ref_nlay3(ang,geo_sin,fr(i));% compute BL
  end
  ref2(ifreq,:) = sqrt(mean(ref, 1));
  ref2(ifreq,:) = -20*log10(abs(ref2(ifreq,:)));
end

subplot('Position',[0.42 .25 .58 .58]);hold on;
pcolor(ang,freq,ref2);shading flat;
set(gca,'XLim',[1 90],'YLim',[100 10000],'layer','top','CLim',[0 35]);
xlabel('Grazing angle (deg.)');
ylabel('Frequency (Hz)');
colorbar;
box on;

save('model_R_from_core.mat','ref2','ang','freq');
saveas(gc1,'model_R_from_core.png','png');

%%
%% Plot to compare..
%%
% idx=[1:5:length(freq)];for i=2:25;subplot(5,5,i);hold on;box on;plot(ang,ref2_c(idx(i),:));plot(ang,ref2_g(idx(i),:),'--k');text(10,1,num2str(freq(idx(i))));set(gca,'xtick',[],'ytick',[]);end;legend('core','grad model');

return;
