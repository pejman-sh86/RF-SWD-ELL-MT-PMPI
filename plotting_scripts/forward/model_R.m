function [] = model_R();
set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);

%S20:
%m = [1.0, 1479.0, 1478.0, 1.32, 1.57, 0.15, 0.03];
%m = [1.0, 1479.0, 1478.0, 1.56, 1.57, 0.5, 0.03];

%S19
%m = [1.0, 1479.0, 1468.0, 1.3, 1.5, 0.15, 0.03];
m = [1.0, 1479.0, 1468.0, 1.49, 1.5, 0.5, 0.03];

load core.mat;

znorm = [0:1./40.:1.];
sz = length(znorm);

fstep = 75.;
freq = [100.:fstep:10000.];
nfreq = length(freq);
ang = [1.:.5:90.];
nang = length(ang);

%----------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

  % rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
  rhos = m(4) + sin(znorm*pi/2).^m(6)*(m(5)-m(4));
  % cs=ct+(cb-ct)*znorm;
  cs=m(2)+(m(3)-m(2))*znorm;

  geo_sin=[    NaN     1511.         0.      1.029;  ...
           m(1)*1/20.*ones(sz,1) cs'  m(7)*ones(sz,1)  rhos'; ...
               NaN     m(3)       m(7)       m(5)];
%geo_sin

gc1 = figure(1);hold on;
%subplot('Position',[left bottom width height])
subplot('Position',[0.05 .61 .3 .3]);hold on;
plot(c1(:,2),c1(:,1),'r');
%plot(c2(:,2),c2(:,1),'r');
plot(geo_sin(2:end-1,2),znorm*m(1),'k');
set(gca,'YDir','reverse','XLim',[1450 1550],'YLim',[0 2],'layer','top');
box on;
text(1450,-.30,'h      c1       c2      r1      r2      nu      a') 
text(1450,-.2,num2str(m)) 


subplot('Position',[0.05 .05 .3 .3]);hold on;
plot(r1(:,2),r1(:,1),'r');
%plot(r2(:,2),r2(:,1),'r');
plot(geo_sin(2:end-1,4),znorm*m(1),'k');
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

save('model_R.mat','ref2','ang','freq');
saveas(gc1,'model_R.png','png');

return;
