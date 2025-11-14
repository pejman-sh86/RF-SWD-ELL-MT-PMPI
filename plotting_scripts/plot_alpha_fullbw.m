function []=plot_alpha_fullbw();

load ../band1/nf_band1.mat;
alfrtmean1=alfrtmean; frt1 = frt;nfa1=nfa;nfv1=nfv;vptmean1=vptmean;
load ../band2/nf_band2.mat;
alfrtmean2=alfrtmean; frt2 = frt;nfa2=nfa;nfv2=nfv;vptmean2=vptmean;
load ../band3/nf_band3.mat;
alfrtmean3=alfrtmean; frt3 = frt;nfa3=nfa;nfv3=nfv;vptmean3=vptmean;
load ../band4/nf_band4.mat;
alfrtmean4=alfrtmean; frt4 = frt;nfa4=nfa;nfv4=nfv;vptmean4=vptmean;

load Biot_halfspace_test.mat;

figure();
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',16)
plot(frt1,vptmean1,'b','Linewidth',2);
plot(frt2,vptmean2,'b','Linewidth',2);
plot(frt3,vptmean3,'b','Linewidth',2);
plot(frt4,vptmean4,'b','Linewidth',2);
plot(frt1,nfv1(:,1),'-k');
plot(frt2,nfv2(:,1),'-k');
plot(frt3,nfv3(:,1),'-k');
plot(frt4,nfv4(:,1),'-k');
plot(frt1,nfv1(:,2),'-k');
plot(frt2,nfv2(:,2),'-k');
plot(frt3,nfv3(:,2),'-k');
plot(frt4,nfv4(:,2),'-k');
plot(biot_f,biot_V(1,:),'r','Linewidth',2);
set(gca,'XScale','log','XLim',[frt1(1)-1 frt4(end)+1]);
xlabel('Frequency (Hz)');
ylabel('Sound velocity (m/s)');

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',16)
h1=plot(frt1,alfrtmean1,'b','Linewidth',2);
plot(frt2,alfrtmean2,'b','Linewidth',2);
plot(frt3,alfrtmean3,'b','Linewidth',2);
plot(frt4,alfrtmean4,'b','Linewidth',2);
h2=plot(frt1,nfa1(:,1),'-k');
plot(frt2,nfa2(:,1),'-k');
plot(frt3,nfa3(:,1),'-k');
plot(frt4,nfa4(:,1),'-k');
plot(frt1,nfa1(:,2),'-k');
plot(frt2,nfa2(:,2),'-k');
plot(frt3,nfa3(:,2),'-k');
plot(frt4,nfa4(:,2),'-k');
h3=plot(biot_f,biot_A(1,:)./biot_f,'r','Linewidth',2);
set(gca,'XScale','log','XLim',[frt1(1)-1 frt4(end)+1]);
xlabel('Frequency (Hz)');
ylabel('Attenuation (nepers/m/Hz)');
legend([h1,h2,h3],'inversion','inversion credibility','true')

return;
