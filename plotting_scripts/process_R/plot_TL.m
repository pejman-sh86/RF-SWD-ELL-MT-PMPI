function []=plot_TL();

NSMP = 1250;
NRNG = 1000;
rep_ar   = dlmread('rjmh_mfi_replica_ar.dat');
rep_noar = dlmread('rjmh_mfi_replica_noar.dat');
range = [.1:.01:10.09];

for i=1:NSMP;
  TL1_400(i,:) = 20*log10(rep_ar(((i-1)*NRNG)+1:i*NRNG,1));
  TL2_400(i,:) = 20*log10(rep_noar(((i-1)*NRNG)+1:i*NRNG,1));
  TL1_800(i,:) = 20*log10(rep_ar(((i-1)*NRNG)+1:i*NRNG,2));
  TL2_800(i,:) = 20*log10(rep_noar(((i-1)*NRNG)+1:i*NRNG,2));
  if(isinf(TL1_400(i,:)) == 1);TL1_400(i,:)=[];end;
  if(isinf(TL2_400(i,:)) == 1);TL2_400(i,:)=[];end;
  if(isinf(TL1_800(i,:)) == 1);TL1_800(i,:)=[];end;
  if(isinf(TL2_800(i,:)) == 1);TL2_800(i,:)=[];end;
end;

%% Compute HPD as function of range:
for i=1:NRNG;

  [nf1_400(i,:)] = hpd(TL1_400(:,i),100,95);
  [nf2_400(i,:)] = hpd(TL2_400(:,i),100,95);
  [nf1_800(i,:)] = hpd(TL1_800(:,i),100,95);
  [nf2_800(i,:)] = hpd(TL2_800(:,i),100,95);

end;
figure(1);hold on; box on;
%for i=1:NSMP;
%  subplot(2,1,1);hold on; box on;
%  plot(range,TL1_400(i,:),'-b');
%  plot(range,TL2_400(i,:),'--r');
%  subplot(2,1,2);hold on; box on;
%  plot(range,TL1_800(i,:),'-b');
%  plot(range,TL2_800(i,:),'--r');
%end;
subplot(2,1,1);hold on; box on;
poly(:,1) = [nf1_400(:,1);flipud(nf1_400(:,2))];

size(range)
size(poly)

poly(:,2) = [range';flipud(range')];
fill(poly(:,2),poly(:,1),[0. 0. 1.],'LineStyle','none');

poly(:,1) = [nf2_400(:,1);flipud(nf2_400(:,2))];
fill(poly(:,2),poly(:,1),[1. 0. 0.],'LineStyle','none');

xlabel('Range (km)');
ylabel('Transmission Loss (dB)');
set(gca,'XLim',[0 10],'YLim',[-130 -40]);
subplot(2,1,2);hold on; box on;
poly(:,1) = [nf1_800(:,1);flipud(nf1_800(:,2))];
fill(poly(:,2),poly(:,1),[0. 0. 1.],'LineStyle','none');

poly(:,1) = [nf2_800(:,1);flipud(nf2_800(:,2))];
fill(poly(:,2),poly(:,1),[1. 0. 0.],'LineStyle','none');

xlabel('Range (km)');
ylabel('Transmission Loss (dB)');
set(gca,'XLim',[0 10],'YLim',[-130 -40]);

return;
