function []=plot_pecan_obsvsrep(filename);

NSRC  = 5;
NBAND = 3;
bands = [50, 100, 200];
filebase     = strrep(filename,'sample.mat','');
plotfiledatr = strcat(filebase,'datfitr.');
plotfiledati = strcat(filebase,'datfiti.');
datfile      = strcat(filebase,'s_03_f_50_200_sim.txt');
repfile      = strcat(filebase,'replica.dat');
%repfile2     = strcat(filebase,'replica2.dat');
%repfile3     = strcat(filebase,'replica3.dat');
plotext2     = 'png';
plotext3     = 'eps';

dat = load(datfile);
rep = load(repfile);
%rep2= load(repfile2);
%rep3= load(repfile3);
z = rep(1,:);
rep = rep(2:end,:);

%% Assemble complex data and compute amplitude scaling
j = 1;
for isrc = 1:NSRC;
for iband = 1:NBAND;
  datc(iband,:,isrc) = complex(dat(j,:),dat(j+1,:));
  repc(iband,:,isrc) = complex(rep(j,:),rep(j+1,:));
  j = j+2;
end;
end;
NHYD = size(datc,2);
% Compute scaling and residuals
for isrc = 1:NSRC;
for iband = 1:NBAND;
  A(iband,isrc) = sum( conj(repc(iband,:,isrc)).*datc(iband,:,isrc) ) / ...
                  sum( conj(repc(iband,:,isrc)).*repc(iband,:,isrc) );
  repc2(iband,:,isrc) = A(iband,isrc)*repc(iband,:,isrc);
  sd(iband,isrc) = sqrt((sum( conj(datc(iband,:,isrc)).*datc(iband,:,isrc) ) - ...
                   sum( conj(conj(repc(iband,:,isrc)).*datc(iband,:,isrc)).* ...
		            (conj(repc(iband,:,isrc)).*datc(iband,:,isrc)) ) / ...
                   sum( conj(repc(iband,:,isrc)).*repc(iband,:,isrc) ))/NHYD);
end;
end;
for isrc = 1:NSRC;
for iband = 1:NBAND;
  res(iband,:,isrc) = datc(iband,:,isrc)-repc2(iband,:,isrc);
  SNR(iband,isrc) = 10*log10(sum(conj(repc2(iband,:,isrc)).*repc2(iband,:,isrc))/sum(conj(res(iband,:,isrc)).*res(iband,:,isrc)));
end;
end;
SNR

sd
nx = NBAND;
ny = NSRC;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
ifr = 0;
isrc = 0;
figdatr=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
figdati=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
i = 0;
for isrc=1:NSRC;
for iband=1:NBAND;
  i = i + 1;

  figure(figdatr);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;

  tmp1 = real(datc(iband,:,isrc))/max(abs(real(datc(iband,:,isrc))));
  %plot(z,tmp1,'xk');
  tmp2 = real(repc2(iband,:,isrc))/max(abs(real(datc(iband,:,isrc))));
  plot(z,tmp2,'--k','LineWidth',1);
  sdvec = std(tmp1-tmp2)*ones(size(tmp1));
  H=errorbar(z,tmp1,sdvec,'ok','MarkerSize',2,'LineWidth',1,'MarkerFaceColor','k');
  errorbar_tick(H,80);

  figure(figdati);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  tmp1 = imag(datc(iband,:,isrc))/max(abs(imag(datc(iband,:,isrc))));
  %plot(z,tmp1,'xk');
  tmp2 = imag(repc2(iband,:,isrc))/max(abs(imag(datc(iband,:,isrc))));
  plot(z,tmp2,'--k','LineWidth',1);
  sdvec = std(tmp1-tmp2)*ones(size(tmp1));
  H=errorbar(z,tmp1,sdvec,'ok','MarkerSize',2,'LineWidth',1,'MarkerFaceColor','k');
  errorbar_tick(H,80);

  figure(figdatr);
  set(gca,'YLim',[-1.2 1.2],'XLim',[0 140],'LineWidth',1);
  set(gca,'YTick',[-.8:.4:.8]);
  if(isrc == NSRC);xlabel('Depth (m)');else;set(gca,'XTickLabel',[]);end;
  if(iband > 1);set(gca,'YTickLabel',[]);
  else;ylabel('Norm. pressure');end;
  text(5,.9,['Source' num2str(isrc) ',' num2str(bands(iband)) 'Hz'],'FontSize',12,'Color',[0,0,0]);
  figure(figdati);
  set(gca,'YLim',[-1.2 1.2],'XLim',[0 140],'LineWidth',1);
  set(gca,'YTick',[-.8:.4:.8]);
  if(isrc == NSRC);xlabel('Depth (m)');else;set(gca,'XTickLabel',[]);end;
  if(iband > 1);set(gca,'YTickLabel',[]);
  else;ylabel('Norm. pressure');end;
  text(5,.9,['Source' num2str(isrc) ',' num2str(bands(iband)) 'Hz'],'FontSize',12,'Color',[0,0,0]);
end;
end;

ifr = 0;
isrc = 0;
figresr=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
figresi=figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
i = 0;
for isrc=1:NSRC;
for iband=1:NBAND;
  i = i + 1;

  figure(figresr);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  tmp1 = real(res(iband,:,isrc))%;/max(abs(real(res(iband,:,isrc))));
  plot(z,tmp1,'.-k');

  figure(figresi);
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  tmp1 = imag(res(iband,:,isrc));%/max(abs(imag(res(iband,:,isrc))));
  plot(z,tmp1,'.-k');

  figure(figresr);
  set(gca,'YLim',[-.03 .03],'XLim',[0 140],'LineWidth',1);
  set(gca,'YTick',[-.03:.01:.03]);
  if(isrc == NSRC);xlabel('Depth (m)');else;set(gca,'XTickLabel',[]);end;
  if(iband > 1);set(gca,'YTickLabel',[]);
  else;ylabel('Norm. pressure');end;
  text(5,.9,['Source' num2str(isrc) ',' num2str(bands(iband)) 'Hz'],'FontSize',12,'Color',[0,0,0]);
  figure(figresi);
  set(gca,'YLim',[-.03 .03],'XLim',[0 140],'LineWidth',1);
  set(gca,'YTick',[-.03:.01:.03]);
  if(isrc == NSRC);xlabel('Depth (m)');else;set(gca,'XTickLabel',[]);end;
  if(iband > 1);set(gca,'YTickLabel',[]);
  else;ylabel('Norm. pressure');end;
  text(5,.9,['Source' num2str(isrc) ',' num2str(bands(iband)) 'Hz'],'FontSize',12,'Color',[0,0,0]);
end;
end;


print(figdatr,'-painters','-r250',strcat(plotfiledatr,plotext3),'-depsc');
print(figdatr,'-painters','-r250',strcat(plotfiledatr,plotext2),'-dpng');
print(figdati,'-painters','-r250',strcat(plotfiledati,plotext3),'-depsc');
print(figdati,'-painters','-r250',strcat(plotfiledati,plotext2),'-dpng');
return;
