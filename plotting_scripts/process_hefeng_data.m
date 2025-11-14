function [] = process_hefeng_gata();

load('VfVpFp_dispersion_PH1-80Greentw4p5secfp68-6.mat');
vp=vp(300:end,:);
fp=fp(:,300:end,:);
vp=vp(1:201,:);
fp=fp(:,1:201,:);


vf1 = vf(find(vf<=2.5))/max(vf(find(vf<=2.5)))*pi

1000.*(cos(vff)+1.)/2.+600.;




fp2 = fp;
for i=1:50;
for j=1:64;
  fp2(j,1:limlo(j),i)=0.;
end;end;

for i=1:50;
for j=1:64;
  fp3(j,:,i)=fp2(j,:,i)/max(fp2(j,:,i));
end;end;


for i=1:50;
  fig1=figure(1);
  pcolor(vf,vp(:,i),fp3(:,:,i)');
  shading flat;
  set(gca,'YLim',[600 1600],'XLim',[0 3.8]);
  saveas(fig1,strcat('data_',num2str(i,'%02i'),'.png'),'png');
  pause(.1);
end;




return;
