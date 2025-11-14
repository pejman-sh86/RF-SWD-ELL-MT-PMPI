load otofc.dat;  load otofedg.dat;
frq=otofc(fr_otonum)
foff=4 - 20; frcl=foff+fr_otonum;

for n=1:numfrq
% LOAD DIRECT PATH TL MODEL RESULTS
  fn= ['bl' num2str(fr_otonum(n)) '.dr' phone]; fname1= [dirnmtl fn];
  load (fname1);  tmp=eval(fn(1:4));
  rangemax=max(tmp(:,6));
  
% CHECK THAT SNR > snr FOR DIRECT PATH DATA on positive ranges & limit ranges to less than maxrange m or the rangemax for modeled TL
[ins2,ysnrall] = fit_noise2(tldsnr,frcl(n),snr,minrange,maxrange,rangemax); %sliding window
%[ins2,ysnrall] = fit_noise1(tldsnr,frcl(n),snr,minrange,maxrange,rangemax); %polynomial fit

  for ii=1:4; %do several passes to remove outliers
     
% COMPUTE ANGLES AT SOURCE
  yang=interp1(tmp(:,6),tmp(:,1),abs(rld(ins2,1)));
  thb=yang*pi/180;
  ths=acos(cos(thb)*cs/cb)*180/pi;
  
   yall=interp1(tmp(:,6),tmp(:,1),abs(rld(min(ins2):max(ins2),1)));
  thball=yall*pi/180; thsall=acos(cos(thball)*cs/cb)*180/pi;

% INTERPOLATE TL VALUES ON DATA GRID AND COMPUTE SL
  ytl=interp1(tmp(:,6),tmp(:,4),abs(rld(ins2,1)));
  %sl_bp=yrl+ytl;
  sl_bp=rld(ins2,frcl(n))+ytl;

  if frq(n)> 50 & frq(n)< 3500; [psl ssl mu]=polyfit(ths,sl_bp,poly_num);  
else; [psl ssl mu]=polyfit(ths,sl_bp,poly_num); end     
 %  [psl ssl mu]=polyfit(ths,sl_bp, poly_num);   
 
  [ysl delta]=polyval(psl,ths,ssl,mu);
   nrat(n)=length(find( abs(sl_bp-ysl)>delta))/length(sl_bp);
  [yslall delta2]=polyval(psl,thsall,ssl,mu);
  stdev(n)=std(sl_bp-ysl);
  nrat_dev(n)=length(find( abs(sl_bp-ysl)>stdev(n)))/length(sl_bp);
  
  % ysl=polyval(psl,[1:90]);
% ysltot=[ysltot ysl];
 sld(min(ins2):max(ins2),frcl(n))=yslall;
 rldscr(ins2,frcl(n))=rld(ins2,frcl(n));

%WINDOW OUT DATA THAT FALL "DETDB" OUTSIDE OF FIT
  outliers= find(abs(ysl-sl_bp)>2*delta);
  ins2=setdiff(ins2,ins2(outliers));
end

 numplot=ceil(numfrq/3);
 figure(20+ii); orient tall; subplot(numplot,3,n); box on; hold on;
  plot(ths,sl_bp,'k.'); hold on;  plot(ths,ysl,colr2)
  plot(ths(outliers),sl_bp(outliers),'g.'); 

if n==floor(numfrq/2);  ylabel('SL (dB re 1 \muPa^2s/Hz  at 1m)'); end 
      if n >(numplot-1)*3 & n <numfrq+1;   
         set(gca,'XTick',[0:15:90]); xlabel('Source Angle (deg)');
      else  set(gca,'XTick',[0:15:90],'XTicklabel',[ ]); end   
      title([num2str(frq(n)),' Hz']) 
%      axis([0 90 80 150])     
      
      
 figure(9); orient tall; subplot(numplot,3,n); box on; hold on;
  plot(ths,sl_bp-ysl,['k.']); hold on;
  plot(thsall,delta2,'b'); plot(thsall,-delta2,'b')
%    if size(intersect(n,ylab))>0;  ylabel('SL (dB)'); end   
    if n==floor(numfrq/2);  ylabel('Residual SL error(dB re 1 uPa^2s/Hz)'); end   
      if n >(numplot-1)*3 & n <numfrq+1;   
         set(gca,'XTick',[0:30:90]); xlabel('Source Angle (deg)');
      else  set(gca,'XTick',[0:30:90],'XTicklabel',[ ]); end   
%      if n==numfrq; legend('measured','piston model','fit'); end
      title([num2str(frq(n)),' Hz']) 
      
%            figure(6)
%       box on; hold on;
%  plot(ths,sl_bp,'k.'); hold on;
%  plot(ths,ysl,colr2);xlabel('Source Angle (deg)');
%ylabel('SL (dB re 1 uPa^2s/Hz)');

%  figure(4); orient tall; subplot(numplot,3,n);  hold on;box on;
%  plot(rld(ins2,1),rld(ins2,frcl(n)),'k.' );
%  plot(rldscr(:,1),-rldscr(:,frcl(n))+rldscr6msfft(:,frcl(n)),'r.' );
%plot(rldscr(:,1),rldscr(:,frcl(n)),'r.' );
%      if size(intersect(n,ylab));  ylabel('RL (dB)'); end   
%      if n >(numplot-1)*3 & n <numfrq+1;  xlabel('Range (m)'); 
%      else  set(gca,'XTicklabel',[ ]); end   
%      if n==numfrq; legend('measured','meas fit'); end
%      title([num2str(frq(n)),' Hz']) 
%       axis([minrange maxrange -5 5]);
      
if iplot==2 
 figure(7); subplot(numplot,3,n); hold on;
  plot(tldsnr(:,1),tldsnr(:,frcl(n)),colr2 );
  plot(tldsnr(:,1),ysnrall,'r')
  plot([-maxrange maxrange],[6 6],'k');
%  plot(rld(ins2,1),sl_bp-40,'b')
      if size(intersect(n,ylab))>0;  ylabel('SNR (dB)'); end   
      if n >(numplot-1)*3 & n <numfrq+1; xlabel('Range (m)'); 
      else  set(gca,'XTicklabel',[ ]); end   
      if n==numfrq; legend('measured','meas fit'); end
      title([num2str(frq(n)),' Hz']) 
       axis([minrange maxrange 0 20]);
       
  figure(3); orient tall; subplot(numplot,3,n);  hold on;
%  plot(rld(ins2,1),rld(ins2,frcl(n)),'k.' );
 % plot(rld(ins2,1),rld(ins2,frcl(n)),'r.' );
  plot(rld(:,1),rld(:,frcl(n)),'k.' );
  plot(rldscr(:,1),rldscr(:,frcl(n)),'c.' );
%  plot(rld(ins2,1),yrl,'g')
%  plot(rld(ins2,1),sl_bp-40,'b')
      if size(intersect(n,ylab));  ylabel('RL (dB)'); end   
      if n >(numplot-1)*3 & n <numfrq+1;  xlabel('Range (m)'); 
      else  set(gca,'XTicklabel',[ ]); end   
      if n==numfrq; legend('measured',' > SNR thresh'); end
      title([num2str(frq(n)),' Hz']) 
       axis([minrange maxrange ymin ymax]);
       
  figure(6); orient tall; subplot(numplot,3,n); hold on;
  plot(rld(ins2,1),rld(ins2,frcl(n)),'k.' ); 
  plot(tmp(:,6),tmp(:,4),colr2)
%  plot(rld(ins2,1),ytl,'g')
      if n==1 | n==4 | n==7;  ylabel('TL (dB)'); end   
      if n >6 & n <numfrq+1; xlabel('Range (deg)');  
      else  set(gca,'XTicklabel',[ ]); end   
      if n==numfrq; legend('measured','omni model','omni fit'); end
      title([num2str(frq(n)),' Hz']) 
      axis('ij'); axis([minrange maxrange 60 110])     
    
    figure(4); subplot(numplot,3,n); hold on;
%  plot(tmp(:,1),tmp(:,4),'b')
	plot(ths,ytl,colr2)
      if n==1 | n==4 | n==7;  ylabel('TL (dB)'); end   
      if n >6 & n <numfrq+1;   xlabel('Angle (deg)');
      else  set(gca,'XTicklabel',[ ]); end   
      if n==numfrq; legend('measured','omni model','omni fit'); end
      title([num2str(frq(n)),' Hz']) 
      axis('ij'); axis([0 90 20 90])     
    
    end  
 end
 
 % for i=1:3, subplot(3,3,i); axis([0 90 120 155]);end
 
 %NOW WRITE OUT SOURCE LEVEL FILE 
 %IF SOME COLUMNS ARE NOT REPRESENTED BECAUSE SNR<snr THEN
 % FILL IN COLUMNS AS (BLSIM_BP expectes columns from band 20-40)
   yangall=interp1(tmp(:,6),tmp(:,1),abs(rld(:,1)));
 sld(:,3)=yangall;
 
 if minrange<0; fsnm=[dirnm 'bl data\hp_results\' fnbl '_sl_' id '_neg.asc']
 else fsnm=[dirnm 'bl data\hp_results\' fnbl '_sl_' id '.asc'];end
    
%  save(fsnm, 'sld', '-ascii')

if iplot==3
   figure(8)
 for i=[1:18]; subplot(7,3,i);axis([0 90 105 140]);end
 for i=[16]; subplot(7,3,i);axis([0 90 100 135]);end
 for i=[17:19]; subplot(7,3,i);axis([0 90 95 130]);end
 end 
 
 figure
semilogx(frq,stdev,'b');hold on
semilogx(frq,1-nrat_dev,'r')
xlabel('Frequency (Hz)');ylabel('Standard Deviation (dB)')
 axis([100 10000 0 2])
 
if 1==2
%Calculate error as a function of angle
 fac=factor(length(ths));
 if max(fac)>8; fac1=max(fac);fac2=length(ths)/fac1;
   else; fac1=prod(fac(end-1:end));fac2=length(ths)/fac1; end    
 for i=1:6
    pp=((1:fac1)+((i-1)*fac1) );mn_ths(i)=mean(ths(pp)); std_bl(i)=std(sl_bp(pp)-ysl(pp));
 end
  plot([mn_ths(1:end-1) ths(end)],-std_bl); plot([mn_ths(1:end-1) ths(end)],std_bl)
  axis([0 90 -6 6]); set(gca,'Xtick',[0:30:90])
  xlabel('Grazing Angle (deg)'); ylabel('Res Error (dB)')
end