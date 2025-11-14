%function [] = auv_prune_data(filename);

filename_raw = 'R_20090522_p175_2752_ALL_psCMP.mat';

filebase = strrep(filename_raw,'.mat','');

i_save = 1;
i_inter = 0;  %% 1 = resample onto even grid
i_frave = 1;  %% 1 = average frequencies 
i_pave  = 1;  %% 1 = apply ping averaging
nfave = 1;    %% total freq to ave is (2*nfave)+1
nskip = 2;   %% stride between pings
npave = 5;   %% No. pings in average

sdcutoff = 0.04;

NANG = 32;
NPING = 3433;
nang_re = 31;
load(filename_raw);

%% Frequencies:
%idx1 = [7, 12, 18, 52, 64, 74];
idx1 = [7, 12, 19, 44, 58, 68];

idx2 = idx1+nfave;
%freq = [975,1100,1250,2100,2400,2700];
freqorig = freq;

freq1 = freq(idx1);
freq2 = freq(idx2);
freq = round((freq1+freq2)/2.);

%idxlo = idx-nfave;
%idxhi = idx+nfave;
idxlo = idx1;
idxhi = idx2;

nf = length(idx1);
pstart = npave+1;
discard = zeros(NPING,1);
pad = npave;
idxd = [346-pad: 372+pad,...
        428-pad: 440+pad,...
       1670-pad:1690+pad,...
       1717-pad:1727+pad,...
       1782-pad:1810+pad,...
       1980-pad:1990+pad,...
       2046-pad:2056+pad,...
       2300-pad:2319+pad,...
       2340-pad:2374+pad,...
       2492-pad:2530+pad,...
       2740-pad:2764+pad,...
       2866-pad:2870+pad,...
       3164-pad:3180+pad,...
       3210-pad:3230+pad,...
       3255-pad:3280+pad];
discard(idxd) = 1;   
col = {'r','g','b','c','m','y','k'};
col = char(col);

iiping = 0;
sd2=zeros((NPING-pad-pstart)/nskip,length(freq),NANG);
Rex=ones(size(sd2));

for iping = pstart:nskip:NPING-pad;
   if(discard(iping) == 1);continue;end;
   disp(iping);
   iiping = iiping + 1;
   idxping(iiping) = iping;
   filename = strcat('p',num2str(iping,'%04i'),'_','pave',num2str(2*npave,'%03i'),...
              '_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.txt');
   filenamea= strcat('p',num2str(iping,'%04i'),'_','pave',num2str(2*npave,'%03i'),...
              '_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.mat');
   plotname1 = strcat('p',num2str(iping,'%04i'),'_','pave',num2str(2*npave,'%03i'),...
              '_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.png');
   plotname2 = strcat('p',num2str(iping,'%04i'),'_','pave',num2str(2*npave,'%03i'),...
              '_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'sd','.png');

   cw = 1512.3;
   rw = 1.029;
   tb_win = 0.0085;
   hmx = tb_win*1600./2.;

   if(i_pave == 1);
     Rptmp = RR_psCMP(:,:,iping-npave:iping+npave);
     Rp = mean(Rptmp,3);
     theta = mean(ANG_all_tilt(:,iping-npave:iping+npave),2)*360./(2.*pi);
     z_t = 2.*mean(alt(iping-npave:iping+npave));
   else;
     Rptmp = RR_psCMP(:,:,iping);
     Rp = RR_psCMP(:,:,iping);
     theta = ANG_all_tilt(:,iping)*360./(2.*pi);
     z_t = 2.*alt(iping);
   end;

   if(i_frave == 1);
      for ifr=1:nf
        Rtmp = Rp(idxlo(ifr):idxhi(ifr),:);
        R(ifr,:) = sqrt(sum(Rtmp.*Rtmp,1)/(nfave+1));
      end;
   else;
      R = Rp(idx1,:);
   end;

   thdiff = diff(theta);
   if(i_inter == 1);
     theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];
     theta_even = theta_even';
     for ifr = 1:nf
       R_even(ifr,:)   = interp1(theta,R(ifr,:),theta_even);
     end
   else
      theta_even = theta;
      R_even = R;
   end;

   nx = 3;
   ny = ceil(nf/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   figw = 14;
   figh = 7;
   fig1 = figure(1);hold on; box on;
   fig2 = figure(2);hold on; box on;
   set(fig1, 'Position', [0 0 1000 500])
   set(fig2, 'Position', [0 0 1000 500])
   for ifr = 1:nf;
      set(0, 'currentfigure',fig1);
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold off;box on;
      idxtmp = [idxlo(ifr):idxhi(ifr)];
      idat = 1;
      for ip = 1:2*npave+1;
        for isub = 1:nfave+1;
          plot(theta_even,Rptmp(idxtmp(isub),:,ip),col(isub,:));hold on;
          res(idat,:) = Rptmp(idxtmp(isub),:,ip);
          idat = idat + 1;
        end;
      end;
      sd2(iiping,ifr,:) = std(res,1);
      sd1(iiping,ifr) = mean(sd2(iiping,ifr,:));
      sd2f(iiping,ifr,:) = (std(res(1:2:end,:),1)+std(res(2:2:end,:),1))/2.;
      sd1f(iiping,ifr) = mean(sd2f(iiping,ifr,:));
      dRdTh = (R_even(ifr,2:end)-R_even(ifr,1:end-1))./...
              ((theta_even(2:end)-theta_even(1:end-1))*pi/180)';
      dRdTh = [dRdTh(1),dRdTh];
      sd2p(iiping,ifr,:) = squeeze(sd2(iiping,ifr,:))./sqrt(1+(dRdTh).^2)';
      sd1p(iiping,ifr) = mean(sd2p(iiping,ifr,:));
      
      idxex = find(sd2p(iiping,ifr,:)>sdcutoff);
      Rex(iiping,ifr,idxex) = 0;
      clear idxex;
      idxex = find(theta > 55.);
      if(ifr < 4);Rex(iiping,ifr,idxex) = 0;end;
      clear idxex;
      idxex = find(theta > 47. & theta < 55);
      if(ifr == 5);Rex(iiping,ifr,idxex) = 0;end;
      clear idxex;      
      plot(theta_even,squeeze(R_even(ifr,:))./squeeze(Rex(iiping,ifr,:))','-k','LineWidth',2);
      plot(theta_even,squeeze(sd2(iiping,ifr,:)),'--k','LineWidth',1);
      plot(theta_even,squeeze(sd2p(iiping,ifr,:)),'--r','LineWidth',1);

      set(gca,'FontSize',14);
      set(gca,'Layer','top');
      set(gca,'YLim',[0 0.65],'XLim',[28 68]);
      text(50,0.58,[num2str(sd1(iiping,ifr),'%6.4f') '  ' num2str(freq(ifr),'%04i') ' Hz'],'FontSize',12)
      text(29,0.58,['Ping' num2str(iping,'%5d')],'FontSize',12)
      if(ifr == 1 | ifr == 4);
          yl=ylabel('Refl. coeff.');
      end;
      if(ifr ~= 1 & ifr ~=4);set(gca,'YTickLabel',[]);end;
      if(ifr < 4);set(gca,'XTickLabel',[]);
      else;xl=xlabel('Angle (deg.)');end;
   
      set(0, 'currentfigure',fig2);
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold off;box on;
      idxtmp = [idxlo(ifr):idxhi(ifr)];
      plot(theta_even,squeeze(sd2(iiping,ifr,:)),'--k','LineWidth',1);hold on;
      plot(theta_even,squeeze(sd2p(iiping,ifr,:)),'--r','LineWidth',1);
      plot(theta_even,squeeze(sd2f(iiping,ifr,:)),'--b','LineWidth',1);

      set(gca,'FontSize',14);
      set(gca,'Layer','top');
      set(gca,'YLim',[0 0.08],'XLim',[28 68]);
      text(50,0.58,[num2str(sd1(iiping,ifr),'%6.4f') '  ' num2str(freq(ifr),'%04i') ' Hz'],'FontSize',12)
      text(29,0.58,['Ping' num2str(iping,'%5d')],'FontSize',12)
      if(ifr == 1 | ifr == 4);
          yl=ylabel('Std. dev.');
      end;
      if(ifr ~= 1 & ifr ~=4);set(gca,'YTickLabel',[]);end;
      if(ifr < 4);set(gca,'XTickLabel',[]);
      else;xl=xlabel('Angle (deg.)');end;
      if(ifr == 6);legend('sd','sdp','sdf');end;
   end;
   if(i_save == 1);
     saveas(fig1,plotname1,'png');
     saveas(fig2,plotname2,'png');
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%
     %%  Save data
     %%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

     %F(1).cw = cw;
     %F(1).rw = rw;

     %F(1).freq = freq;
     %for ifr = 1:nf;
     %    F(ifr).dat = R_even(ifr,:);
     %    F(ifr).ang = theta_even;
     %end
     %F(1).Rex = Rex;

     %save(filenamea,'F');

     filebase = strrep(filenamea,'.mat','');
     fid = fopen('filebase.txt','w');
     filebase
     fprintf(fid,'%i\n',length(filebase));
     fprintf(fid,'%s',filebase);
     fclose(fid);

     save(filename,'z_t','-ascii');
     save(filename,'cw','-ascii','-append');
     save(filename,'rw','-ascii','-append');
     save(filename,'hmx','-ascii','-append');
     for ifr = 1:nf
       tmp = R_even(ifr,:);
       save(filename,'tmp','-ascii','-append');
     end
     tmp = theta_even';
     save(filename,'tmp','-ascii','-append');

     for ifr = 1:nf
       tmp = squeeze(Rex(iiping,ifr,:));
       tmp = tmp';
       save(filename,'tmp','-ascii','-append');
     end
     for ifr = 1:nf
%       tmp = squeeze(sd2p(iiping,ifr,:))';
       tmp = (squeeze(sd2f(iiping,ifr,:)))';
%       tmp = (squeeze(sd2(iiping,ifr,:)))';
       save(filename,'tmp','-ascii','-append');
     end
%     tmp = sd1p(iiping,:);
     tmp = sd1f(iiping,:);
%     tmp = sd1(iiping,:);
     save(filename,'tmp','-ascii','-append');
   end;
%   pause(.1);
%   k = waitforbuttonpress;
%   if k == 0
%    disp('Button click')
%   else;
%    disp('Key press: Discard')
%   end
   %close all;

end

fid = fopen('ping_list.txt','wt');
fprintf(fid,'%d\n',length(idxping));
for iping = 1:length(idxping);
  fprintf(fid,'%d\n',idxping(iping));
end;
fclose(fid);

NPING2 = length(idxping);
save sdpar_track.mat sd1 sd2 sd2p idxping;

figure();
nx = 1;
ny = 6;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
angidx = [1,3,30,32];
for ifr = 1:nf;
   subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold off;box on;
   isub = 1;
   for iang = 1:length(angidx);
       plot(idxping,sd2(1:NPING2,ifr,angidx(iang)),col(isub,:));hold on;
       isub = isub + 1;
   end;
   set(gca,'FontSize',14);
   set(gca,'Layer','top');
   set(gca,'YLim',[0 0.1],'XLim',[0 3500]);
   text(45,0.55,[num2str(freq(ifr),'%04i') ' Hz'],'FontSize',12)
   xlabel('Ping')
   if(ifr == 1 | ifr == 4);ylabel('Std. dev.');end;
   if(ifr ~= 1 & ifr ~= 4);set(gca,'YTickLabel',[]);end;
   grid on;
end;
legend('31.4','32.7','63.2','67.1');


nx = 6;
ny = 6;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
figure();
angidx = [1,3,30,32];
for iang = 1:NANG;
   subplot('Position',[loc(1,iang) loc(2,iang) spw sph]);
   hold off;box on;grid on;
   isub = 1;
   plot(idxping,sd2(1:NPING2,4,iang));hold on;
   set(gca,'FontSize',14);
   set(gca,'Layer','top');
   set(gca,'YLim',[0 0.1],'XLim',[0 3500]);
   text(50,0.01,[num2str(theta(iang),'%04i') ' Hz'],'FontSize',12)
   xlabel('Ping')
   if(iang == 1 | iang == 7);ylabel('Std. dev.');end;
   if(iang ~= 1 & iang ~= 7);set(gca,'YTickLabel',[]);end;
end;



figure();
nx = 1;
ny = 6;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
angidx = [1,3,30,32];
for ifr = 1:nf;
   subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold off;box on;
   isub = 1;
   for iang = 1:length(angidx);
       plot(idxping,sd2p(1:NPING2,ifr,angidx(iang)),col(isub,:));hold on;
       isub = isub + 1;
   end;
   set(gca,'FontSize',14);
   set(gca,'Layer','top');
   set(gca,'YLim',[0 0.1],'XLim',[0 3500]);
   text(45,0.55,[num2str(freq(ifr),'%04i') ' Hz'],'FontSize',12)
   xlabel('Ping')
   if(ifr == 1 | ifr == 4);ylabel('Std. dev.');end;
   if(ifr ~= 1 & ifr ~= 4);set(gca,'YTickLabel',[]);end;
   grid on;
end;
legend('31.4','32.7','63.2','67.1');


nx = 6;
ny = 6;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
figure();
angidx = [1,3,30,32];
for iang = 1:NANG;
   subplot('Position',[loc(1,iang) loc(2,iang) spw sph]);
   hold off;box on;grid on;
   isub = 1;
   plot(idxping,sd2p(1:NPING2,4,iang));hold on;
   set(gca,'FontSize',14);
   set(gca,'Layer','top');
   set(gca,'YLim',[0 0.1],'XLim',[0 3500]);
   text(50,0.01,[num2str(theta(iang),'%04i') ' Hz'],'FontSize',12)
   xlabel('Ping')
   if(iang == 1 | iang == 7);ylabel('Std. dev.');end;
   if(iang ~= 1 & iang ~= 7);set(gca,'YTickLabel',[]);end;
end;



%return;
