function [] = auv_new_processing(filename);

filebase = strrep(filename,'.mat','');

i_save = 0;
i_inter = 0;  %% 1 = resample onto even grid
i_frave = 1;  %% 1 = average frequencies 
i_pave  = 1;  %% 1 = apply ping averaging
nfave = 1;    %% total freq to ave is (2*nfave)+1
nskip = 1;   %% stride between pings
npave = 0;   %% stride between pings

nang_re = 31;
load(filename);

%% Frequencies:
idx = [7, 12, 18, 52, 64, 76];
freq = [975,1100,1250,2100,2400,2700];

idxlo = idx-nfave;
idxhi = idx+nfave;

nf = length(idx);
pstart = 1;

for iping = pstart:nskip:31;
   disp(iping);
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
   hmx = tb_win*1700./2.;

   if(i_pave == 1);
     Rptmp = RR_psCMP(:,:,iping-npave:iping+npave);
     Rp = mean(Rptmp,3);
     theta = mean(ANG_all_tilt(:,iping-npave:iping+npave),2)*360./(2.*pi);
     z_t = 2.*mean(alt(iping-npave:iping+npave));
   else;
     Rp = RR_psCMP(:,:,iping);
     theta = ANG_all_tilt(:,iping)*360./(2.*pi);
     z_t = 2.*alt(iping);
   end;

   if(i_frave == 1);
      for ifr=1:nf
        Rtmp = Rp(idxlo(ifr):idxhi(ifr),:);
        R(ifr,:) = sqrt(sum(Rtmp.*Rtmp,1)/((2*nfave)+1));
      end;
   else;
      R = Rp(idx,:);
   end;

   thdiff = diff(theta);
   if(i_inter == 1);
     theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];
     theta_even = theta_even';
     Rex = ones(nf,length(theta_even));
     for ifr = 1:nf
       R_even(ifr,:)   = interp1(theta,R(ifr,:),theta_even);
     end
   else
      theta_even = theta;
      R_even = R;
      Rex = ones(nf,length(theta_even));
   end;

   nx = 3;
   ny = ceil(nf/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   figw = 10.5;
   figh = 4.;
   fig1 = figure('visible','off');hold on; box on;
   fig2 = figure('visible','off');hold on; box on;
   set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   set(fig2,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   for ifr = 1:nf;
      set(0, 'currentfigure',fig1);
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      set(gca,'Layer','top');
      set(gca,'YLim',[0 0.65],'XLim',[28 68]);
      idx = [idxlo(ifr):idxhi(ifr)];
      idat = 1;
      for ip = 1:2*npave+1;
        for isub = 1:2*nfave+1;
          plot(theta_even,Rptmp(idx(isub),:,ip),':');
          res(idat,:) = Rptmp(idx(isub),:,ip);
          idat = idat + 1;
        end;
      end;
      sd(iping,ifr)=mean(std(res,1));
%      sdtmp=std(res,1);
      plot(theta_even,R_even(ifr,:),'-r','LineWidth',2);
      text(45,0.55,[num2str(sd(iping,ifr),'%6.4f') '  ' num2str(freq(ifr),'%04i') ' Hz'],'FontSize',12)
      xlabel('Angle (deg.)')
      if(ifr == 1);ylabel('Refl. coeff.');end;
      if(ifr ~= 1);set(gca,'YTickLabel',[]);end;
%      set(0, 'currentfigure',fig2);
%      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
%      plot(theta_even,sdtmp,'-r','LineWidth',2);
   end;
   if(i_save == 1);
     saveas(fig1,plotname1,'png');
%     saveas(fig2,plotname2,'png');
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

       tmp = Rex(ifr,:);
       save(filename,'tmp','-ascii','-append');

     end
     tmp = sd(iping,:);
     save(filename,'tmp','-ascii','-append');
   end;
   k = waitforbuttonpress;
   close all;

end
save sdpar_track.mat sd;
return;
