%%
%% Performs ping averaging and computed sigma for AUV
%%
function [] = auv_pingave_compute_sigma_track();

%bands = [925, 1075, 1275];
%bands = [1000, 1250, 2400];
%bands = [950, 1300, 1850, 2400];
%bands = [2400,2850,3200];
bands = [1000,1250,2400,2850,3200];
nf = length(bands);
pstgl  = 1;
pendgl = 100;
%pendgl = 3433;
%pstgl  = 2735;
%pendgl = 2780;
NAVE   = 10;
NSKIP  = 5;
NPINGGL = pendgl-pstgl+1;
fileext  = '.txt';

%for i=1:NAVE/2:NPINGGL-1;
for i=1:NSKIP:NPINGGL-1;

   pst  = pstgl+(i-1);
   pend = pstgl+(i-1)+NAVE-1;
   filebase1 = 'p';
   filebase2 = strcat('_',num2str(bands(1),'%04i'),...
                      '_',num2str(bands(end),'%04i'));
   outfile1 = strcat(filebase1,num2str(pst,'%04i'),'_pave_',num2str(NAVE),...
                     filebase2,fileext);
   plotfile = strcat(filebase1,num2str(pst,'%04i'),'_pave_',num2str(NAVE),...
                     filebase2,'.png');
   outfile2 = strcat(filebase1,num2str(pst,'%04i'),filebase2,'_sigma',fileext);

   for i = 1:pend-pst+1;
      ping = pst+i-1;
      infile = strcat(filebase1,num2str(ping,'%04i'),filebase2,fileext);
      disp(infile);
      fp(:,:, i) = dlmread(infile);
   end

   dat = fp(5:5+nf-1,:,:);
   ang = squeeze(fp(5+nf,:,:));
   ang_cntr = mean(ang,2);
   NANG = size(ang,1);

   for i = 1:size(fp,3);
      for j = 1:nf;
         dat_interp(j,:,i) = interp1(ang(:,i),dat(j,:,i),ang_cntr,...
         'nearest','extrap');
      end;
   end;

   sigma = std(dat_interp,0,3);
   for i = 1:nf;
      for j=1:NANG
         if(j<=ceil(NAVE/2))
            std_filt(i,j) = sum(sigma(i,1:j+floor(NAVE/2)))/(NAVE-ceil(NAVE/2)+j);
         elseif(j>=NANG-floor(NAVE/2))
            std_filt(i,j) = sum(sigma(i,j-ceil(NAVE/2)+1:NANG))/(NANG-j+ceil(NAVE/2));
         else
            std_filt(i,j) = sum(sigma(i,j-floor(NAVE/2):j+floor(NAVE/2)))/NAVE;
         end;
      end
   end;

   %%
   %% Need to divide by sqrt(N) to properly represent sigma for the mean!
   %%
   sigma2 = sigma/sqrt(size(fp,3));
   std_filt2 = std_filt/sqrt(size(fp,3));
   %sigma2 = sigma%/sqrt(size(fp,3));
   %std_filt2 = std_filt%/sqrt(size(fp,3));

%   fig1=figure();hold on;box on;
%   fig2=figure();hold on;box on;
%   for i = 1:nf;
%      figure(fig1);
%      subplot(2,3,i);hold on;box on;
%      text(55,0.055,[num2str(bands(i)) ' Hz'],'FontSize',12)
%      set(gca,'YLim',[0 0.06],'XLim',[28 68]);
%      plot(ang_cntr,sigma(i,:),'xk');
%      plot(ang_cntr,std_filt(i,:),'--b');
%      xlabel('Angle (deg.)')
%      ylabel('Std. dev.')
%
%      figure(fig2);
%      subplot(2,3,i);hold on;box on;
%      text(55,0.028,[num2str(bands(i)) ' Hz'],'FontSize',12)
%      set(gca,'YLim',[0 0.03],'XLim',[28 68]);
%      plot(ang_cntr,sigma2(i,:),'xk');
%      plot(ang_cntr,std_filt2(i,:),'--b');
%      xlabel('Angle (deg.)')
%      ylabel('Std. dev. of mean')
%   end;

   nx = 4;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   figw = 10.5;
   figh = 4.;
   fig1 = figure('visible','off');hold on; box on;
   set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   for i = 1:nf;
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      set(gca,'Layer','top');
      set(gca,'YLim',[0 0.65],'XLim',[28 68]);
      text(55,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12)
%      for j = 1:size(fp,3);
%%         plot(ang(:,j),dat(i,:,j),'--b');
%         plot(ang_cntr,dat_interp(i,:,j),'-b','LineWidth',0.5);
%      end;
%      dat_ave(i,:) = mean(dat(i,:,:),3);
      dat_ave(i,:) = median(dat(i,:,:),3);
      plot(ang_cntr,dat_ave(i,:),'-r','LineWidth',2);
%      plot(ang(:,6),dat(i,:,6),'-g','LineWidth',2);
%      plot(ang(:,1),dat(i,:,1),'-k','LineWidth',2);
      xlabel('Angle (deg.)')
      if(i==1);ylabel('Refl. coeff.');end;
      if(i~=1);set(gca,'YTickLabel',[]);end;
   end;
   sd = mean(sigma,2);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% Save data
   %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   outfile1
   z_t = mean(fp(1,1,:));
   cw  = mean(fp(2,1,:));
   rw  = mean(fp(3,1,:));
   hmx = max(fp(4,1,:));
   save(outfile1,'z_t','-ascii');
   save(outfile1,'cw','-ascii','-append');
   save(outfile1,'rw','-ascii','-append');
   save(outfile1,'hmx','-ascii','-append');
   for i = 1:nf

      tmp = dat_ave(i,:);
      save(outfile1,'tmp','-ascii','-append');
      clear tmp;

   end
   tmp = ang_cntr';
   save(outfile1,'tmp','-ascii','-append');
   clear tmp;

   for i = 1:nf

      tmp = fp(nf+4+1+i,:,1);
      save(outfile1,'tmp','-ascii','-append');
      clear tmp;

   end
   tmp = 5*mean(std_filt2,2)';
   save(outfile1,'tmp','-ascii','-append');
   clear tmp;

   tmp = std_filt2(1,:);
   save(outfile2,'tmp','-ascii');
   for i = 2:nf
      tmp = std_filt2(i,:);
      save(outfile2,'tmp','-ascii','-append');
      clear tmp;
   end
   saveas(fig1,plotfile,'png');
   close(fig1);

end;
return;
