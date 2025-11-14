%%
%% Performs ping averaging and computed sigma for AUV
%%
function [] = auv_pingave_compute_sigma();

pst  = 175;
pend = 185;
NAVE = 9;
nf = 3;
bands = [1012, 1212, 2362]
filebase = 'R_20090522_p0175_1000_2350';
filebase1 = 'R_20090522_p0';
filebase2 = '_1000_2350';
fileext  = '.txt';
outfile1 = strcat(filebase1,num2str(pst),'_p0',num2str(pend),filebase2,fileext);
outfile2 = strcat(filebase1,num2str(pst),'_p0',num2str(pend),filebase2,'_sigma',fileext);

for i = 1:pend-pst+1;
   ping = pst+i-1;
   infile = strcat(filebase1,num2str(ping),filebase2,fileext);
   disp(infile)
   fp(:,:, i) = dlmread(infile);
end

dat = fp(5:5+nf-1,:,:);
ang = squeeze(fp(5+nf,:,:));
idx_cntr = round(size(fp,3)/2)
%ang_cntr = ang(:,idx_cntr);
ang_cntr = mean(ang,2);
NANG = size(ang,1);

for i = 1:size(fp,3);
   for j = 1:nf;
      dat_interp(j,:,i) = interp1(ang(:,i),dat(j,:,i),ang_cntr,...
      'nearest','extrap');
   end;
end;

%figure();hold on;box on;
%for i = 1:size(fp,3);
%   plot(fp(11,:,i),'xk');
%end;

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

fig1=figure();hold on;box on;
fig2=figure();hold on;box on;
for i = 1:nf;
   figure(fig1);
   subplot(2,3,i);hold on;box on;
   text(55,0.055,[num2str(bands(i)) ' Hz'],'FontSize',12)
   set(gca,'YLim',[0 0.06],'XLim',[28 68]);
   plot(ang_cntr,sigma(i,:),'xk');
   plot(ang_cntr,std_filt(i,:),'--b');
   xlabel('Angle (deg.)')
   ylabel('Std. dev.')

   figure(fig2);
   subplot(2,3,i);hold on;box on;
   text(55,0.028,[num2str(bands(i)) ' Hz'],'FontSize',12)
   set(gca,'YLim',[0 0.03],'XLim',[28 68]);
   plot(ang_cntr,sigma2(i,:),'xk');
   plot(ang_cntr,std_filt2(i,:),'--b');
   xlabel('Angle (deg.)')
   ylabel('Std. dev. of mean')
end;


figure();hold on;box on;
for i = 1:nf;
   subplot(2,3,i);hold on;box on;
   set(gca,'Layer','top');
   set(gca,'YLim',[0 0.65],'XLim',[28 68]);
   text(55,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12)
   for j = 1:size(fp,3);
%      plot(ang(:,j),dat(i,:,j),'--b');
      plot(ang_cntr,dat_interp(i,:,j),'-b','LineWidth',0.5);
   end;
%   dat_ave(i,:) = mean(dat(i,:,:),3);
   dat_ave(i,:) = median(dat(i,:,:),3);
   plot(ang_cntr,dat_ave(i,:),'-r','LineWidth',2);
%   plot(ang(:,6),dat(i,:,6),'-g','LineWidth',2);
%   plot(ang(:,1),dat(i,:,1),'-k','LineWidth',2);
   xlabel('Angle (deg.)')
   ylabel('Refl. coeff.')
end;
sd = mean(sigma,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outfile1
z_t = mean(fp(1,1,:))
cw  = mean(fp(2,1,:))
rw  = mean(fp(3,1,:))
hmx = max(fp(4,1,:))
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

tmp = std_filt2(1,:);
save(outfile2,'tmp','-ascii');
for i = 2:nf

   tmp = std_filt2(i,:);
   save(outfile2,'tmp','-ascii','-append');
   clear tmp;

end

return;
