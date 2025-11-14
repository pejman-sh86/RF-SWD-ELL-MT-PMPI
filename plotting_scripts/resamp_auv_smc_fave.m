function [] = resamp_auv_smc_fave(filename);

filebase = strrep(filename,'.mat','');

i_inter = 0;
i_frave = 1;
nfave = 1;   %% total freq to ave is (2*nfave)+1
nang_re = 31;
load(filename);

%% Frequencies:
%%1000 1200 2000
%idx = [5, 11, 19];
%freq = [925, 1075, 1275];
%%1000 1200 2000
%idx = [8, 18, 64];
%freq = [1000, 1250, 2400];
idx = [64, 82, 96];
freq = [2400,2850,3200];

idxlo = idx-nfave;
idxhi = idx+nfave;

nf = length(idx);
pst = 1;

for iping = 1:3433;
%for iping = 1:3:101;
%for iping = 1:101;
   iiping = (iping-1) + pst;
   filename = strcat('p',num2str(iiping,'%04i'),'_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.txt');
   filenamea= strcat('p',num2str(iiping,'%04i'),'_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.mat');
   plotname = strcat('p',num2str(iiping,'%04i'),'_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.png');

%   theta = ANG_all_tilt(:,iping)*360./(2.*pi);
%   R = Rsmooth(idx,:,iping);
%   theta = ANG(:,iping)*360./(2.*pi);
%   R = Rmap(idx,:,iping);
%   theta = ANG(:,iping)*360./(2.*pi);
%   R = RR(idx,:,iping);
   theta = ANG_all_tilt(:,iping)*360./(2.*pi);
   z_t = 2.*alt(iping);
%   z_t = 2.*ALT_t(iping);
   cw = 1512.3;
   rw = 1.029;
   tb_win = 0.0085;
   hmx = tb_win*1700./2.;

   if(i_frave == 1);

      for ifr=1:nf
         Rtmp = RR_psCMP(idxlo(ifr):idxhi(ifr),:,iping);
         R(ifr,:) = sqrt(sum(Rtmp.*Rtmp,1)/((2*nfave)+1));
      end;

   else;
      R = RR_psCMP(idx,:,iping);
   end;

   thdiff = diff(theta);
   if(i_inter == 1);
      theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];
   else
      theta_even = theta;
   end;

   theta_even = theta_even';
   Rex = ones(nf,length(theta_even));
   for ifr = 1:nf

      R_even(ifr,:)   = interp1(theta,R(ifr,:),theta_even);

   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %%  Save data
   %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   F(1).cw = cw;
   F(1).rw = rw;

   F(1).freq = freq;
   for ifr = 1:nf;
       F(ifr).dat = R_even(ifr,:);
       F(ifr).ang = theta_even;
   end
   F(1).Rex = Rex;

   save(filenamea,'F');

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
   for i = 1:nf

      tmp = F(i).dat;
      save(filename,'tmp','-ascii','-append');

   end
   save(filename,'theta_even','-ascii','-append');

   for i = 1:nf

       tmp = F(1).Rex(i,:);
       save(filename,'tmp','-ascii','-append');

   end

end
return;
