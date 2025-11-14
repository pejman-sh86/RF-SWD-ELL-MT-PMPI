function [] = resamp_auv_smc(filename);

filebase = strrep(filename,'.mat','');

i_inter = 0;
i_frave = 3;
nang_re = 31;
load(filename);

%% Frequencies:
%%1000 1200 2000
%idx = [8, 18, 48];
%freq = [1000, 1250, 2000];
%%1000 1200 2000
idx = [8, 18, 64];
freq = [1000, 1250, 2400];
%idx = [6, 20, 42, 64];
%freq = [950, 1300, 1850, 2400];
%%1000 2000
%idx = [8, 48];
%freq = [1000, 2000];
%%2000 3200
%idx = [48, 64, 80, 96];
%freq = [2000 2400 2800 3200];

nf = length(idx);
pst = 1;

for iping = 1:3433;
%for iping = 1:3:101;
%for iping = 1:101;
   iiping = (iping-1) + pst;
   filename = strcat('p',num2str(iiping,'%04i'),'_1000_2400.txt');
   filenamea= strcat('p',num2str(iiping,'%04i'),'_1000_2400.mat');
   plotname = strcat('p',num2str(iiping,'%04i'),'_1000_2400.png');

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
   hmx = tb_win*1750./2.;

   if(i_frave == 1);

      for ifr=1:nf
         idxup = idx(ifr) + 1;
         idxdn = idx(ifr) - 1;
         Rtmp = RR_psCMP(idxdn:idxup,:,iping);
         R(ifr,:) = sqrt(sum(Rtmp.*Rtmp,1)/nf);
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
