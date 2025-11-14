function [] = resamp();

outfile  = strcat('x_s07_1_4_25_4lay.mat');
outfilea = strcat('x_s07_1_4_25_4lay.txt');

nx = 5;
ny = 5;

nang_re = 80;
load blmn7_b_bl_ping_id1neg_nb_neg.mat;

idx   = [1,1,2,3,4,5,7,9,12,15,19,24,31,39,49,62,79, 99,125,159,199];
idxlo = [1,1,2,3,4,5,7,9,12,15,19,24,31,38,48,60,77, 97,123,156,194];
idxhi = [1,2,2,3,4,5,7,9,12,16,20,25,32,40,50,63,81,101,128,161,199];
freq = [100., 125, xbl.pref(idx(3:end),1)'];
freq

nf = length(idx);
theta = xbl.ang(1:52)';

%% Freq Ave:
for i = 1:length(idx);
   R(i,:) = mean(xbl.bl(idxlo(i):idxhi(i),1:52),1);
end;
R = R';
size(R)

R2 = xbl.bl(idx,1:52)';

%minlim = [0.20, 1450, 1.2, 0.001, ...
%                1450, 1.2, 0.001];
%maxlim = [2.60, 1750, 2.2, 1.000, ...
%                1750, 2.2, 1.000];
%minlim = [0.30, 1450, 1.2, 0.001, ...
%          0.30, 1450, 1.2, 0.001, ...
%                1450, 1.2, 0.001];
%maxlim = [2.50, 1750, 2.2, 1.000, ...
%          2.50, 1750, 2.2, 1.000, ...
%                1750, 2.2, 1.000];
%minlim = [0.30, 1450, 1.2, 0.001, ...
%          0.30, 1450, 1.2, 0.001, ...
%          0.30, 1450, 1.2, 0.001, ...
%                1450, 1.2, 0.001];
%maxlim = [2.50, 1750, 2.2, 1.000, ...
%          2.50, 1750, 2.2, 1.000, ...
%          2.50, 1750, 2.2, 1.000, ...
%                1750, 2.2, 1.000];
minlim = [0.20, 1450, 1.2, 0.001, ...
          0.20, 1450, 1.2, 0.001, ...
          0.20, 1450, 1.2, 0.001, ...
          0.20, 1450, 1.2, 0.001, ...
                1450, 1.2, 0.001];
maxlim = [2.50, 1750, 2.2, 1.000, ...
          2.50, 1750, 2.2, 1.000, ...
          2.50, 1750, 2.2, 1.000, ...
          2.50, 1750, 2.2, 1.000, ...
                1750, 2.2, 1.000];
mtrue  = (maxlim-minlim)/2.;

thdiff = diff(theta);
idx2 = find(thdiff == 0);
theta(idx2) = [];
R(idx2,:) = [];

theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];

Rex = ones(length(theta_even),nf);
for ifr = 1:nf

   Rex1 = -1*(isnan(R(:,ifr))-1);
   Rtmp = R(find(Rex1 == 1),ifr);
   Rtmp2 = R(:,ifr);
   Rtmp2(find(Rex1 == 0)) = 0;
   thetatmp = theta(find(Rex1 == 1));

   size(thetatmp)
   size(Rtmp)

   R_even(:,ifr)   = interp1(thetatmp,Rtmp,theta_even);
   R_even2(:,ifr)   = interp1(theta,Rtmp2,theta_even);
   Rex(find(R_even2(:,ifr) == 0),ifr) = 0;

   clear Rex1 Rtmp Rtmp2 thetatmp;
end

dang = diff(theta);
for i=2:length(theta)-1;
   if(dang(i)>1.2*dang(i-1) | dang(i) < 0.);
      idx3 = find(theta_even>theta(i) & theta_even < theta(i+1));
      Rex(idx3,:) = 0;
      clear idx3;
   end;
end;
R_even(find(Rex == 0)) = 0.;
R_even(isnan(R_even)) = 0.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = 1511.37;
F(1).rw = 1.029;

%F(1).freq = xbl.pref(idx,1)
F(1).freq = freq;
for ifr = 1:nf;
    F(ifr).dat = R_even(:,ifr);
    F(ifr).ang = theta_even;
end


F(1).minlim = minlim;
F(1).maxlim = maxlim;
F(1).mtrue  = mtrue;
F(1).Rex    = Rex;
F(1).nmod   = length(F(1).maxlim);

save(outfile,'F');
save(outfilea,'minlim','-ascii');
save(outfilea,'maxlim','-ascii','-append');
save(outfilea,'mtrue','-ascii','-append');
for i = 1:nf

    tmp = F(i).dat';
    save(outfilea,'tmp','-ascii','-append');

end
ang = F(1).ang(1:size(Rex,1));
save(outfilea,'ang','-ascii','-append');

for i = 1:nf

    tmp = F(1).Rex(:,i)';
    save(outfilea,'tmp','-ascii','-append');

end

for i = 1:nf

   figure(1);
   subplot(nx,ny,i);hold on;box on;
   plot(theta,R(:,i),'.k');
   plot(theta_even(:),R_even(:,i),'or');
   set(gca,'YLim',[0 20]);
   figure(2);
   subplot(nx,ny,i);hold on;box on;
   idx = find(Rex(:,i)==1);
   plot(theta,R(:,i),'.k');
   plot(theta_even(idx),R_even(idx,i),'or');
   set(gca,'YLim',[0 20]);
   figure(3);
   subplot(nx,ny,i);hold on;box on;
   idx = find(Rex(:,i)==1);
   plot(theta_even(idx),R_even(idx,i),'.r');
   set(gca,'YLim',[0 20]);
   figure(4);
   subplot(nx,ny,i);hold on;box on;
   plot(theta,R(:,i),'.--r');
   plot(theta,R2(:,i),'.k');
   set(gca,'YLim',[0 20]);

end;

return;
