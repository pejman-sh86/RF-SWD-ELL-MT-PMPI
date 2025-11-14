function [] = resamp();

outfile  = strcat('x_s19_1_6_25_0lay.mat');
outfilea = strcat('x_s19_1_6_25_0lay.txt');

nang_re = 130;
load blme19_2_bl_ping_ida_nb.mat;

%idx = [5,6,8,10,13,17,21,27,33,43];
%idx = [5,6,8,10,13,17,21,27,33];
%idx = [6,8,10,13,17,21,27,33];
idx = [8,10,13,17,21,27,33];
%idx = [10,13,17,21,27,33];
nf = length(idx);
theta = xbl.ang(1:90)';
R = xbl.bl(idx,1:90);
R = R';

minlim = [0.50, 1450, 1450, 1.2, 1.3, 0.001, 0.001];
maxlim = [1.76, 1550, 1550, 1.6, 1.8, 1.500, 1.000];
mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.500, 0.5005];

%minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 1450, 1.3, 0.001];
%maxlim = [2.0, 1550, 1550, 1.6, 1.8, 1.000, 1550, 1.8, 1.000];
%mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.500, 1500, 1.5, 0.5005];

%minlim = [];
%maxlim = [];
%mtrue  = [];

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

F(1).freq = xbl.pref(idx,1)
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
   subplot(3,3,i);hold on;box on;
   plot(theta,R(:,i),'.k');
   plot(theta_even(:),R_even(:,i),'or');
   figure(2);
   subplot(3,3,i);hold on;box on;
   idx = find(Rex(:,i)==1);
   plot(theta,R(:,i),'.k');
   plot(theta_even(idx),R_even(idx,i),'or');
   figure(3);
   subplot(3,3,i);hold on;box on;
   idx = find(Rex(:,i)==1);
   plot(theta_even(idx),R_even(idx,i),'.r');


end;

return;
