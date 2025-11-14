function [] = resamp();

outfile  = strcat('x_s20_1_4_32_0lay.mat');
outfilea = strcat('x_s20_1_4_32_0lay.txt');

nang_re = 130;
load blmn20_t_bl_ping_idaf_NB_neg_smooth1.mat;

idx  = [5,6,8,10,13,17,21,27,33,43];
idx1 = [5,6,8,10,12,15,19,24,30,38];
idx2 = [5,7,9,11,14,18,23,29,37,47];

nf = length(idx);
theta = thav1(1:139,13);
theta
for i = 1:length(idx);
  R(:,i) = mean(xavg1(1:139,idx1(i):idx2(i)),2);
end;

minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 0.001];
maxlim = [2.0, 1550, 1550, 1.6, 1.8, 1.000, 1.000];
mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.500, 0.5005];
%minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 1450, 1.3, 0.001];
%maxlim = [2.0, 1550, 1550, 1.6, 1.8, 1.000, 1550, 1.8, 1.000];
%mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.500, 1500, 1.5, 0.5005];

%minlim = [];
%maxlim = [];
%mtrue  = [];

thdiff = diff(theta);
idx3 = find(thdiff == 0);
theta(idx3) = [];
R(idx3,:) = [];

theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];

Rex = -1*(isnan(R)-1);

for ifr = 1:nf
   R(find(Rex(:,ifr) == 0),ifr) = 0;
end
clear Rex;
Rex = ones(length(theta_even),nf);
for ifr = 1:nf

   ifr
   R_even(:,ifr)   = interp1(theta,R(:,ifr),theta_even);
   Rex(find(R_even(:,ifr) == 0),ifr) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = 1511.37;
F(1).rw = 1.029;

F(1).freq = freq(idx)
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

return;
