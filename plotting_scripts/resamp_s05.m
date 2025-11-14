function [] = resamp();

outfile  = strcat('x_s05_1_8_32_0lay.mat');
outfilea = strcat('x_s05_1_8_32_0lay.txt');

nang_re = 130;
load blje5_5_bl_ping_ida_nb_smooth.mat;
idx = [9 13 18 25 33 43 57];
nf = length(idx);
theta = thav_ping(1:244,end);
R = xavg_ping(1:244,idx);

minlim = [0.5, 1450., 1450., 1.2, 1.3, 0.001, 0.001];
maxlim = [2.0, 1550., 1550., 1.6, 1.8,   1.5, 1.];
mtrue  = [1, 1500, 1500, 1.4, 1.55, 0.5005, 0.5005];

thdiff = diff(theta);
idx2 = find(thdiff == 0);
theta(idx2) = [];
R(idx2,:) = [];

theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];

Rex = -1*(isnan(R)-1);

for ifr = 1:nf
   R(find(Rex(:,ifr) == 0),ifr) = 0;
end
clear Rex;
Rex = ones(length(theta_even),nf);
for ifr = 1:nf

   R_even(:,ifr)   = interp1(theta,R(:,ifr),theta_even);
   Rex(find(R_even(:,ifr) == 0),ifr) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = 1510.7;
F(1).rw = 1.029;

F(1).freq = freq(idx);
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
