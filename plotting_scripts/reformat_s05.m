function [] = reformat();

outfile  = strcat('x_s05_1_8_32_0lay.mat');
outfilea = strcat('x_s05_1_8_32_0lay.txt');
outfile2  = strcat('x_s05_1_8_32_0lay2.mat');

%%
%% 50 Hz bandwidth:
%%
load blje5_5_bl_ping_ida_nb_smooth.mat;
idx = [9 13 18 25 33 43 57];

%%
%% 75 Hz bandwidth:
%%
%load blje5_5_sl_idaf_nb_smooth.mat;
%idx = [10,13,17,21,27,33,43,53,67];
%idx = [10,13,17,21,27,33,43,53];
%idx = [5,6,8,10,13,17,21,27,33,43];


nf = length(idx);

%%
%% 50 Hz data:
%%
theta = thav_ping(1:244,idx);
R = xavg_ping(1:244,idx);

%%
%% 75 Hz data:
%%
%theta = thav1(1:244,idx);
%R = xavg1(1:244,idx);

minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 0.001];
maxlim = [1.9, 1550, 1550, 1.6, 1.8, 1.500, 1.000];

%minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 1450, 1.2, 0.001];
%maxlim = [2.0, 1550, 1550, 1.6, 1.8, 1.000, 1550, 1.8, 1.000];

%minlim = [0.05   1450.000   1.200   0.001 ...
%          0.05   1450.000   1.200   0.001 ...
%          0.05   1450.000   1.400   0.001 ...
%                 1450.000   1.400   0.001];
%maxlim = [1.00   1520.000   1.700   1.000 ...
%          1.00   1520.000   1.800   1.000 ...
%          1.00   1520.000   1.800   1.000 ...
%                 1520.000   1.900   1.000];

%minlim = [0.05,  1450.000,  1.200, 0.001,...
%          0.05,  1450.000,  1.200, 0.001,...
%          0.05,  1450.000,  1.400, 0.001,...
%          0.05,  1450.000,  1.400, 0.001,...
%                 1450.000,  1.400, 0.001];
%maxlim = [0.75,  1520.000,  1.600, 1.000,...
%          0.75,  1520.000,  1.700, 1.000,...
%          0.75,  1520.000,  1.800, 1.000,...
%          0.75,  1520.000,  1.800, 1.000,...
%                 1520.000,  1.900, 1.000];

%minlim = [0.05,  1450.000,  1.200, 0.001,...
%          0.05,  1450.000,  1.200, 0.001,...
%          0.05,  1450.000,  1.400, 0.001,...
%          0.05,  1450.000,  1.400, 0.001,...
%          0.05,  1450.000,  1.400, 0.001,...
%                 1450.000,  1.400, 0.001];
%maxlim = [0.75,  1520.000,  1.600, 1.000,...
%          0.75,  1520.000,  1.700, 1.000,...
%          0.75,  1520.000,  1.800, 1.000,...
%          0.75,  1520.000,  1.800, 1.000,...
%          0.75,  1520.000,  1.800, 1.000,...
%                 1520.000,  1.900, 1.000];

mtrue  = minlim+(maxlim-minlim)/2.;

for ifr = 1:length(idx);

  %%
  %% Find NaN data
  %%
  Rex = -1*(isnan(R(:,ifr))-1);

  R_tmp = R(find(Rex == 1),ifr);
  theta_tmp = theta(find(Rex == 1),ifr);

  F(ifr).dat = R_tmp(1:1:end);
  F(ifr).ang = theta_tmp(1:1:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = 1510.7;
F(1).rw = 1.029;

  
F(1).freq = freq(idx);

F(1).minlim = minlim;
F(1).maxlim = maxlim;
F(1).mtrue  = mtrue;
F(1).nmod   = length(F(1).maxlim);

F2 = F;

save(outfile,'F');

for i = 1:nf
    
    F2(i).dat = F(i).dat(1:1:end);
    F2(i).ang = F(i).ang(1:1:end);
    tmp(i) = length(F(i).dat);

end
save(outfile2,'F2');
save(outfilea,'tmp','-ascii');
save(outfilea,'minlim','-ascii','-append');
save(outfilea,'maxlim','-ascii','-append');
save(outfilea,'mtrue','-ascii','-append');
for i = 1:nf

    tmp = F(i).dat';
    save(outfilea,'tmp','-ascii','-append');

end

for i = 1:nf

    tmp = F(i).ang';
    save(outfilea,'tmp','-ascii','-append');

end

figure(1);
for i=1:length(idx);
  subplot(ceil(length(idx)/4),4,i);hold on;box on;
  plot(F(i).ang,F(i).dat,'--k');
  plot(F2(i).ang,F2(i).dat,':k');
%  plot(thav_ping(:,idx(i)),xavg_ping(:,idx(i)),'or');

  text(70,30,num2str(F(1).freq(i)));
  set(gca,'XLim',[0 90],'YLim',[0 40]);
end;

return;
