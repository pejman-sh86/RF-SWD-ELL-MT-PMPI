function []=mfi2D_putnoise();

rep = load('rjmh_pecan_replica.dat');
z = rep(1,:);
rep = rep(2:end,:);
pct = 0.1;  %% Percent error on data

%%
%% Compute noise
%%
for idat = 1:size(rep,1);
  sigma = pct * max(abs(rep(idat,:)));
  noise(idat,:) = sigma*randn(1,size(rep,2));
end;
sim = rep+noise;

%%
%% Plot
%%
figure();
nx = 6;
ny = 5;
xim = 0.01;
yim = 0.05/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for idat = 1:size(rep,1);
  subplot('Position',[loc(1,idat) loc(2,idat) spw sph]);hold on;box on;
  plot(z,rep(idat,:),'b');
  plot(z,sim(idat,:),'--r');
%  set(gca,'XTickLabel',[],'YTickLabel',[])
end;

%%
%% Save as ascii
%%
save('rjmh_pecan_s_03_f_50_200_sim.txt','sim','-ascii')

return;
