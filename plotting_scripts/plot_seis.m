function []=plot_seis(seis_file);

RANGE = 0;
load(seis_file);

seis = double(DATA);
ntr = size(DATA,2)
ndt = size(DATA,1)
for i=1:ndt
  dt(i) = DELTAT * i;
end
r = double(RANGE);

figure(1);
hold on; box on; grid on;
s1max = 0.5*max(abs(seis(:,1)));
for i=1:100%ntr
  plot(seis(:,i)+i*s1max,dt);
end
set(gca,'YDir','reverse')

return;
