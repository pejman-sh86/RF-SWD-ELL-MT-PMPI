function [] = xraw(seis_file);

RANGE = 0;
load(seis_file);
%load(x.mat);

t_lim = [0 1.2];
fact = 5e3
seis = double(DATA');
ntr = size(DATA,2)
ndt = size(DATA,1)
for i1=1:ndt
  t(i1) = DELTAT * i1;
end
r = double(RANGE);

i = find(t_lim(1) <= t & t <= t_lim(2));
x1 = seis(:, i);
for i1 = 1:length(r)
     x1(i1, :) = fact*x1(i1, :)+r(i1);
end
plot(x1, t(i), 'b')
axis([min(r), max(r), t(i(1)), t(i(end))])
xlabel('Offset (m)')
ylabel('Time re trigger (s)')

return;
