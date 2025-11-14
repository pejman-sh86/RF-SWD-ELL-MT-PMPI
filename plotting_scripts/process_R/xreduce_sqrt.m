function xreduce_sqrt(seis_file)

RANGE = 0;
load(seis_file);

c=1525.;
tsq_lim=[0.00 0.10];
fact=2e2
seis = double(DATA');
ntr = size(DATA,2)
ndt = size(DATA,1)
for i=1:ndt
  t(i) = DELTAT * i;
end
r = double(RANGE);

%clf
for i = 1:length(r)-1
  dr(i)=r(i+1)-r(i);
end
for i = 1:length(r)
     tsq = sqrt(t.^2-r(i)^2/c^2);
     i1 = find(tsq_lim(1) <= tsq & tsq <= tsq_lim(2));
     line('Color', 'b', ...
          'XData', fact*seis(i, i1)+r(i), ...
          'YData', tsq(i1), 'LineWidth',.1)
end

xlabel('Range (m)'); ylabel('Reduced Time (s)');axis([min(r) max(r) tsq_lim]); box on 
