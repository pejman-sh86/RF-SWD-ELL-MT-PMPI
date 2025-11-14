function [] = pick_source()

outfile = 'source_s16_low.dat';
ibeg = 3500;
iend = 3900;
ibegt = floor((iend-ibeg)/40);
fact = 100;
i_offset = 4;
sd = 0.35;
rd = 122;
dt = 0.000041667;
f_cut = 500;
order = 8;
f_ny = 1/(2*dt);

load x;
fid = fopen(outfile,'w');

R = sqrt(x.r(i_offset)^2+(abs(sd-rd))^2);
source(:,1) = x.t_ax(ibeg:iend)-x.t_ax(ibeg);
source(:,2) = R .* x.t_s(i_offset,ibeg:iend);

%figure(1)
%plot(source(:,1),source(:,2));
%hold on;


taper = ones(size(source(:,1)));
x = [1:length(taper)-ibegt+1]';
taper(ibegt:end) = taper(ibegt:end) .* cos(pi/2.*x/x(end)).^fact;

figure(2)
plot(source(:,1),taper());
hold on;

source(:,2) = source(:,2).*taper;
figure(1)
hold on;
plot(source(:,1),source(:,2),'--k');

%
% Low-pass filter the source signal:
% (causes phase delay of a few points...)
%
butter(order,f_cut/f_ny)
[b,a] = butter(order,f_cut/f_ny);
Hd = dfilt.df2t(b,a);
phasedelay(Hd)
y = filter(Hd,source(:,2));
source(:,2) = y;

%source(22:33,2) = - source(11:22,2);

figure(1)
plot(source(:,1),source(:,2),'--r');

for i = 1:length(source)
  fprintf(fid,'%14.8f',source(i,1),source(i,2));fprintf(fid,'\n');
end

fclose(fid);

return;
