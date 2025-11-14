function [] = pick_source()

infile = 'source_13_raw.dat';
outfile = 'source_13_low.dat';
ibegt = 20;
fact = 400;
i_offset = 4;
dt = 0.000041667;
f_cut = 5500;
order = 8;
f_ny = 1/(2*dt);

source = load(infile);
fid = fopen(outfile,'w');

figure(1)
plot(source(:,1),source(:,2));
hold on;

taper = ones(size(source(:,1)));
x = [1:length(taper)-ibegt+1]';
ibegt
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
