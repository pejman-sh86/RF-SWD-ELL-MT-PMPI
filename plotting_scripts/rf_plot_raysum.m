function []=rf_plot_raysum(trfile);

%trfile = 'sample.tr';

nmx = 200;

fid = fopen(trfile);
tline = fgets(fid);
tline = fgets(fid);
tmp = sscanf(tline, '%d');

ntr = tmp(1);
nsmp = tmp(2);

ismp = 1;
itr = 0;
while ~feof(fid);
  tline = fgets(fid);
  if(tline(1)=='#');
    if(ismp > 0)itr = itr + 1; end;
    ismp = 0;
    continue;
  end;
  ismp = ismp + 1;
  tr(itr,ismp,1:3) = sscanf(tline, '%f');
 
end;

R = tr(:,:,1);
T = tr(:,:,2);
V = tr(:,:,3);
save tmp.mat V T R;

Vmaxamp = max(max(V));
Rmaxamp = max(max(R));
Tmaxamp = max(max(T));
V = V/Vmaxamp;
R = R/Rmaxamp;
T = T/Tmaxamp;

figure();
hold on;box on;
for irf=1:ntr;
  plot(V(irf,:) + (irf/2),'-k');
end;
set(gca,'XLim',[0,nmx],'YLim',[0 ntr/2+Vmaxamp])

figure();
hold on;box on;
for irf=1:ntr;
  plot(R(irf,:) + (irf/2),'-k');
end;
set(gca,'XLim',[0,nmx],'YLim',[0 ntr/2+Vmaxamp])

figure();
hold on;box on;
for irf=1:ntr;
  plot(T(irf,:) + (irf/2),'-k');
end;
set(gca,'XLim',[0,nmx],'YLim',[0 ntr/2+Vmaxamp])

return;
