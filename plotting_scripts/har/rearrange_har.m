%
% Rearrange Harrison data using
% parameters.dat
%
function [] = rearrange_har();

datafile = 'har_data_750-1265.mat';
load V0Th0.mat;

%i = 9;
%for ifreq = 1:12
%i = 24;
i = 17;
for ifreq = 1:12
  B(ifreq).dat = V0(i,1:2:91);
  fr(ifreq) = freq(i);
  i = i+1;
end

save(datafile,'B','fr');

return;
