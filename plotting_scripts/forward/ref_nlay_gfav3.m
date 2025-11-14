function [ref] = ref_nlay_gfav3(thd,geo,freq,frbw)
%    ref = ref_nlay_fav3(1:90,geo,25:33);
%
% compute the plane wave reflection coeffiecnt for arbitrary number of isospeed layers
% averages over a 1/3 octave band
%
% thd are grazing angles in degrees
%
% geo contains environmental data, e.g., 
% geo=[NaN       1511          0            1.03;
%     1.7         1500          0.015         1.45;
%     2.8         1555          0.1          1.7 ;
%     0.6         1600          0.1            2;
%     NaN         1565          0.01          1.9];
%
% freq are the frequencies (Hz)
%
% Charles W. Holland - Jan 2004
% Jan Dettmer Version Aug 2005

% assume frequencies are spaced in 1/3 octave centers
% divide frequencies into logarithmic spacing; 8 should suffice to get smooth average

nave = 6;
flo = freq - freq*frbw;
fhi = freq + freq*frbw;
fr = [flo:(fhi-flo)/(nave-1):fhi];
rf = ref_nlay3(thd,geo,fr);
rf = abs(rf);
nang = length(rf);

%
% Gaussian average (Harrison and Harrison)
%
df = mean(diff(fr));
f0 = freq;
for iang = 1:nang

    r_ave(iang) = sum(rf(:,iang) .* ...
                  exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
                  sum(exp(-(fr'-f0).^2/(frbw*f0)^2)*df);

end;

ref = r_ave;

return;
