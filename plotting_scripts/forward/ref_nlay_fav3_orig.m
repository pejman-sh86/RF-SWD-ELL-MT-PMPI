function [ref] = ref_nlay_fav3(thd,geo,freq)
%    ref = ref_nlay_fav3(1:90,geo,25:33);
%
% compute the plane wave reflection coeffiecnt for arbitrary number of isospeed layers
% averages over a 1/3 octave band
%
% thd are grazing angles in degrees
%
% geo contains environmental data, e.g., 
% geo=[NaN       1511          0            1.03
%     1.7         1500          0.015         1.45
%     2.8         1555          0.1          1.7 
%     0.6         1600          0.1            2
%     NaN         1565          0.01          1.9];
%
% freq are the frequencies (Hz)
%
% Charles W. Holland - Jan 2004

% assume frequencies are spaced in 1/3 octave centers
% divide frequencies into logarithmic spacing; 8 should suffice to get smooth average

otofc=...
[100 125 160 200 250 315  400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000];


for k=1:length(freq);
    j=find(freq(k)==otofc);
    otonum(k)=j+19;
end

n=(otonum(1)*8:(otonum(end)+1)*8)-4
fr=10.^(n/80);
rf = ref_nlay3(thd, geo,fr);
rf=rf.^2;

for k=1:length(otonum)
 kp=(k-1)*8+1;
    ref(k,:)= sqrt(mean(rf(kp:kp+8,:), 1)); 
end
