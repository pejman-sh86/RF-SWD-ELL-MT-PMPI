function [ref] = ref_nlay_fav3(thd,geo,freq)
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

nave = 8;
otofc=...
[1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5 16 20 25 31.5 40 50 63 80 100 125 160 200 250 315  400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000];

for k=1:length(freq);
    j=find(freq(k)==otofc);
    otonum(k)=j; end

n=(otonum(1)*nave:(otonum(end)+1)*nave)-(nave/2);
fr=10.^(n/(nave*10));
%fr

%save bla geo;
%nave = 8;
%freq_low = freq.*2.^(-1/6);
%freq_up  = freq.*2.^(1/6);
%fr = [freq_low:(freq_up-freq_low)/nave:freq_up];
%f = 100:1:2400;
%idx = find(f <= freq_up & f >= freq_low);
%fr = f(idx);

rf = ref_nlay3(thd,geo,fr);

rf = abs(rf);

ref = sqrt(mean(rf.^2, 1));

return;
