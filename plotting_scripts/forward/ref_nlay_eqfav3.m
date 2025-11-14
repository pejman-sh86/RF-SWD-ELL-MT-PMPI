function [ref] = ref_nlay_eqfav3(thd,geo,freq)
%    ref = ref_nlay_fav3(1:90,geo,25:33);
%
% compute the plane wave reflection coeffiecnt for arbitrary number of isospeed layers
% averages over 100 Hz bands
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

spacing = [-50 -33.33 -16.67 0 16.67 33.33 50];

fr = spacing + freq;
rf = ref_nlay3(thd, geo,fr);
rf=rf.^2;

ref = sqrt(mean(rf(1:7,:), 1)); 


return;
