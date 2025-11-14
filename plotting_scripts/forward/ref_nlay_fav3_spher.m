function [spher_r_num] = ref_nlay_fav3(thd,geo,freq)
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
% Jan Dettmer Version Aug 2005

% assume frequencies are spaced in 1/3 octave centers
% divide frequencies into logarithmic spacing; 8 should suffice to get smooth average

nave = 8;

d3 = .01;
thdm = thd-d3;
thdp = thd+d3;

freq_low = freq.*2.^(-1/6)
freq_up  = freq.*2.^(1/6)
fr = [freq_low:(freq_up-freq_low)/nave:freq_up];

rf=ref_nlay3(thd,geo,fr);
rfm=ref_nlay3(thdm,geo,fr);
rfp=ref_nlay3(thdp,geo,fr);

%rf = ref_nlay3(thd, geo,fr);
rf = rf.^2;
ref = sqrt(mean(rf, 1));
rfp = rfp.^2;
refp = sqrt(mean(rfp, 1));
rfm = rfm.^2;
refm = sqrt(mean(rfm, 1));

[thd_deg,spher_r_num]=spherical_ref_Brek_num_b(thd,...
                     thdm,thdp,1511,ref,refm,refp,freq,140);

spher_r_num = -20*log10(abs(spher_r_num));

return;
