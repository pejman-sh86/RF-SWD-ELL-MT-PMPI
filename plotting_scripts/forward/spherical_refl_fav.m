function [Rs, Rp] = spherical_refl_fav(geo,z_t,freq,Th)
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

nave = 10;
frbw = 1/20;
%
%freq_low = freq.*2.^(-1/6);
%freq_up  = freq.*2.^(1/6);
%fr = [freq_low:(freq_up-freq_low)/nave:freq_up];
%
%otofc=...
%[1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5 16 20 25 31.5 40 50 63 80 100 125 160 200 250 315  400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000];

%for k=1:length(freq);
%    j=find(freq(k)==otofc);
%    otonum(k)=j; 
%end;
%ifr = 1
%for k=1:length(freq);
%    frlo = freq(1)-10.;
%    frup = freq(1)+10.;
%    fr = freq;
%    fr(ifr) = freq(k)-freq(k)/40.;
%    fr(ifr+1) = freq(k);
%    fr(ifr+2) = freq(k)+freq(k)/40.;
%    ifr = ifr+3;
%end;

%n=(otonum(1)*nave:(otonum(end)+1)*nave)-(nave/2);
%fr=10.^(n/(nave*10))

frlo = freq(1)-freq(1)*frbw;
frup = freq(1)+freq(1)*frbw;
fr = [frlo:(frup-frlo)/(nave-1):frup];

[Rsmf, Rpmf] = spherical_refl(geo,z_t,fr,Th);

Rpmf = abs(Rpmf');
Rsmf = abs(Rsmf');
%Rp = sqrt(mean(Rpmf.^2, 2));
%Rs = sqrt(mean(Rsmf.^2, 2));

%
% Gaussian average (Harrison and Harrison)
%
nang = length(Th);
df = mean(diff(fr));
f0 = freq;
for iang = 1:nang
    Rs(iang) = sum(Rsmf(:,iang) .* ...
               exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
               sum(exp(-(fr'-f0).^2.0/(frbw*f0)^2.0)*df);
    Rp(iang) = sum(Rpmf(:,iang) .* ...
               exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
               sum(exp(-(fr'-f0).^2/(frbw*f0)^2)*df);
end;

return;

