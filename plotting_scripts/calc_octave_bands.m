function [fr] = calc_octave_bands(freq)
nave = 8;
otofc=...
[1.25 1.6 2 2.5 3.15 4 5 6.3 8 10 12.5 16 20 25 31.5 40 50 63 80 100 125 160 200 250 315  400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000];

for k=1:length(freq);
    j=find(freq(k)==otofc);
    j
    otonum(k)=j; 
end

n=(otonum(1)*nave:(otonum(end)+1)*nave)-(nave/2)
fr=10.^(n/(nave*10));

return;

