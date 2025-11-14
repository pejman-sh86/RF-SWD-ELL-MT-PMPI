freq = [315 400 500 630 800 1000 1250 1600 2000];

for i = 1:length(freq)
 
  band(i,1) = freq(i)/(2^(1./6.));
  band(i,2) = freq(i);
  band(i,3) = freq(i)*(2^(1./6.));
  band(i,4) = band(i,3)-band(i,1);

end

band

