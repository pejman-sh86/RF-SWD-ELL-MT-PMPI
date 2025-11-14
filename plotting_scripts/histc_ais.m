function [n1,lim] = histc_ais(y_in,logwt,bins);

%
% Computes histrograms corrected for annealed importance sampling. 
%

wt = exp(logwt);
nbin = length(bins)-1;
n1 = zeros(size(bins));

for i = 1:length(y_in)
   for j = 1:nbin
      if ((y_in(i) > bins(j)) & (y_in(i) <= bins(j+1)))
         n1(j) = n1(j) + abs(wt(i));
      end
   end
end

n1 = n1/(sum(n1)*(max(y_in)-min(y_in))/nbin);

return;
