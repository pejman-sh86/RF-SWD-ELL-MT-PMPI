function [n1,lim] = hist_ais(y_in,logwt,nbin);

%
% Computes histrograms corrected for annealed importance sampling. 
%

wt = exp(logwt);
bins = [min(y_in):(max(y_in)-min(y_in))/nbin:max(y_in)];
lim = bins(2:end)-(max(y_in)-min(y_in))/nbin/2;
dx = mean(diff(lim));
n1 = zeros(size(lim));

for i = 1:length(y_in)
   for j = 1:nbin
      if ((y_in(i) > bins(j)) & (y_in(i) <= bins(j+1)))
         n1(j) = n1(j) + abs(wt(i));
      end
   end
end

n1 = n1/(sum(n1)*dx);
%n1 = n1/(sum(n1)*(max(y_in)-min(y_in))/nbin)*(max(y_in)-min(y_in));

return;
