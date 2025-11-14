function [n1] = histc_tstar(y_in,E,bins,Tstar);

const = min(E);
nbin = length(bins)-1;
n1 = zeros(size(bins));

for i = 1:length(y_in)
   for j = 1:nbin
      if ((y_in(i) > bins(j)) & (y_in(i) <= bins(j+1)))
         n1(j) = n1(j) + exp(-(1.-1./Tstar)*(E(i)-const));
      end
   end
end

n1 = n1/(sum(n1)*(max(y_in)-min(y_in))/nbin)*(max(y_in)-min(y_in));

return;
