function [] = calc_critical_values();
%
% Calculates the critical values for the KS test 
% for a certain number of data points. The number of data 
% matters for small numbers.
%
% Number of data at the different frequencies:
%N = [52 58 63 80 106 121 130 131];
%N = [58 63 80 106 121 130 131 132 132];
%N = [62 76 84 87 91 92 92 92 92 92];
%N = [90 98 121 140 140 140 140 141 141 141];
N = [51];
%
% Number of realisations:
%
M = 15000;

for j = 1:length(N)
  k = j*2;		% plot parameter
%
% genrate M realisations of noise with length N(j) 
%
  G = randn(M,N(j));

  sd = std(G,0,2);
  %mu = mean(G,2);
  %sd = 1;
  mu = 0;


  for i = 1:N(j)
    y(i) = i * 1/N(j);
  end

  fprintf(1,'%4i',j);
  fprintf(1,' sorting and calc theoretical value...');
  for i = 1:M
    Gs(i,:) = sort(G(i,:));
    yth(i,:) = normcdf(Gs(i,:),mu,sd(i));
    diff(i,:) = y - yth(i,:);
  end
  fprintf(1,'done.\n');

  mxdiff = sqrt(N(j))*max(abs(diff),[],2);

  fprintf(1,'calculating percentiles...');
  percentiles = [80:0.5:100];
  z =  prctile(mxdiff,percentiles);
  fprintf(1,'done.\n');
  foo(:,k:k+1) = [percentiles; z]';
  clear Gs G y yth diff percentiles z;
end
save critical_values.mat foo;

%figure(1);
%hold on;
%stairs(Gs(1,:),y);
%stairs(Gs(1,:),yth(1,:),'r');

%figure(2);
%hist(mxdiff);
%s_mxdiff = sort(mxdiff);

%save bla s_mxdiff;

%
% 
%
Pr = 0;
t = 1.07
for i = 1:10000000
  Pr = Pr + (-1)^(i-1) * exp(-2*i^2*t^2);
end
Pr = 1 - 2 * Pr

return;
%
%This is the end my fiend
%EOF
