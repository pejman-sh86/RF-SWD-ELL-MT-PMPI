%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [logE,meanL,maxL,occ,BIC,AIC,C] = evidence(filebase,len,len2);

isd = 0;
ppd   = strcat(filebase,'_sample.mat');
sdfile  = strcat(filebase(1:len),'_exstd.txt');
data  = strcat(filebase,'.mat');
covd  = strcat(filebase(1:len2),'cov_cov.mat');
%covd  = strcat(filebase,'_cov.mat');

%N  = 6;  % Number of data per frequency
%Nf = 1;	% Number of data per frequency
ppd2 = ppd;
data2 = data;
covd2 = covd;
ii=1;
load(ppd2);
if(isd == 0)
   load(covd2);
   D = F;
else
   sd = load(sdfile);
%   sd = 0.05*ones(1,N);
end

load(data2);
C(ii).ppd = ppd2;

m = A(:,2:end)';
phi = A(:,1)';

Q  = length(phi);	% Size of sample
M  = size(m,1);		% Number of model parameters
N  = size(F(1).ang,2);  % Number of data per frequency
Nf = length(F(1).freq);	% Number of data per frequency


%
% Values from sample:
% -----------------------------------------------------------------------------
meanL(ii) = -mean(A(:,1));
maxL(ii) = -min(A(:,1));

mu = mean(m')';
Cm = cov(m');
invCm = inv(Cm);
detCm = det(Cm);

%
% Volume of the model space:
% -----------------------------------------------------------------------------
V(ii) = prod(F(1).maxlim-F(1).minlim);
%V(ii) = 1.;
logV(ii) = log(V(ii));

%
% Normalization for f(m_i)
% -----------------------------------------------------------------------------
fnorm(ii)  = 0;
fnorm2(ii) = 0;
for i = 1:Nf

  if(isd == 0)
     [L,U] = lu(D(i).csave);
     logdetCd(i)  = sum(log(abs(diag(U))));
  else
     [L,U] = lu(sd(i)*sd(i)*eye(N));
     logdetCd(i)  = sum(log(abs(diag(U))));
  end
  idx = find(F(1).Rex(1:N,i) == 1);
  N2(i) = length(idx);
%  N2(i) = N;

  fnorm(ii)  = fnorm(ii)  + 0.5*N2(i)*log(2*pi) + 0.5*logdetCd(i);
  fnorm2(ii) = fnorm2(ii) + 0.5*N2(i)*log(2*pi) + 0.5*logdetCd(i);

end


fnorm3(ii) = -fnorm(ii);
fnorm(ii) = -logV(ii) - fnorm(ii);
fnorm2(ii) = - fnorm2(ii);

%
% Normalization for g(m_i)
% -----------------------------------------------------------------------------
gnorm(ii) = -log( (2*pi)^(M/2) ) - log( sqrt(detCm) );

%
% Compute gamma (exponent of g(m_i))
% -----------------------------------------------------------------------------
for i = 1:Q;

   gamma(i) = 0.5*(m(:,i)-mu)'*invCm*(m(:,i)-mu);

end;

%
% Compute log of evidence 
% -----------------------------------------------------------------------------
z = -1e308;
for i = 1:Q;
   z = logplus(z,-gamma(i)+phi(i));
end;
logE(ii) = log(Q) - gnorm(ii) + fnorm(ii) - z

%
% Ockham Factor (for comparison)
% -----------------------------------------------------------------------------
occ(ii) = (fnorm3(ii) + maxL(ii)) + log(1/V(ii)*sqrt(detCm)*(2*pi)^(M/2));

%
% Bayesian Information Criterion (for comparison)
% -----------------------------------------------------------------------------
BIC(ii) = -2 * (fnorm2(ii) + maxL(ii)) + M * log(sum(N2));

%
%  Akaike Information Criterion (for comparison)
% -----------------------------------------------------------------------------
AIC(ii) = -2 * (fnorm2(ii) + maxL(ii)) + 2 * M;

fprintf(1,'%s,\t %4i, %4i, %12.4f, %4i \n',ppd2, M,sum(N2),log(sum(N2)),Nf);

clear m A gamma phi F D ppd2 covd2;

return;
