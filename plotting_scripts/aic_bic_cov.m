%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [minL,BIC,AIC,C] = aic_bic(filebase,len,len2);

ifirst = 0;
isd = 0;
resfile = strcat(filebase,'_res.txt');
data    = strcat(filebase(1:len2),'.mat');
covfile = strcat(filebase(1:len),'_cov.mat');
%covfile = 'x_s21_1_5_25_5lay_cov.mat';

res = load(resfile);
load(covfile);
D = F; clear F;
ii=1;

load(data);
C(ii).file = filebase;

M  = length(F(1).minlim);	% Number of model parameters
N  = size(F(1).dat,1);		% Number of data per frequency
Nf = length(F(1).freq);		% Number of data per frequency

%
% Compute misfit:
% --------------------------------------------------------------------------------
res = F(1).Rex(:,1:Nf)'.*res(1:Nf,:);

%
% Normalization for f(m_i)
% --------------------------------------------------------------------------------
fnorm(ii)  = 0;
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

end

fnorm(ii) = -fnorm(ii);

minL(ii) = 0.;
for i = 1:length(F(1).freq)
  
  idx = find(F(1).Rex(1:N,i) == 1);
  N2(i) = length(idx);

%  sd2(i) = sqrt(1/(N2(i)-(M/Nf))*sum(res(i,:).^2));

  minL(ii) = minL(ii) + 0.5*sum(res(i,:)*inv(D(i).csave)*res(i,:)');
  minL(ii)


end
minL(ii)

%
% Bayesian Information Criterion (for comparison)
% --------------------------------------------------------------------------------
BIC(ii) = -2 * (fnorm(ii) - minL(ii)) + M * log(sum(N2));

%
%  Akaike Information Criterion (for comparison)
% --------------------------------------------------------------------------------
AIC(ii) = -2 * (fnorm(ii) - minL(ii)) + 2 * M;


fprintf(1,'%s\n',resfile);
fprintf(1,'%s,\t %4i, %12.4f, %12.4f, %4i, %12.4f \n',filebase, M,minL(ii),M*log(sum(N2)),Nf,fnorm(ii));

return;
