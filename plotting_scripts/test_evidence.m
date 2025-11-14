%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [log10E,L,minL,occ,BIC,AIC,V,Q,M,N,Cm,detCm,detCd2,fnorm,fnorm2,fnorm3,C] = evidence(filebase,len);

ifirst = 0;
isd = 1;
ppd   = strcat(filebase,'_sample.mat');
%covd  = strcat(filebase,'_cov.mat');
%covd  = 'sim_C_1_7_40_4layrg_cov.mat';
%covd  = 'sim_B_1_10_50_4lay_cov.mat';
%covd   = 'x_s13_2_10_50_5layb_cov.mat';
%sdfile = 'x_s13_2_10_50_exstd.txt';
%covd  = 'x_s02_1_3_16_4lay_cov.mat';
covd   = 'x_s01_1_3_20_5layb_cov.mat';
sdfile = 'x_s01_1_3_20_exstd.txt';
data  = strcat(filebase(1:len),'.mat');

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
end

   load(data2);
   C(ii).ppd = ppd2;
   disp(ppd2);

   m = A(:,2:end)';
   phi = A(:,1)';

   Q  = length(phi);		% Size of sample
   M  = size(m,1);		% Number of model parameters
   N  = size(F(1).ang,2)	% Number of data per frequency
   Nf = length(F(1).freq);	% Number of data per frequency

%
% Values from sample:
% --------------------------------------------------------------------------------
   L(ii) = mean(A(:,1));
   minL(ii) = min(A(:,1));

   mu = mean(m')';
   Cm = cov(m');
   invCm = inv(Cm);
   detCm = det(Cm);

%
% Volume of the model space:
% --------------------------------------------------------------------------------
   V(ii) = prod(F(1).maxlim-F(1).minlim);
   log10V(ii) = log10(V(ii));

%
% Scale data covariance matrix for numerical stability
% --------------------------------------------------------------------------------
if(isd ==0)
   for i = 1:length(F(1).freq)
     D(i).csave = N*D(i).csave;
     tmp(i) = det(D(i).csave);
   end

   if(ifirst == 1)
      mindetCd = max(tmp);
      save mindetCd.mat mindetCd;
   end
   if(ifirst == 0)
      load('mindetCd.mat')
   end
end
%
% Normalization for f(m_i)
% --------------------------------------------------------------------------------
   fnorm(ii)  = 0;
   fnorm2(ii) = 0;
   for i = 1:length(F(1).freq)

  if(isd == 0)
     detCd(i)  = det(D(i).csave)/mindetCd;
     detCd2(i) = det(D(i).csave);
  else
     detCd(i)  = det(sd(i) * eye(N));
     detCd2(i)  = det(sd(i) * eye(N));
  end


     idx = find(F(1).Rex(1:N,i) == 1);
     N2(i) = length(idx);

     fnorm(ii)  = fnorm(ii) + log10((2*pi)^(N/2)) + log10(sqrt(detCd(i)));
     fnorm2(ii) = fnorm2(ii) + log((2*pi)^(N/2)) + log(sqrt(detCd(i)));

   end

   fnorm3(ii) = - fnorm(ii);
   fnorm(ii) = -log10V(ii) - fnorm(ii);
   fnorm2(ii) = - fnorm2(ii);

%
% Normalization for g(m_i)
% --------------------------------------------------------------------------------
   gnorm(ii) = -log10( (2*pi)^(M/2) ) - log10( sqrt(detCm) );

%
% Compute gamma (exponent of g(m_i))
% --------------------------------------------------------------------------------
%   for i = 1:Q;
%
%      gamma(i) = 0.5*(m(:,i)-mu)'*invCm*(m(:,i)-mu);
%
%   end;

%
% Compute log10 of evidence 
% --------------------------------------------------------------------------------
%   B = min(phi);
%   log10E(ii) = log10(Q) - gnorm(ii) + fnorm(ii) - (log10(exp(1))*B) - ...
%                log10(sum( exp( - gamma + phi - B ) ));
   log10E(ii) = 0;

%
% Ockham Factor (for comparison)
% --------------------------------------------------------------------------------
   occ(ii) = (fnorm3(ii) - minL(ii)/log(10)) + log10(1/V(ii)*sqrt(detCm)*(2*pi)^(M/2));

%
% Bayesian Information Criterion (for comparison)
% --------------------------------------------------------------------------------
   BIC(ii) = -2 * (fnorm2(ii) - minL(ii)) + M * log(sum(N2));

%
%  Akaike Information Criterion (for comparison)
% --------------------------------------------------------------------------------
   AIC(ii) = -2 * (fnorm2(ii) - minL(ii)) + 2 * M;

fprintf(1,'%s,\t %4i, %4i, %12.4f, %4i \n',ppd2, M,sum(N2),log(sum(N2)),Nf);

   clear hpd m A gamma phi B F D ppd2 covd2;

return;
