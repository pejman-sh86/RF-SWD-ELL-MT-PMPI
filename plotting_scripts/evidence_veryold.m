%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [E,E2,logE,log10E,L,minL,occ] = evidence();

ifirst = 0;
ppd   = 'sim_C_1_7_40_6lay_sample.mat';
covd  = 'sim_C_1_7_40_6lay_cov.mat';
hpdf  = 'sim_C_1_7_40_6lay_hpds.txt';

%ppd   = 'x_s16_1_5_40_6layrgcov_sample.mat';
%covd  = 'x_s16_1_5_40_6layrgcov_cov.mat';

%ppd   = 'x_s13_2_10_50_5lay2_sample.mat';
%covd  = 'x_s13_2_10_50_5lay_cov.mat';

%ppd   = 'x_s02_1_3_16_5laycov_sample.mat';
%covd  = 'x_s02_1_3_16_5laycov_cov.mat';
%hpdf  = 'x_s02_1_3_16_5laycov_hpds.txt';

if(ifirst == 0)
   load('mindetCd.mat')
end

load(ppd);
load(covd);
hpd = load(hpdf);

m = A(:,2:end)';
phi = A(:,1)';

Q = length(phi);		% Size of sample
M = size(m,1);			% Number of model parameters
N = size(F(1).csave,1);		% Number of data per frequency

%
% Values from sample:
%
L = mean(A(:,1));
minL = min(A(:,1));

mu = mean(m')';
Cm = cov(m');
invCm = inv(Cm);
detCm = det(Cm);

%
% Ockham Factor
%
occ = log10(minL*prod((hpd(:,2)-hpd(:,1))./(F(1).maxlim-F(1).minlim)'));

%
% Volume of the model space:
%
V = prod(F(1).maxlim-F(1).minlim);

%
% Normalization for f(m_i)
%
for i = 1:length(F(1).freq)
  
  
%  F(i).csave = 0.05*eye(N,N);
  F(i).csave = N*F(i).csave;
  tmp(i) = det(F(i).csave);

end

if(ifirst == 1)
   mindetCd = max(tmp);
   save mindetCd.mat mindetCd;
end

fnorm = 0;
for i = 1:length(F(1).freq)

  detCd = det(F(i).csave)/mindetCd;
  sqrt(detCd)
  fnorm = fnorm + log10(1/sqrt(detCd));

end

fnorm = fnorm/V;

gnorm = 1/((2*pi)^(M/2)*sqrt(detCm));
for i = 1:Q;

   g2(i) = ((0.5*(m(:,i)-mu)'*invCm*(m(:,i)-mu)));
   psi(i) = -(g2(i)-phi(i));

end;

mean(psi);
%B = min(psi);
%log10E = log10(Q) + fnorm - B - log10(gnorm) - log10(sum(exp(psi-B)));
log10E = log10(Q) + fnorm - log10(gnorm) - log10(sum(exp(psi)));

logE = log10E/log10(exp(1));

E = 10^(logE);

E2 = Q*fnorm/gnorm*(1/sum(exp(psi)));

return;
