function [E,occ,L,fnorm,sigma_wd] = occam();

ifirst = 0;

covd   = '../sim_C_1_7_40_4layrg_cov.mat';
ppd    =    'sim_C_1_7_40_4layrg_sample.mat';
hpd=   load('sim_C_1_7_40_4layrg_hpds.txt');

if(ifirst == 0)
   load('mindetCd.mat')
end

load(ppd);
load(covd);

m = A(:,2:end)';
phi = A(:,1)';
L = min(phi);

Q = length(phi);		% Size of sample
M = size(m,1);			% Number of model parameters
N = size(F(1).csave,1);		% Number of data per frequency

sigma_w  = F(1).maxlim'-F(1).minlim';
sigma_wd = hpd(:,2)-hpd(:,1);

occ = prod(sigma_wd./sigma_w);

%
% Normalization for f(m_i)
%
for i = 1:length(F(1).freq)
  
  F(i).csave = 0.05*eye(N,N);
  F(i).csave = N*F(i).csave;
  tmp(i) = det(F(i).csave);

end

if(ifirst == 1)
   mindetCd = min(tmp);
   save mindetCd.mat mindetCd;
end

fnorm = 1;
for i = 1:length(F(1).freq)

  detCd = det(F(i).csave)/mindetCd;
  fnorm = fnorm * 1/sqrt(detCd);

end

f = fnorm *exp(-L);

E = f* occ;


return;
