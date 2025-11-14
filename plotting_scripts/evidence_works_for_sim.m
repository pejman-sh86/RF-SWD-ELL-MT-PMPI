%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [log10E,L,minL,occ] = evidence(filebase);

ifirst = 0;
ppd   = strcat(filebase,'_sample.mat');
covd  = strcat(filebase,'_cov.mat');

%ppd   = 'x_s16_1_5_40_6layrgcov_sample.mat';
%covd  = 'x_s16_1_5_40_6layrgcov_cov.mat';

%ppd   = 'x_s13_2_10_50_2lay_sample.mat';
%covd  = 'x_s13_2_10_50_2lay_cov.mat';

%ppd   = 'x_s02_1_3_16_5laycov_sample.mat';
%covd  = 'x_s02_1_3_16_5laycov_cov.mat';

ii=1;
%for ii=1:7;

%   ppd2  = strcat(ppd(1:13),num2str(ii),ppd(15:end));
%   covd2 = strcat(covd(1:13),num2str(ii),covd(15:end));

   if(ifirst == 0)
      load('mindetCd.mat')
   end

   load(ppd);
   load(covd);

   m = A(:,2:end)';
   phi = A(:,1)';

   Q = length(phi);		% Size of sample
   M = size(m,1);			% Number of model parameters
   N = size(F(1).csave,1);		% Number of data per frequency

%
% Values from sample:
%
   L(ii) = mean(A(:,1));
   minL(ii) = min(A(:,1));

   mu = mean(m')';
   Cm = cov(m');
   invCm = inv(Cm);
   detCm = det(Cm);

%
% Ockham Factor
%
   occ(ii) = log10(minL(ii)*prod((max(m')-min(m'))./(F(1).maxlim-F(1).minlim)));

%
% Volume of the model space:
%
   V = prod(F(1).maxlim-F(1).minlim);

%
% Normalization for f(m_i)
%
   for i = 1:length(F(1).freq)
  
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
     fnorm = fnorm + log10(1/((2.*pi)^(N/2)*sqrt(detCd)));
     fnorm = fnorm + log10(1/sqrt(detCd));
   end

   fnorm = fnorm - log10(V);

   gnorm = log10(1/((2*pi)^(M/2)*sqrt(detCm)));
   for i = 1:Q;

      g2(i) = ((0.5*(m(:,i)-mu)'*invCm*(m(:,i)-mu)));
      psi(i) = -(g2(i)-phi(i));

   end;

   B = min(psi);
   log10E(ii) = log10(Q) + fnorm - gnorm - log10(sum(exp(psi-B)))+(log10(exp(1))*B);

   clear hpd m A psi B fnorm gnorm g2 psi;
%end;

return;
