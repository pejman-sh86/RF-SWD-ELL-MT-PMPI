function [N] = corr_noise(sd,nfreq,nang,icorr)
%
% loop over freqs 315-1600 Hz
%
for ifreq = 1:nfreq

  %
  % G(0,1 noise):
  %
  noise = randn(nang(ifreq),1);
  
  %
  % Set up laplacian cov mat:
  %
  lag = [-nang(ifreq):nang(ifreq)];
  b = 4.;
  laplace = 1./(2.*b)*exp(-abs(lag)/b);
  laplace = laplace/max(laplace);
  laplace = sd^2*laplace;

  C = zeros(nang(ifreq),nang(ifreq));
  for i=1:nang(ifreq)
    nang(ifreq) - i + 2 
    length(laplace)-i+1
    C(i,:) = laplace(nang(ifreq) - i + 2: end-i);
  end
  imagesc(C);

  % 
  % Decompose: C = L'*L
  %
  L = chol(C);

  %
  % Draw samples with cov C:
  %
  cnoise = L*noise;
 
  N(ifreq).covnn = C;
  N(ifreq).cn = cnoise;
  N(ifreq).laplace = laplace;

end

return;

