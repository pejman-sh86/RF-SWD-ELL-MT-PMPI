function [N] = corr_noise(sd,nfreq,nang,icorr)
%
% loop over freqs 315-1600 Hz
%
for ifreq = 1:nfreq

  randn('state',137)
  noise = zeros(nang(ifreq));
  noise = randn(nang(ifreq),1);
  
%
% Set up art. cov. mat:
%

  Cdiag = zeros(nang(ifreq),9);
  Cdiag(:,1) = (sd(ifreq)^2)/5*ones(nang(ifreq),1);
  Cdiag(:,2) = (sd(ifreq)^2)/4*ones(nang(ifreq),1);
  Cdiag(:,3) = (sd(ifreq)^2)/3*ones(nang(ifreq),1);
  Cdiag(:,4) = (sd(ifreq)^2)/2*ones(nang(ifreq),1);
  Cdiag(:,5) = sd(ifreq)^2*ones(nang(ifreq),1);
  Cdiag(:,6) = (sd(ifreq)^2)/2*ones(nang(ifreq),1);
  Cdiag(:,7) = (sd(ifreq)^2)/3*ones(nang(ifreq),1);
  Cdiag(:,8) = (sd(ifreq)^2)/4*ones(nang(ifreq),1);
  Cdiag(:,9) = (sd(ifreq)^2)/5*ones(nang(ifreq),1);

  C = zeros(nang(ifreq),nang(ifreq));
  C = spdiags(Cdiag,-4:4,C);
  clear Cdiag;
  % 
  % Decompose: C = L'*L
  %
  L = chol(C);
  %
  % Draw samples with cov C:
  %
  cnoise = L'*noise;
  if(icorr == 1)
    %
    % Zero padding:
    %
%    npad = floor(nang(ifreq)/6);
    npad = floor(nang(ifreq)-4);
    %
    % Calc autocovariance:
    %
    Antnt = xcov(cnoise)/nang(ifreq);
    %
    % Set up covariance matrix in a diagonal form:
    %
    Cdiag = zeros(nang(ifreq),(2*nang(ifreq))-1);
    for i = 1:(2*nang(ifreq))-1
      Cdiag(:,i) = Antnt(i)*ones(nang(ifreq),1);
    end
    Cdiag(:,1:npad) = 0;
    Cdiag(:,(2*nang(ifreq))-npad:(2*nang(ifreq))-1) = 0;
    E = diag(ones(nang(ifreq),1));
    Chat = zeros(nang(ifreq),nang(ifreq));
    Chat = spdiags(Cdiag,-nang(ifreq)+1:nang(ifreq)-1,Chat);
    Csave = E*Chat;
 
    N(ifreq).cov = Csave;
    N(ifreq).covnn = C;

  else
    E = diag(ones(nang(ifreq),1));
    N(ifreq).cov = sd(ifreq)^2 * diag(ones(nang(ifreq),1));
%    N(ifreq).cn = sd(ifreq) * noise;
  end
  N(ifreq).cn = cnoise;

end

return;

