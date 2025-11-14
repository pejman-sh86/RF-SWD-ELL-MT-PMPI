function est_cov_mat
%
% Test real data:
%
nfreq = 8;
nang = 135;
A = load('residuals.dat');

%
% loop over freqs 315-1600 Hz
%
ttl = [{'315 Hz'},{'400 Hz'},{'500 Hz'},{'630 Hz'},{'800 Hz'},...
       {'1000 Hz'},{'1250 Hz'},{'1600 Hz'},{'2000 Hz'}];
k = 1;
for ifreq = 1:nfreq
  res = zeros(nang);
  res = A(ifreq,:);
%
% Get rid of NaNs:
%
  j = 1;
  for i=1:nang
    if(isnan(res(i)) ~= 1)
      tmp(j) = res(i);
      j = j+1;
    end
  end
  newnang = j-1;
  res = tmp;
%
% Zero padding:
%
  npad = floor(newnang/6)
%
% Calc autocovariance:
%
  Aresres = xcov(res)/newnang;
%
% Set up covariance matrix in a diagonal form:
%
  sd3 = zeros(newnang,(2*newnang)-1);
  for i = 1:(2*newnang)-1
    sd3(:,i) = Aresres(i)*ones(newnang,1);
  end
  sd3(:,1:npad) = 0;
  sd3(:,(2*newnang)-npad:(2*newnang)-1) = 0;
  E = diag(ones(newnang,1));
  Chat = zeros(newnang,newnang);
  Chat = spdiags(sd3,-newnang+1:newnang-1,Chat);
  Csave = E*Chat;
%  save cov_mat_est.mat Csave;

  F(ifreq).csave = Csave;
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(Chat);
%
% "Uncorrelate" the residuals:
%
  reshat = inv(L3')*res';
  Aresres2 = xcov(reshat)/newnang;
  
  figure(1);
  hold on;
  subplot(4,4,k)
  plot(Aresres,'-k.')
  text(100,1,ttl(ifreq))
  if k>12; xlabel('lag');end;
  if ((k==1) | (k==5) | (k==9) | (k==13)); ylabel('A_{xx}');end;
%  if k==5; ylabel('A_{xx}');end;
%  if k==9; ylabel('A_{xx}');end;
%  if k==13; ylabel('A_{xx}');end;
  set(gca,'YLim',[-0.5 1.5]);
  subplot(4,4,k+1)
  plot(Aresres2,'-k.')
  text(100,1,ttl(ifreq))
  if k>12; xlabel('lag');end;
  set(gca,'YLim',[-0.5 1.5]);

  x = -5:1:5;
  figure(2);
  hold on;
  subplot(4,4,k);
  [n1,xout] = hist(res,x);
  n1 = n1/(newnang*(xout(2)-xout(1)));
  bar(xout,n1);
  title(ttl(ifreq));
  if k>12; xlabel('lag');end;
  if k==1; ylabel('A_{xx}');end;
  if k==5; ylabel('A_{xx}');end;
  if k==9; ylabel('A_{xx}');end;
  if k==13; ylabel('A_{xx}');end;
  axis([-5 5 0 0.8]);
  subplot(4,4,k+1);
  [n1,xout] = hist(reshat,x);
  n1 = n1/(newnang*(xout(2)-xout(1)));
  bar(xout,n1);
  title(ttl(ifreq));
  if k>12; xlabel('lag');end;
  axis([-5 5 0 0.8]);
  k = k+2;

  z = runtest(res');
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(1,'%10.4f',abs(z));
  if abs(z)<1.96
    fprintf(1,'\t passed');
  else
    fprintf(1,'\t failed');
  end
  fprintf(1,'\n');
  z = runtest(reshat);
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(1,'%10.4f',abs(z));
  if abs(z)<1.96
    fprintf(1,'\t passed');
  else
    fprintf(1,'\t failed');
  end
  fprintf(1,'\n');
end
save cov_mat_est.mat F;

return;

