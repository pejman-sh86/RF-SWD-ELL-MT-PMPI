function est_cov_mat_har
%
% Test real data:
%
% UNCOMMENT SAVE STATEMENT!!!!
%
dex = 0;
npar = 7;

tab = readtextfile('measured_parameters_3b.dat');
[outfile,F,converg_cov,converg_marg] =...
convert_parfile(tab);

sample = 'sample.dat';
residual_file = strrep(sample,'sample.dat','residuals_3b_0.mat');
cov_mat_file = strrep(sample,'sample.dat','cov_mat_est_3b_0.mat');

load(residual_file);
nfreq = length(F(1).freq);
ang = F(1).ang(F(1).astart:F(1).aend);
nang = length(ang);


k = 1;
for ifreq = 1:nfreq
  res = (B(ifreq).diff(F(1).astart:F(1).aend))';
%
% Zero padding:
%
%  npad = floor(newnang/6)
  npad = 1;
%
% Calc autocovariance:
%
  [Aresres,lag1] = xcov(res,'biased');
%  [Aresres,lag1] = xcov(res,'coeff');
%
% Set up covariance matrix in a diagonal form:
%
  sd3 = zeros(nang,(2*nang)-1);
  for i = 1:(2*nang)-1
    sd3(:,i) = Aresres(i)*ones(nang,1);
  end
  sd3(:,1:npad) = 0;
  sd3(:,(2*nang)-npad:(2*nang)-1) = 0;
  E = diag(ones(nang,1));
  Chat = zeros(nang,nang);
  Chat = spdiags(sd3,-nang+1:nang-1,Chat);
  Csave = E*Chat;

  C(ifreq).csave = Csave;
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(Chat);
%
% "Uncorrelate" the residuals:
%
  reshat = inv(L3')*res;
  [Aresres2,lag2] = xcov(reshat,'biased');
%  [Aresres2,lag2] = xcov(reshat,'coeff');
  
  figure(1);
  hold on;
  subplot(6,4,k)
  hold on, box on;
  plot(lag1,Aresres,'-k.')
  plot([-40 40],[0 0],'k')
%  text(100,1,ttl(ifreq))
  title([num2str(F(1).freq(ifreq)) ' Hz'])

%  if k>12; xlabel('lag');end;
%  if ((k==1) | (k==5) | (k==9) | (k==13)); ylabel('A_{xx}');end;
%  if k==5; ylabel('A_{xx}');end;
%  if k==9; ylabel('A_{xx}');end;
%  if k==13; ylabel('A_{xx}');end;
%  set(gca,'YLim',[-0.5 1.5]);
  subplot(6,4,k+1)
  hold on, box on;
  plot(lag2,Aresres2,'-k.')
  plot([-40 40],[0 0],'k')
%  text(100,1,ttl(ifreq))
  title([num2str(F(1).freq(ifreq)) ' Hz'])
%  if k>12; xlabel('lag');end;
%  set(gca,'YLim',[-0.5 1.5]);

  x = -.1:.02:.1;
  figure(2);
  hold on;
  subplot(6,4,k);
  hold on, box on;
  [n1,xout] = hist(res,x);
  n1 = n1/(nang*(xout(2)-xout(1)));
  bar(xout,n1);
%  title(ttl(ifreq));
  title([num2str(F(1).freq(ifreq)) ' Hz'])
%  if k>12; xlabel('lag');end;
%  if k==1; ylabel('A_{xx}');end;
%  if k==5; ylabel('A_{xx}');end;
%  if k==9; ylabel('A_{xx}');end;
%  if k==13; ylabel('A_{xx}');end;
%  axis([-5 5 0 0.8]);
  subplot(6,4,k+1);
  hold on, box on;
  x = -3:0.5:3;
  [n1,xout] = hist(reshat,x);
  n1 = n1/(nang*(xout(2)-xout(1)));
  bar(xout,n1);
%  title(ttl(ifreq));
  title([num2str(F(1).freq(ifreq)) ' Hz'])
%  if k>12; xlabel('lag');end;
%  axis([-5 5 0 0.8]);


  figure(3);
  hold on;
  subplot(6,4,k);
  hold on, box on;
  plot(ang,res);
  plot([0 90],[0 0],'k')
  set(gca,'XLim',[0 90],'YLim',[-.1 .1 ]);
%  title(ttl(ifreq));
  title([num2str(F(1).freq(ifreq)) ' Hz'])
%  if k>12; xlabel('lag');end;
%  if k==1; ylabel('A_{xx}');end;
%  if k==5; ylabel('A_{xx}');end;
%  if k==9; ylabel('A_{xx}');end;
%  if k==13; ylabel('A_{xx}');end;
%  axis([-5 5 0 0.8]);
  subplot(6,4,k+1);
  hold on, box on;
  plot(ang,reshat);
  plot([0 90],[0 0],'k')
  set(gca,'XLim',[0 90],'YLim',[-3 3]);
%  title(ttl(ifreq));
  title([num2str(F(1).freq(ifreq)) ' Hz'])
%  if k>12; xlabel('lag');end;
%  axis([-5 5 0 0.8]);





  k = k+2;

  z = runtest(res);
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

save(cov_mat_file,'C');

return;

