function gen_corr_err

r = randn(100,1);
Arr = xcov(r)/100;

%
% Set up art. cov. mat:
%
sd(:,1) = 1*ones(100,1);
sd(:,2) = 2*ones(100,1);
sd(:,3) = 3*ones(100,1);
sd(:,4) = 4*ones(100,1);
sd(:,5) = 3*ones(100,1);
sd(:,6) = 2*ones(100,1);
sd(:,7) = 1*ones(100,1);

E = diag(ones(100,1));
C = zeros(100,100);
C = spdiags(sd,-3:3,C);

%
% Decompose: C = L'*L
%
L = chol(C);
%
% Draw samples with cov C:
%
rtilde = L'*r;
Arr2 = xcov(rtilde)/100;
%
% Check:
%
rhat = inv(L')*rtilde;
Arr3 = xcov(rhat)/100;

%
% A) Put full xcov in matrix:
%
for i = 1:199
  sd2(:,i) = Arr2(i)*ones(100,1);
end
Ctilde = zeros(100,100);
Ctilde = spdiags(sd2,-99:99,Ctilde);

save bla C Ctilde;

L2 = chol(Ctilde);

rhat2 = inv(L2')*rtilde;
Arr4 = xcov(rhat2)/100;

figure(1);
subplot(4,1,1)
plot(Arr,'-+k')
subplot(4,1,2)
plot(Arr2,'-+k')
subplot(4,1,3)
plot(Arr3,'-+k')
subplot(4,1,4)
plot(Arr4,'-+k')

return;
