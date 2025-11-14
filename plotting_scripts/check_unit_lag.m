clear;close all;

x=load('fort.77');
mw=load('fort.78');
C1=load('fort.81');

NFP = size(mw,2);
Cx=cov(x);
Cmd=cov(mw);
ncov = size(mw,1);
ncov2 = 100;

mcsum1 = zeros(NFP,NFP);
for ismp=1:ncov2;
for ipar=1:NFP;
  for jpar=1:NFP;
    mcsum1(ipar,jpar) = mcsum1(ipar,jpar)+mw(ismp,ipar)*mw(ismp,jpar);
  end;
end;
end;
Cov1=mcsum1/ncov2;

figure();hold on;
plot(diag(Cx));
figure();hold on;
plot(diag(Cmd),'-b');
plot(diag(Cov1),'--r');
plot(diag(C1),'--k');
