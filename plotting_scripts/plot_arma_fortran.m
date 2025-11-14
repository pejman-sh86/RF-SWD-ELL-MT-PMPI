function []=plot_arma_fortran();

load fort.73
res1 = fort;
load fort.74
res2 = fort;
load fort.75
obs  = fort;
load fort.76
rep  = fort;
load fort.77
ar   = fort;
load fort.78
res3 = fort;
load fort.79
ari = fort;
N= length(res1);
f1 = res1(1,1);
f2 = res1(1,end);

figure();hold on; box on;
plot(obs(1,:),obs(2,:),'-k');
plot(rep(1,:),rep(2,:),'-r');
legend('obs','rep');
plot([f1 f2],[0 0],':k');

figure();hold on; box on;
plot(res1(1,:),res1(2,:),'-k');
plot(ar(1,:),ar(2,:),'-r');
plot(ari(1,:),ari(2,:),'-b');
legend('res','ar','ari');
plot([f1 f2],[0 0],':k');

figure();hold on; box on;
plot(res1(1,:),res1(2,:),'-k');
plot(res2(1,:),res2(2,:),'-r');
plot(res3(1,:),res3(2,:),'-b');
legend('res1','res2','res3');
plot([f1 f2],[0 0],':k');

%%
%% Autocorrelation
%%
figure();hold on;box on;
[acs1,lags1] = xcorr(res1(2,:),'coeff');   % ACS of prediction error
[acs2,lags2] = xcorr(res2(2,:),'coeff');   % ACS of prediction error
[acs3,lags3] = xcorr(res3(2,:),'coeff');   % ACS of prediction error
plot(lags1,acs1,'-k')
plot(lags2,acs2,'-r')
plot(lags3,acs3,'-b')
plot([lags1(1) lags1(end)],[0 0],':k');
set(gca,'XLim',[-N N]);
xlabel('Lag');ylabel('ACF');
legend('Raw Signal','I applied','ARI applied')

return;
