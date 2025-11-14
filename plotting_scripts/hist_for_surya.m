function []=res_analysis_surya();

xx = -6.25:.01:6.25;
nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);

rep1 = dlmread('noisy_data_rep_with_AR.dat');
rep2 = dlmread('noisy_data_rep_without_AR.dat');
obs  = dlmread('noisy_data_obs.dat');
obs=obs(39:91,:);
rep1=rep1(39:91,:);
rep2=rep2(39:91,:);

%% Remove mean from residuals:
res1=obs(:,2)-rep1(:,2);
res2=obs(:,2)-rep2(:,2);
res1=res1-mean(res1);
res2=res2-mean(res2);

%% Standardize with standard deviation:
sd1=std(res1);
sd2=std(res2);
res1st = res1/sd1;
res2st = res2/sd2;

%% bins for histograms:
dn = 0.3;
lim = [-3:dn:3];

[n1,lim]=hist(res1st(:),lim);
n1=n1/sum(n1)/dn;
[n2,lim]=hist(res2st(:),lim);
n2=n2/sum(n2)/dn;

%% Need to shift bins 1/2 to left for stairs command:
lim2 =lim - dn/2.;

figure();hold on; box on;

%% Axx
subplot(2,2,1);hold on;box on;
plot([-52:52],xcorr(res2st,'coeff'),'o-k')
plot([-52 52],[0 0],'--k','LineWidth',1);
xlim([-52 52]);
ylim([-0.45 1.1]);
subplot(2,2,3);hold on;box on;
plot([-52:52],xcorr(res1st,'coeff'),'o-k')
plot([-52 52],[0 0],'--k','LineWidth',1);
xlim([-52 52]);
ylim([-0.45 1.1]);
xlabel('Lag (s)');

%% Res hists
subplot(2,2,2);hold on;box on;
stairs(lim2,n2,'-k');
plot(xx,nd,'--k')
xlim([-4 4]);
ylim([ 0 0.65]);
subplot(2,2,4);hold on;box on;
stairs(lim2,n1,'-k');
plot(xx,nd,'--k')
xlim([-4 4]);
ylim([ 0 0.65]);
xlabel('Standard deviations');

return;
