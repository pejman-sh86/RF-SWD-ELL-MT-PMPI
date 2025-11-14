function chck_est
%
% Test real data:
% Take ML model, calc residuals, use C_d estimate to 
% decorrelate, plot histograms of that and Autocovariance.
% Then do runs and KS test with residuals.
%
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz
%load map.mat;

resfile1     = '../run_assa/x_jan_2_sph_res.txt';
resfile2     = 'x_jan_2_sph_ci2w_res.txt';
cov_mat_file = '../run_assa/x_jan_2_sph_ci_cov.mat';
plotfile1    = 'x_jan_2_sph_ci2w_autocorr_compare.eps';
plotfile2    = 'x_jan_2_sph_ci2w_res_hist.eps';
plotfile3    = 'x_jan_2_sph_ci2w_ks_test.eps';
logfile      = 'x_jan_2_sph_ci2w_chck_log2.txt';

fid = fopen(logfile,'w');

load('critical_values.mat');
load(cov_mat_file);
%load('bla');
res2_tmp = load(resfile2);
res1_tmp = load(resfile1);
%load(stdfile);
%ml_sd = stdv_dB(:,2);
aincr    =  1;          %
a_start  =  1;          %
nang_tmp =135;          % Nr of angles
fincr    =  1;          %
f_start  =  6;          %

freq = F(1).freq;
nfreq    =  length(F(1).freq);
for i = 1:nfreq
  nang(i) = length(F(i).csave);
  ml_sd(i) = sqrt(F(i).csave(1,1));
end

llw = [0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ...
       0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ];
llh1 = [0.77 0.58 0.34 0.15];
llh2 = [0.765 0.56 0.355 0.15];

i = 1; j = 1; k = 1; l = 1;

for ifreq = 1:nfreq
  res1 = res1_tmp(ifreq,:)';
  res2 = res2_tmp(ifreq,:)';
  Aresres = xcov(res1)/nang(ifreq);
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(F(ifreq).csave);
%
% "Uncorrelate" the residuals:
%
  reshat = inv(L3')*res2;
  Aresres2 = xcov(reshat)/nang(ifreq);
 
%
% Plot Autocorrelations
%
  gc1 = figure(1);
  hold on;
   if(k <= 8)
     subplot('Position',[llw(j) llh1(1) 0.18 0.18]);
   else
     subplot('Position',[llw(j) llh1(3) 0.18 0.18]);
   end
  lag = -floor(length(Aresres)/2):1:floor(length(Aresres)/2);
  plot(lag,Aresres/max(Aresres),'-k.');
  text(length(Aresres2)/6,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8);
  set(gca,'FontSize',8);
  set(gca,'XLim',[lag(1) lag(end)],'YLim',[-0.5 1.1]);
  if((k==1)|(k==9)); 
%    ylabel('Axx');
     ylabel('A')
  end;
  set(gca,'XTick',[-100 -50 0 50 100],'XTickLabel',[]);
  set(gca,'YTick',[-0.5 0 0.5 1],'YTickLabel',{'';'';'';''});
  if(k==1); set(gca,'YTickLabel',{'';'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'';'0';'';'1'});end;


  if k+1 <= 8
    subplot('Position',[llw(j) llh1(2) 0.18 0.18]);
    j = j + 1;
  else
    subplot('Position',[llw(j) llh1(4) 0.18 0.18]);
    j = j + 1;
  end
  lag2 = -floor(length(Aresres2)/2):1:floor(length(Aresres2)/2);
  plot(lag2,Aresres2/max(Aresres2),'-k.');
  text(length(Aresres2)/6,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8);
%  set(gca,'FontSize',8);
  set(gca,'FontSize',10);
  set(gca,'XLim',[lag2(1) lag2(end)],'YLim',[-0.5 1.1]);
  if(k>=9);xlabel('lag');end;
  if ((k==1) | (k==9)); ylabel('Axx');end;
  set(gca,'XTick',[-100 -50 0 50 100],'XTickLabel',[-100 -50 0 50 100]);
  set(gca,'YTick',[-0.5 0 0.5 1],'YTickLabel',{'';'';'';''});
  if(k==1); set(gca,'YTickLabel',{'';'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'';'0';'';'1'});end;


%
% Plot histograms
%
  gc2 = figure(2);
  hold on;

   if(k <= 8)
     subplot('Position',[llw(i) llh1(1) 0.18 0.18]);
   else
     subplot('Position',[llw(i) llh1(3) 0.18 0.18]);
   end
  hold on;box on;
  x = -3:1/(nang(ifreq))*20:3;
  xx = -3:.01:3;
  [n1,xout] = hist(res1/ml_sd(ifreq),x);
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  area = sum(n1) * (x(2)-x(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(0.8,1.5,[num2str(freq(ifreq)) ' Hz'],'FontSize',8)
  set(gca,'FontSize',8);
  set(gca,'XTick',[-2 0 2]);
  set(gca,'XTickLabel',[]);
  set(gca,'YTick',[0 0.5 1 1.5]);
  set(gca,'YTickLabel',{'';'';''});
  if(k==1); set(gca,'YTickLabel',{'0';'';'1';''});end;
  if(k==9); set(gca,'YTickLabel',{'0';'';'1';''});end;
  axis([-3 3 0 1.7]);

  if k+1 <= 8
    subplot('Position',[llw(i) llh1(2) 0.18 0.18]);
    i = i + 1;
  else
    subplot('Position',[llw(i) llh1(4) 0.18 0.18]);
    i = i + 1;
  end
  hold on;box on;
  x = -3:1/(nang(ifreq))*20:3;
  xx = -3:.01:3;
  [n1,xout] = hist(reshat,x);
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  area = sum(n1) * (x(2)-x(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(0.8,1.5,[num2str(freq(ifreq)) ' Hz'],'FontSize',8)
  set(gca,'FontSize',8);
  if k>=9; xlabel('Residual [\sigma]');end;
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2';'0';'2'});
  set(gca,'YTick',[0 0.5 1 1.5],'YTickLabel',{'';'';'';''});
  if(k==1); set(gca,'YTickLabel',{'0';'';'1';''});end;
  if(k==9); set(gca,'YTickLabel',{'0';'';'1';''});end;
  axis([-3 3 0 1.7]);


%
% Test uncorrected
%
  z = runtest(res1);
  p = normcdf(z,0,1);
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(1,'%10.4f\t%10.4f',abs(z),p);
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(fid,'%10.4f\t%10.4f',abs(z),p);
  if abs(z)<1.96
    fprintf(1,'\t run passed');
    fprintf(fid,'\t run passed');
  else
    fprintf(1,'\t run failed');
    fprintf(fid,'\t run failed');
  end
  [pass,value,bereik,opp,dplot] = kolmogorov(res1/ml_sd(ifreq),k);
  if(value < foo(end,3))
      idx2 = find(value<foo(:,3));
      idx2 = idx2(1);
  else
      idx2 = length(foo);
  end
  p1 = (100-foo(idx2,2))/100;
  fprintf(fid,'   ks:');
  fprintf(fid,' %7.4f ',value,p1);
  if p1 >= 0.05
    fprintf(1,'  ks passed');
    fprintf(fid,'  ks passed');
  else
    fprintf(1,'  ks failed');
    fprintf(fid,'  ks failed');
  end
  fprintf(1,'%8.4f',value,p1);
  fprintf(1,'\n');
  fprintf(fid,'\n');
%
% Plot uncorrected!
%
  gc3 = figure(3);
  hold on;
  if(k <= 8)
    subplot('Position',[llw(l) llh1(1) 0.18 0.18]);
  else
    subplot('Position',[llw(l) llh1(3) 0.18 0.18]);
  end
  hold on; box on;
  plot(bereik,opp,'--k');plot(bereik,dplot,'k');
  if((k == 1) | (k == 9));ylabel('cum. prob.');end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-4 -2 0 2 4],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if (k == 1);set(gca,'YTickLabel',{'0';'';'1'});end;
  if (k == 9);set(gca,'YTickLabel',{'0';'';'1'});end;
  text(1.4,1.1,[num2str(freq(ifreq)) ' Hz'],'FontSize',8)
  set(gca,'FontSize',8);
%
% Test corrected
%
  z = runtest(reshat);
  p = normcdf(z,0,1);
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(1,'%10.4f\t%10.4f',abs(z),p);
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(fid,'%10.4f\t%10.4f',abs(z),p);
  if abs(z)<1.96
    fprintf(1,'\t run passed');
    fprintf(fid,'\t run passed');
  else
    fprintf(1,'\t run failed');
    fprintf(fid,'\t run failed');
  end
  [pass2,value2,bereik2,opp2,dplot2] = kolmogorov(reshat,k+1);
  if(value2 < foo(end,3))
      idx2 = find(value2<foo(:,3));
      idx2 = idx2(1);
  else
      idx2 = length(foo);
  end
  p2 = (100-foo(idx2,2))/100;
  fprintf(fid,'  ks:');
  fprintf(fid,' %7.4f ',value2,p2);
  if p2 >= 0.05
    fprintf(1,'\t ks passed');
    fprintf(fid,'\t ks passed');
  else
    fprintf(1,'\t ks failed');
    fprintf(fid,'\t ks failed');
  end
  fprintf(1,'%8.4f',value2,p2);
  fprintf(1,'\n');
  fprintf(fid,'\n');
%
% Plot corrected!
%
  if k+1 <= 8
    subplot('Position',[llw(l) llh1(2) 0.18 0.18]);
    l = l + 1;
  else
    subplot('Position',[llw(l) llh1(4) 0.18 0.18]);
    l = l + 1;
  end
  hold on; box on;
  plot(bereik2,opp,'--k');plot(bereik2,dplot,'k');
  if k >= 9; xlabel('Residual [\sigma]');end;
  if ((k == 1) | (k == 9));ylabel('cum. prob.');end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-4 -2 0 2 4],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if k >= 9;set(gca,'XTickLabel',{'-4';'-2';'0';'2';'4'});end;
  if k == 1;set(gca,'YTickLabel',{'0';'';'1'});end;
  if k == 9;set(gca,'YTickLabel',{'0';'';'1'});end;
  text(1.4,1.1,[num2str(freq(ifreq)) ' Hz'],'FontSize',8)
  set(gca,'FontSize',8);

  k = k+2;
end

fclose(fid);
saveas(gc1,plotfile1,'epsc2');
saveas(gc2,plotfile2,'epsc2');
saveas(gc3,plotfile3,'epsc2');


return;

% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF 

