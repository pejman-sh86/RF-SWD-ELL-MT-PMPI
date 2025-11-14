function chck_est_dexp_non
%
% Test real data:
% Take ML model, calc residuals, use C_d estimate to 
% decorrelate, plot histograms of that and Autocovariance.
% Then do runs and KS test with residuals.
%
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz
%load map.mat;

sample = 'real2_sample.dat';
resfile1 = strrep(sample,'sample.dat','residuals_dexp0.mat');
resfile2 = strrep(sample,'sample.dat','residuals_dexp2.mat');
cov_mat_file = strrep(sample,'sample.dat','cov_mat_est_dexp2.mat');
plotfile1 = strrep(sample,'sample.dat','autocorr_compare_dexp.eps');
plotfile2 = strrep(sample,'sample.dat','res_hist_dexp.eps');
plotfile3 = strrep(sample,'sample.dat','ks_test_dexp.eps');
logfile = strrep(sample,'sample.dat','chck_log2_dexp.txt');
%stdfile = strrep(sample,'sample.dat','ml_std.mat');
fid = fopen(logfile,'w');

npar = 7;
widths = 4;
widthe = 10;

load(cov_mat_file);
load(resfile2);
C = B;
load(resfile1);
%load(stdfile);
%ml_sd = stdv_dB(:,2);
aincr    =  1;          %
a_start  =  1;          %
nang_tmp =135;          % Nr of angles
fincr    =  1;          %
f_start  =  6;          %
nfreq    =  8;  % Nr of frequencies

freq = [315 400 500 630 800 1000 1250 1600];
for i = 1:nfreq
  nang(i) = length(F(i).csave);
  ml_sd(i) = sqrt(F(i).csave(1,1));
end

llw = [0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ...
       0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ];
llh1 = [0.77 0.58 0.34 0.15];
llh2 = [0.765 0.56 0.355 0.15];

i = 1; j = 1; k = 1; l = 1; jj = 1;

for ifreq = 1:nfreq
  ang = B(ifreq).ang;
  res1 = B(ifreq).diff;
  [scsm1,dscsm1,sdang1,dC1] = acov_non(res1,ang,widths,widthe,ifreq);
  Arr1 = dscsm1;
%
% Decompose Cov Mat, L3 is UPPER triangular matrix:
%
  L3 = chol(F(ifreq).csave);
%
% "Uncorrelate" the residuals:
%
  res2 = C(ifreq).diff;
  reshat = inv(L3')*res2;
  [scsm2,dscsm2,sdang2,dC2] = acov_non(reshat,ang,widths,widthe,ifreq);
  Arr2 = dscsm2;
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
  hold on;box on;
  for itmp = 1:length(Arr1)-1
    tmp1(itmp) = Arr1(end-(itmp-1));
    tmp2(itmp) = -sdang1(end-(itmp-1));
  end
  Arr1 = [tmp1 Arr1];
  sdang1 = [tmp2 sdang1];
  plot(sdang1,Arr1/max(Arr1),'-k.');
  plot([-90 90],[0 0],'k');
  text(sdang1(end)-25,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8);
%  set(gca,'FontSize',8);
  set(gca,'FontSize',10);
  set(gca,'XLim',[sdang1(1) sdang1(end)],'YLim',[-1.1 1.1]);
  if ((k==1) | (k==9)); ylabel('A_{xx}');end;
  set(gca,'XTick',[-100 -50 0 50 100],'XTickLabel',[]);
  set(gca,'YTick',[-1 -0.5 0 0.5 1],'YTickLabel',{'';'';'';'';''});
  if(k==1); set(gca,'YTickLabel',{'-1';'';'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'-1';'';'0';'';'1'});end;


  if k+1 <= 8
    subplot('Position',[llw(j) llh1(2) 0.18 0.18]);
    j = j + 1;
  else
    subplot('Position',[llw(j) llh1(4) 0.18 0.18]);
    j = j + 1;
  end
  hold on;box on;
  for itmp = 1:length(Arr2)-1
    tmp1(itmp) = Arr2(end-(itmp-1));
    tmp2(itmp) = -sdang2(end-(itmp-1));
  end
  Arr2 = [tmp1 Arr2];
  sdang2 = [tmp2 sdang2];
  plot(sdang2,Arr2/max(Arr2),'-k.');
  plot([-90 90],[0 0],'k');
  text(sdang2(end)-25,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8);
%  set(gca,'FontSize',8);
  set(gca,'FontSize',10);
  set(gca,'XLim',[sdang2(1) sdang2(end)],'YLim',[-1.1 1.1]);
  if(k>=9);xlabel('Angle diff. [deg.]');end;
  if ((k==1) | (k==9)); ylabel('A_{xx}');end;
  set(gca,'XTick',[-100 -50 0 50 100],'XTickLabel',[-100 -50 0 50 100]);
  set(gca,'YTick',[-1 -0.5 0 0.5 1],'YTickLabel',{'';'';'';'';''});
  if(k==1); set(gca,'YTickLabel',{'-1';'';'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'-1';'';'0';'';'1'});end;

%
% Plot Residuals
%
  gc4 = figure(4);
  hold on;
   if(k <= 8)
     subplot('Position',[llw(jj) llh1(1) 0.18 0.18]);
   else
     subplot('Position',[llw(jj) llh1(3) 0.18 0.18]);
   end
  hold on;box on;
  plot(ang,res1,'-k.');
  plot([0 90],[0 0],'k');
  set(gca,'FontSize',8);

  if k+1 <= 8
    subplot('Position',[llw(jj) llh1(2) 0.18 0.18]);
    jj = jj + 1;
  else
    subplot('Position',[llw(jj) llh1(4) 0.18 0.18]);
    jj = jj + 1;
  end
  hold on;
  plot(ang,reshat,'-k.');
  plot([0 90],[0 0],'k');
  set(gca,'FontSize',8);


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
  nd = 1/sqrt(2)*exp(-abs(xx));
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
  nd = 1/sqrt(2)*exp(-abs(xx));
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
  [pass,value,bereik2,opp,dplot] = kolmogorov_dexp(res1/ml_sd(ifreq),k);
  if pass == 1
    fprintf(1,'\t ks passed');
%    fprintf(fid,'\t ks passed');
  else
    fprintf(1,'\t ks failed');
%    fprintf(fid,'\t ks failed');
  end
  fprintf(1,'%10.4f',value);
  fprintf(1,'\n');
  fprintf(fid,'\t ks %10.4f',value);
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
  plot(bereik2,opp,'--k');plot(bereik2,dplot,'k');
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
  [pass2,value2,bereik2,opp,dplot] = kolmogorov_dexp(reshat,k+1);
  if pass2 == 1
    fprintf(1,'\t ks passed');
%    fprintf(fid,'\t ks passed');
  else
    fprintf(1,'\t ks failed');
%    fprintf(fid,'\t ks failed');
  end
  fprintf(1,'%10.4f',value2);
  fprintf(1,'\n');
  fprintf(fid,'\t ks %10.4f',value2);
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

