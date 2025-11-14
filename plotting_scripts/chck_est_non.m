function [] = chck_est_non();
%
% Test real data:
% Take ML model, calc residuals, use C_d estimate to 
% decorrelate, plot histograms of that and Autocovariance.
% Then do runs and KS test with residuals.
%
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz
%load map.mat;

base = 'x_s05_1_8_32_0lay_';
data_file     = 'x_s05_1_8_32_0lay.mat';
data_file_txt = 'x_s05_1_8_32_0lay.txt';
cov_mat_file  = 'x_s05_1_8_32_0lay_cov.mat';
logfile       = 'x_s05_1_8_32_0lay_cov.log';
sdfile        = 'x_s05_1_8_32_0lay_sdev.txt';
replicafile   = 'x_s05_1_8_32_0lay_rep.dat';
residualfile  = 'x_s05_1_8_32_0lay_res.txt';
plotfile1     = 'x_s05_1_8_32_0lay_axx.png';
plotfile2     = 'x_s05_1_8_32_0lay_ks.png';
plotfile3     = 'x_s05_1_8_32_0lay_cum.png';
ibl = 1;

fid = fopen(logfile,'w');

widths = 4;
widthe = 10;

load('critical_values');
%foo = foo(:,2:end);
load(cov_mat_file);
[B] = read_r_replica(replicafile);
[C] = read_r_replica(data_file_txt);

freq =[800, 1000, 1250, 1600, 2000, 2500, 3200];
nfreq = length(freq);
clear tmp;
pksu = [0.00 0.09 0.09 0.00 0.02 0.13 0.20 0.20];
pksc = [0.04 0.07 0.20 0.00 0.00 0.03 0.20 0.16];
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
  ang = F(ifreq).ang;
  F(ifreq).res = C(ifreq).dat'-B(ifreq).dat';

%  median(res1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Compute the Autocovariance function to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [scsm1,dscsm1,sdang1,dC1] = acov_non_b(F(ifreq).res,ang,widths,widthe,ifreq);
  Arr1 = dscsm1;
%
% Decompose Cov Mat, L3 is UPPER triangular matrix:
%
  L3 = chol(F(ifreq).csave);
%
% "Uncorrelate" the residuals:
%
%  median(res2)
  F(ifreq).reshat = inv(L3')*F(ifreq).res;
%  median(reshat)
  [scsm2,dscsm2,sdang2,dC2] = acov_non_b(F(ifreq).reshat,ang,widths,widthe,ifreq);
  Arr2 = dscsm2;
%
% Runs tests:
%
  zu = runtest(F(ifreq).res);
  pu = normcdf(zu,0,1);
  zc = runtest(F(ifreq).reshat);
  pc = normcdf(zc,0,1);
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

  if(k == 1)
    text(-80,.9,['(a)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-80,-1.467,['(b)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-80,-4.367,['(a)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-80,-6.7,['(b)'],'FontSize',12,...
    'HorizontalAlignment','Right');
  end

  for itmp = 1:length(Arr1)-1
    tmp1(itmp) = Arr1(end-(itmp-1));
    tmp2(itmp) = -sdang1(end-(itmp-1));
  end
  Arr1 = [tmp1 Arr1];
  sdang1 = [tmp2 sdang1];
  plot(sdang1,Arr1/max(Arr1),'-k.');
  plot([-90 90],[0 0],'k');
  text(sdang1(end)-5,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,...
  'HorizontalAlignment','Right');
  if round(1000*pu)/1000 > 10^(-3)
    text(sdang1(end)-5,0.7,['p=' num2str(round(1000*pu)/1000)],'FontSize',8,...
    'HorizontalAlignment','Right');
  else
    text(sdang1(end)-5,0.7,['p<' num2str(10^(-3))],'FontSize',8,...
    'HorizontalAlignment','Right');
  end
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
  text(sdang2(end)-5,0.9,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,...
  'HorizontalAlignment','Right');
  if round(1000*pc)/1000 > 10^(-3)
    text(sdang1(end)-5,0.7,['p=' num2str(round(1000*pc)/1000)],'FontSize',8,...
    'HorizontalAlignment','Right');
  else
    text(sdang1(end)-5,0.7,['p<' num2str(10^(-3))],'FontSize',8,...
    'HorizontalAlignment','Right');
  end
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
  gc4 = figure(2);
  hold on;
   if(k <= 8)
     subplot('Position',[llw(jj) llh1(1) 0.18 0.18]);
   else
     subplot('Position',[llw(jj) llh1(3) 0.18 0.18]);
   end
  hold on;box on;
  plot(ang,F(ifreq).res,'-k.');
  plot([0 90],[0 0],'k');
  set(gca,'FontSize',8);

  if k+1 <= 8
    subplot('Position',[llw(jj) llh1(2) 0.18 0.18]);
    jj = jj + 1;
  else
    subplot('Position',[llw(jj) llh1(4) 0.18 0.18]);
    jj = jj + 1;
  end
  hold on;box on;
  plot(ang,F(ifreq).reshat,'-k.');
  plot([0 90],[0 0],'k');
  set(gca,'FontSize',8);


%
% Plot histograms
%
  gc2 = figure(3);
  hold on;

   if(k <= 8)
     subplot('Position',[llw(i) llh1(1) 0.18 0.18]);
   else
     subplot('Position',[llw(i) llh1(3) 0.18 0.18]);
   end
  hold on;box on;
  x = -3:1/(nang(ifreq))*20:3;
  xx = -3:.01:3;
  [n1,xout] = hist(F(ifreq).res'/ml_sd(ifreq),x);
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  area = sum(n1) * (x(2)-x(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(2.5,1.0,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,'FontSize',8,...
  'HorizontalAlignment','Right')
  set(gca,'FontSize',8);
  set(gca,'XTick',[-2 0 2]);
  set(gca,'XTickLabel',[]);
  set(gca,'YTick',[0 0.5 1]);
  set(gca,'YTickLabel',{'';'';''});
  if(k==1); set(gca,'YTickLabel',{'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'0';'';'1'});end;
  axis([-3 3 0 1.1]);

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
  [n1,xout] = hist(F(ifreq).reshat,x);
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  area = sum(n1) * (x(2)-x(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(2.5,1.0,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,'FontSize',8,...
  'HorizontalAlignment','Right')
  set(gca,'FontSize',8);
  if k>=9; xlabel('Residual [\sigma]');end;
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2';'0';'2'});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if(k==1); set(gca,'YTickLabel',{'0';'';'1'});end;
  if(k==9); set(gca,'YTickLabel',{'0';'';'1'});end;
  axis([-3 3 0 1.1]);


%
% KS Test uncorrected
%
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(1,'%10.4f\t%10.4f',abs(zu),pu);
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(fid,'%10.4f\t%10.4f',abs(zu),pu);
  if abs(zu)<1.96
    fprintf(1,'\t run passed');
    fprintf(fid,'\t run passed');
  else
    fprintf(1,'\t run failed');
    fprintf(fid,'\t run failed');
  end
  [pass,value,bereik2,opp,dplot] = kolmogorov(F(ifreq).res/ml_sd(ifreq),k);
  if(value < foo(end,3))
      idx2 = find(value<foo(:,3));
      idx2 = idx2(1);
  else
      idx2 = length(foo);
  end
  p1(ifreq) = (100-foo(idx2,2))/100;

  if p1(ifreq) >= 0.05
    fprintf(1,'\t ks passed');
%    fprintf(fid,'\t ks passed');
  else
    fprintf(1,'\t ks failed');
%    fprintf(fid,'\t ks failed');
  end
  fprintf(1,'%10.4f',p1(ifreq));
  fprintf(1,'\n');
  fprintf(fid,'\t ks %10.4f',p1(ifreq));
  fprintf(fid,'\n');
%
% Plot uncorrected!
%
  gc3 = figure(4);
  hold on;
  if(k <= 8)
    subplot('Position',[llw(l) llh1(1) 0.18 0.18]);
  else
    subplot('Position',[llw(l) llh1(3) 0.18 0.18]);
  end
  hold on;box on;

  if(k == 1)
    text(-6,0.975,['(a)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-6,-0.2,['(b)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-6,-1.66,['(a)'],'FontSize',12,...
    'HorizontalAlignment','Right');
    text(-6,-2.83,['(b)'],'FontSize',12,...
    'HorizontalAlignment','Right');
  end

  hold on; box on;
  plot(bereik2,opp,'--k');plot(bereik2,dplot,'k');
  if((k == 1) | (k == 9));ylabel('cum. prob.');end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.1]);
  set(gca,'XTick',[-4 -2 0 2 4],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if (k == 1);set(gca,'YTickLabel',{'0';'';'1'});end;
  if (k == 9);set(gca,'YTickLabel',{'0';'';'1'});end;
  text(-3.5,1.0,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,...
  'HorizontalAlignment','Left')
  if p1(ifreq) == 0.2
    text(-3.5,0.9,['p > ' num2str(p1(ifreq))],'FontSize',8,...
    'HorizontalAlignment','Left')
  elseif p1(ifreq) < 0.01
    text(-3.5,0.9,['p < ' num2str(10^(-2))],'FontSize',8,...
    'HorizontalAlignment','Left')
  else
    text(-3.5,0.9,['p = ' num2str(p1(ifreq))],'FontSize',8,...
    'HorizontalAlignment','Left')
  end


%
% KS Test corrected
%
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(1,'%10.4f\t%10.4f',abs(zc),pc);
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(fid,'%10.4f\t%10.4f',abs(zc),pc);
  if abs(zc)<1.96
    fprintf(1,'\t run passed');
    fprintf(fid,'\t run passed');
  else
    fprintf(1,'\t run failed');
    fprintf(fid,'\t run failed');
  end
  [pass2,value2,bereik2,opp,dplot] = kolmogorov(F(ifreq).reshat,k+1);
  if(value2 < foo(end,3))
      idx2 = find(value2<foo(:,3));
      idx2 = idx2(1);
  else
      idx2 = length(foo);
  end
  p2(ifreq) = (100-foo(idx2,2))/100;

  if p2(ifreq) >= 0.05
    fprintf(1,'\t ks passed');
%    fprintf(fid,'\t ks passed');
  else
    fprintf(1,'\t ks failed');
%    fprintf(fid,'\t ks failed');
  end
  fprintf(1,'%10.4f',p2(ifreq));
  fprintf(1,'\n');
  fprintf(fid,'\t ks %10.4f',p2(ifreq));
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
  set(gca,'XLim',[-4 4],'YLim',[0 1.1]);
  set(gca,'XTick',[-4 -2 0 2 4],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if k >= 9;set(gca,'XTickLabel',{'-4';'-2';'0';'2';'4'});end;
  if k == 1;set(gca,'YTickLabel',{'0';'';'1'});end;
  if k == 9;set(gca,'YTickLabel',{'0';'';'1'});end;
  text(-3.5,1.0,[num2str(freq(ifreq)) ' Hz'],'FontSize',8,...
  'HorizontalAlignment','Left')
  if p2(ifreq) == 0.2
    text(-3.5,0.9,['p > ' num2str(p2(ifreq))],'FontSize',8,...
    'HorizontalAlignment','Left')
  elseif p2(ifreq) < 0.01
    text(-3.5,0.9,['p < ' num2str(10^(-2))],'FontSize',8,...
    'HorizontalAlignment','Left')
  else
    text(-3.5,0.9,['p = ' num2str(p2(ifreq))],'FontSize',8,...
    'HorizontalAlignment','Left')
  end

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

