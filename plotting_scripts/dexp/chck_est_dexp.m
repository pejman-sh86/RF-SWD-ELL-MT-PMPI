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

sample = 'real1_sample.dat';
cov_mat_file = strrep(sample,'sample.dat','cov_mat_est_dexp.mat');
plotfile1 = strrep(sample,'sample.dat','autocorr_dexp.eps');
plotfile2 = strrep(sample,'sample.dat','res_hist_dexp.eps');
plotfile3 = strrep(sample,'sample.dat','ks_test_dexp.eps');
logfile = strrep(sample,'sample.dat','chck_log2_dexp.txt');
stdfile = strrep(sample,'sample.dat','ml_std_dexp.mat');
fid = fopen(logfile,'w');

npar = 7;
xml = [1.9486875      1.3543391      1.4835551 0.85430818       1472.925 1465.4089      0.2609109];
%xml = [2.0650185      1.3553093 1.4392162     0.85836735      1473.0071 1465.6735     0.21551576];
%xml = xmap;

load ../blde4_2_ping_ida_inv.mat;
load(cov_mat_file);
load(stdfile);
ml_sd = stdv_dB(:,2);
aincr    =  1;          %
a_start  =  1;          %
nang_tmp =135;          % Nr of angles
fincr    =  1;          %
f_start  =  6;          %
nfreq    =  8;  % Nr of frequencies
%
% Forward model parameters:
%
c1 = 1511; rho1 = 1.029;

%---------------------------------------------------------------------------
%  CODE STARTS
%---------------------------------------------------------------------------
lay_thick = 0.1;    % discretization for layer thickness - controls 
                      % upper frequency limit
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
sz=length(znorm);

%---------------------------------------------------------------------------
%  Load synthetic or real data.
%---------------------------------------------------------------------------

%
% REAL DATA
%
ang = zeros(1,nang_tmp);
freq = zeros(1,nfreq);

j = a_start;
for i = 1:nang_tmp
  ang(i) = xde4.ang(j);
  j = j + aincr;
end
j = f_start;
for i = 1:nfreq
  freq(i) = xde4.pref(j,1);
  j = j + fincr;
end
tmp = xde4.bl(f_start:fincr:f_start+(fincr*nfreq)-1,a_start:aincr:...
              a_start+(aincr*nang_tmp)-1);
for i = 1:nfreq
  jj = 1;
  for j = 1:nang_tmp
    if(isnan(tmp(i,j)) ~= 1)
      dat(jj) = tmp(i,j);
      tmp_ang(jj) = ang(j);
      jj = jj+1;
    end
  end
  F(i).dat = dat';
  nang(i) = jj-1;
  F(i).ang = tmp_ang;
end
disp('I work on real data now!');

[Frep] = forward(xml,F,freq,nfreq,nang);

for i = 1:nfreq
  F(i).diff = (Frep(i).dat-F(i).dat);
end
%
% loop over freqs 315-1600 Hz
%
ttl = [{'315 Hz'},{'400 Hz'},{'500 Hz'},{'630 Hz'},{'800 Hz'},...
       {'1000 Hz'},{'1250 Hz'},{'1600 Hz'},{'2000 Hz'}];
llw = [0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ...
       0.15 0.355 0.56 0.765 0.15 0.355 0.56 0.765 ];
llh1 = [0.77 0.58 0.34 0.15];
llh2 = [0.765 0.56 0.355 0.15];

i = 1; j = 1; k = 1; l = 1;

for ifreq = 1:nfreq
  res = F(ifreq).diff;
  Aresres = xcov(res)/nang(ifreq);
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(F(ifreq).csave);
%
% "Uncorrelate" the residuals:
%
  reshat = inv(L3')*res;
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
  if ((k==1) | (k==9)); ylabel('A_{xx}');end;
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
  set(gca,'FontSize',8);
  set(gca,'XLim',[lag2(1) lag2(end)],'YLim',[-0.5 1.1]);
  if(k>=9);xlabel('lag');end;
  if ((k==1) | (k==9)); ylabel('A_{xx}');end;
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
  [n1,xout] = hist(res/ml_sd(ifreq),x);
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
  z = runtest(res);
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(1,'%10.4f',abs(z));
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(fid,'%10.4f',abs(z));
  if abs(z)<1.96
    fprintf(1,'\t run passed');
    fprintf(fid,'\t run passed');
  else
    fprintf(1,'\t run failed');
    fprintf(fid,'\t run failed');
  end
  [pass,value,bereik2,opp,dplot] = kolmogorov_dexp(res/ml_sd(ifreq),k);
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
  fprintf(1,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(1,'%10.4f',abs(z));
  fprintf(fid,'%10.4f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(fid,'%10.4f',abs(z));
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

%===========================================================================
%   FORWARD.M
%===========================================================================

function [Fm] = forward(m,F,freq,nfreq,nang);

%--------------- GLOBAL VARIABLES ------------------------------------------
global c1 rho1 znorm lay_thick sz
%----------------------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

% rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
rhos = m(2) + sin(znorm*pi/2).^m(4)*(m(3)-m(2));
% cs=ct+(cb-ct)*znorm;
cs=m(5)+(m(6)-m(5))*znorm;

geo_sin=[NaN c1 0 rho1; ...
lay_thick*m(1)*ones(sz,1) cs' m(7)*ones(sz,1) rhos'; ...
        NaN m(6) m(7) m(3)];
%
% Compute reflectivity:
%
for i=1:nfreq
  [ref] = ref_nlay3(F(i).ang,geo_sin,freq(i));% compute Reflection
  Fm(i).dat = ref';
end 

return

% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF 

