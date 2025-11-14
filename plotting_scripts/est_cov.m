function [] = est_cov_mat();
set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);
%
% Test real data:
%
% UNCOMMENT SAVE STATEMENT!!!!
%
dex = 0;
fact = 1;
istat = 1;
isave  = 1;
isavep = 1;
nx = 8;
ny = 6;

%data          = 'x_s21_1_5_40_6lay.mat';
%cov_mat_file  = 'x_s21_1_5_40_6laycov_cov.txt';
%cov_mat_file2 = 'x_s21_1_5_40_6laycov_cov.mat';
%logfile       = 'x_s21_1_5_40_6laycov_cov.log';
%sdfile        = 'x_s21_1_5_40_6laycov_sdev.txt';
%residualfile  = 'x_s21_1_5_40_6laycov_res.txt';
%plotfile1     = 'x_s21_1_5_40_6laycov_axx.png';
%plotfile2     = 'x_s21_1_5_40_6laycov_ks.png';
%plotfile3     = 'x_s21_1_5_40_6laycov_cum.png';
%ibl = 0;

%data          = 'x_s20_1_8_40_0lay.mat';
%cov_mat_file  = 'x_s20_1_8_40_0lay_cov.txt';
%cov_mat_file2 = 'x_s20_1_8_40_0lay_cov.mat';
%logfile       = 'x_s20_1_8_40_0lay_cov.log';
%sdfile        = 'x_s20_1_8_40_0lay_sdev.txt';
%residualfile  = 'x_s20_1_8_40_0lay_res.txt';
%plotfile1     = 'x_s20_1_8_40_0lay_axx.png';
%plotfile2     = 'x_s20_1_8_40_0lay_ks.png';
%plotfile3     = 'x_s20_1_8_40_0lay_cum.png';
%ibl = 1;

%data          = 'x_s19_1_6_25_0lay.mat';
%cov_mat_file  = 'x_s19_1_6_25_0lay_c1_cov.txt';
%cov_mat_file2 = 'x_s19_1_6_25_0lay_c1_cov.mat';
%logfile       = 'x_s19_1_6_25_0lay_c1_cov.log';
%sdfile        = 'x_s19_1_6_25_0lay_c1_sdev.txt';
%residualfile  = 'x_s19_1_6_25_0lay_c1_res.txt';
%plotfile1     = 'x_s19_1_6_25_0lay_c1_axx.png';
%plotfile2     = 'x_s19_1_6_25_0lay_c1_ks.png';
%plotfile3     = 'x_s19_1_6_25_0lay_c1_cum.png';
%ibl = 1;

%data          = 'x_s16_1_5_25_7lay.mat';
%cov_mat_file  = 'x_s16_1_5_25_7laycov_cov.txt';
%cov_mat_file2 = 'x_s16_1_5_25_7laycov_cov.mat';
%logfile       = 'x_s16_1_5_25_7laycov_cov.log';
%sdfile        = 'x_s16_1_5_25_7laycov_sdev.txt';
%residualfile  = 'x_s16_1_5_25_7laycov_res.txt';
%plotfile1     = 'x_s16_1_5_25_7laycov_axx.png';
%plotfile2     = 'x_s16_1_5_25_7laycov_ks.png';
%plotfile3     = 'x_s16_1_5_25_7laycov_cum.png';
%ibl = 0;

%data          = 'x_s13_2_10_50_5layb.mat';
%cov_mat_file  = 'x_s13_2_10_50_5layb_cov.txt';
%cov_mat_file2 = 'x_s13_2_10_50_5layb_cov.mat';
%logfile       = 'x_s13_2_10_50_5layb_cov.log';
%sdfile        = 'x_s13_2_10_50_5layb_sdev.txt';
%residualfile  = 'x_s13_2_10_50_5layb_res.txt';
%plotfile1     = 'x_s13_2_10_50_5layb_axx.png';
%plotfile2     = 'x_s13_2_10_50_5layb_ks.png';
%plotfile3     = 'x_s13_2_10_50_5layb_cum.png';
%ibl = 0;

data          = 'x_s07_1_1_30_nodisp.mat';
cov_mat_file  = 'x_s07_1_1_30_nodisp_cov.txt';
cov_mat_file2 = 'x_s07_1_1_30_nodisp_cov.mat';
logfile       = 'x_s07_1_1_30_nodisp_cov.log';
sdfile        = 'x_s07_1_1_30_nodisp_sdev.txt';
residualfile  = 'x_s07_1_1_30_nodisp_res.txt';
plotfile1     = 'x_s07_1_1_30_nodisp_axx.png';
plotfile2     = 'x_s07_1_1_30_nodisp_ks.png';
plotfile3     = 'x_s07_1_1_30_nodisp_cum.png';
ibl = 0;

%data          = 'x_s07_1_1_100_alf_2lay.mat';
%cov_mat_file  = 'x_s07_1_1_100_alf_2lay_cov.txt';
%cov_mat_file2 = 'x_s07_1_1_100_alf_2lay_cov.mat';
%logfile       = 'x_s07_1_1_100_alf_2lay_cov.log';
%sdfile        = 'x_s07_1_1_100_alf_2lay_sdev.txt';
%residualfile  = 'x_s07_1_1_100_alf_2lay_res.txt';
%plotfile1     = 'x_s07_1_1_100_alf_2lay_axx.png';
%plotfile2     = 'x_s07_1_1_100_alf_2lay_ks.png';
%plotfile3     = 'x_s07_1_1_100_alf_2lay_cum.png';
%ibl = 1;

%data          = 'x_s05_1_8_32_0lay.mat';
%cov_mat_file  = 'x_s05_1_8_32_0lay_cov.txt';
%cov_mat_file2 = 'x_s05_1_8_32_0lay_cov.mat';
%logfile       = 'x_s05_1_8_32_0lay_cov.log';
%sdfile        = 'x_s05_1_8_32_0lay_sdev.txt';
%residualfile  = 'x_s05_1_8_32_0lay_res.txt';
%plotfile1     = 'x_s05_1_8_32_0lay_axx.png';
%plotfile2     = 'x_s05_1_8_32_0lay_ks.png';
%plotfile3     = 'x_s05_1_8_32_0lay_cum.png';
%ibl = 1;

%data          = 'x_s05_2_8_25_4lay.mat';
%cov_mat_file  = 'x_s05_2_8_25_4laycov_cov.txt';
%cov_mat_file2 = 'x_s05_2_8_25_4laycov_cov.mat';
%logfile       = 'x_s05_2_8_25_4laycov_cov.log';
%sdfile        = 'x_s05_2_8_25_4laycov_sdev.txt';
%residualfile  = 'x_s05_2_8_25_4laycov_res.txt';
%plotfile1     = 'x_s05_2_8_25_4laycov_axx.png';
%plotfile2     = 'x_s05_2_8_25_4laycov_ks.png';
%plotfile3     = 'x_s05_2_8_25_4laycov_cum.png';
%ibl = 0;

%data          = 'x_s04_1_3_32_0lay.mat';
%cov_mat_file  = 'x_s04_1_3_32_0lay_cov.txt';
%cov_mat_file2 = 'x_s04_1_3_32_0lay_cov.mat';
%logfile       = 'x_s04_1_3_32_0lay_cov.log';
%sdfile        = 'x_s04_1_3_32_0lay_sdev.txt';
%residualfile  = 'x_s04_1_3_32_0lay_res.txt';
%plotfile1     = 'x_s04_1_3_32_0lay_axx.png';
%plotfile2     = 'x_s04_1_3_32_0lay_ks.png';
%plotfile3     = 'x_s04_1_3_32_0lay_cum.png';
%ibl = 1;

%data          = 'x_s02_1_3_25_7lay.mat';
%cov_mat_file  = 'x_s02_1_3_25_7laycov_cov.txt';
%cov_mat_file2 = 'x_s02_1_3_25_7laycov_cov.mat';
%logfile       = 'x_s02_1_3_25_7laycov_cov.log';
%sdfile        = 'x_s02_1_3_25_7laycov_sdev.txt';
%residualfile  = 'x_s02_1_3_25_7laycov_res.txt';
%plotfile1     = 'x_s02_1_3_25_7laycov_axx.png';
%plotfile2     = 'x_s02_1_3_25_7laycov_ks.png';
%plotfile3     = 'x_s02_1_3_25_7laycov_cum.png';
%ibl = 0;

%data          = 'x_s01_1_3_20_5lay.mat';
%cov_mat_file  = 'x_s01_1_3_20_5lay_cov.txt';
%cov_mat_file2 = 'x_s01_1_3_20_5lay_cov.mat';
%logfile       = 'x_s01_1_3_20_5lay_cov.log';
%sdfile        = 'x_s01_1_3_20_5lay_sdev.txt';
%residualfile  = 'x_s01_1_3_20_5lay_res.txt';
%plotfile1     = 'x_s01_1_3_20_5lay_axx.png';
%plotfile2     = 'x_s01_1_3_20_5lay_ks.png';
%plotfile3     = 'x_s01_1_3_20_5lay_cum.png';
%ibl = 0;

%data          = 'sim_T_1_4_17_0lay.mat';
%cov_mat_file  = 'sim_T_1_4_17_0lay_cov.txt';
%cov_mat_file2 = 'sim_T_1_4_17_0lay_cov.mat';
%logfile       = 'sim_T_1_4_17_0lay_cov.log';
%sdfile        = 'sim_T_1_4_17_0lay_sdev.txt';
%residualfile  = 'sim_T_1_4_17_0lay_res.txt';
%plotfile1     = 'sim_T_1_4_17_0lay_axx.png';
%plotfile2     = 'sim_T_1_4_17_0lay_ks.png';
%plotfile3     = 'sim_T_1_4_17_0lay_cum.png';
%ibl = 0;

%data          = 'sim_C_1_7_40_5layrg.mat';
%cov_mat_file  = 'sim_C_1_7_40_5layrg_cov.txt';
%cov_mat_file2 = 'sim_C_1_7_40_5layrg_cov.mat';
%logfile       = 'sim_C_1_7_40_5layrg_cov.log';
%sdfile        = 'sim_C_1_7_40_5layrg_sdev.txt';
%residualfile  = 'sim_C_1_7_40_5layrg_res.txt';
%plotfile1     = 'sim_C_1_7_40_5layrg_axx.png';
%plotfile2     = 'sim_C_1_7_40_5layrg_ks.png';
%plotfile3     = 'sim_C_1_7_40_5layrg_cum.png';
%ibl = 0;

%data          = 'sim_B_1_10_50_7lay.mat';
%cov_mat_file  = 'sim_B_1_10_50_7lay_cov.txt';
%cov_mat_file2 = 'sim_B_1_10_50_7lay_cov.mat';
%logfile       = 'sim_B_1_10_50_7lay_cov.log';
%sdfile        = 'sim_B_1_10_50_7lay_sdev.txt';
%residualfile  = 'sim_B_1_10_50_7lay_res.txt';
%plotfile1     = 'sim_B_1_10_50_7lay_axx.png';
%plotfile2     = 'sim_B_1_10_50_7lay_ks.png';
%plotfile3     = 'sim_B_1_10_50_7lay_cum.png';
%ibl = 0;

%data          = 'sim_A_4.mat';
%cov_mat_file  = 'sim_A_4_ci_cov.txt';
%cov_mat_file2 = 'sim_A_4_ci_cov.mat';
%logfile       = 'sim_A_4_ci_cov.log';
%sdfile        = 'sim_A_4_sdev.txt';
%residualfile  = 'sim_A_4_res.txt';
%ibl = 0;

fid = fopen(logfile,'w');

load(residualfile);
load('critical_values.mat');

sd2 = load(sdfile);
sd = mean(sd2,1);
res = load(residualfile);
res = res';
res(:,1:end-1) = res(:,1:end-1) - repmat(mean(res(:,1:end-1),1),size(res,1),1);

sd2 = sd2';
load(data);
nfreq = length(F(1).freq(:));
freq = F(1).freq;

M = length(F(1).minlim);
N = size(res,1);
Nf = length(freq);
disp('    M,  N,   Nf, M/Nf: '),fprintf(1,'%5i, %f10.4',M,N,Nf,M/Nf);

if(istat == 1)
    res2 = res(:,1:nfreq)./sd2;
else
    res2 = res;
end

%
% loop over freqs 315-1600 Hz
%
k = 1;
res(1:N,1:nfreq) = res(1:N,1:nfreq).*F(1).Rex(1:N,1:nfreq);

%ny = ceil(2*nfreq/nx)
xim = 0.01;
yim = 0.1/ny;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

for ifreq = 1:nfreq

  idx = find(F(1).Rex(1:N,ifreq) == 1);
  N2 = length(idx)-(M/Nf);
  disp(N2);
%
% Zero padding:
%
%  npad = floor(N/6)
  npad = 1;
%
% Calc autocovariance:
%
  [Aresres,lag1] = xcov(res2(1:N,ifreq));
  Aresres = Aresres/(N2);
%------------------------------------------------
%  Damping
%------------------------------------------------
  x = cos(pi*[-N:N]'./(2*N)).^fact;
  gc0 = figure(1);
  plot(x);
  Aresres = Aresres .* x(2:end-1);

%
% Set up covariance matrix in a diagonal form:
%
  sd3 = zeros(N,(2*N)-1);
  for i = 1:(2*N)-1
    sd3(:,i) = Aresres(i)*ones(N,1);
  end
  sd3(:,1:npad) = 0;
  sd3(:,(2*N)-npad:(2*N)-1) = 0;
  E = diag(ones(N,1));
  Chat = zeros(N,N);
  Chat = spdiags(sd3,-N+1:N-1,Chat);

  if(istat == 1)
      for i = 1:length(Chat)
      for j = 1:length(Chat)
          Chat(i,j) = Chat(i,j)*sd2(i,ifreq)*sd2(j,ifreq);
      end;end;
  end;
  Csave = E*Chat;

  F(ifreq).csave = Csave;
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(Chat);
%
% "Uncorrelate" the residuals:
%
  reshat(:,ifreq) = inv(L3')*res(1:N,ifreq);
  [Aresres2,lag2] = xcov(reshat(:,ifreq));
  Aresres2 = Aresres2/N2;
  
  gc1 = figure(2);
  hold on;
%  subplot(4,4,k)
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);

  plot(lag1,Aresres/(max(Aresres)),'-k.')
  text(20,.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-0.5 0 0.5 1],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==nx+1) | (k==(2*nx)+1)); 
      ylabel('A_{xx}','FontSize',14);
      set(gca,'YTickLabel',[-0.5 0 0.5 1],'FontSize',14);
  end;
  set(gca,'XTick',[-100 0 100],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-100 0 100],'FontSize',14);
  end;
  if(ibl == 0)
     set(gca,'YLim',[-0.5 1.2]);
  end
  set(gca,'XLim',[-N N]);
  set(gca,'FontSize',14);

%  subplot(4,4,k+1)
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);

  plot(lag2,Aresres2/(max(Aresres2)),'-k.')
  text(20,0.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-0.5 0 0.5 1],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  set(gca,'XTick',[-100 0 100],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-100 0 100],'FontSize',14);
  end;
  if(ibl == 0)
     set(gca,'YLim',[-0.5 1.2]);
  end;
  set(gca,'XLim',[-N N]);
  set(gca,'FontSize',14);

  x = -6.375:.75:6.375;
  xx = -6.25:.01:6.25;
  gc2 = figure(3);
  hold on;

  subplot('Position',[loc(1,k) loc(2,k) spw sph]);

  hold on;box on;
  [n1,xout] = hist(res(idx,ifreq)./sd2(idx,ifreq),x);
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  plot(xx,nd,'--k')
  text(1,0.55,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
%  if (k>12 | k == 11 | k == 12); xlabel('Res. (std. dev.)','FontSize',14);end;
%  if ((k==1) | (k==5) | (k==9) | (k==13)); ylabel('','FontSize',14);end;
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  ylabel('','FontSize',14);
  set(gca,'YTick',[0 0.5],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==nx+1) | (k==2*nx+1) | (k==3*nx+1)); 
      set(gca,'YTickLabel',[0 0.5],'FontSize',14);
  end;
  axis([-6 6 0 0.6]);
  set(gca,'FontSize',14);

  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);

  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  set(gca,'YTick',[0 0.5],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  hold on;box on;
%  [n1,xout] = hist(reshat(idx,ifreq),x);
  [n1,xout] = hist(reshat(idx,ifreq),x);
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
%  n1 = n1/(N2*(xout(2)-xout(1)));
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(1,0.55,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  axis([-6 6 0 0.6]);
  set(gca,'FontSize',14);

  if(length(res(idx,ifreq))>20)
    z = runtest(res(idx,ifreq));
    p = normcdf(z,0,1);
  else
    z = 4.;
    p = 0;
  end;
  fprintf(1,'%2.0f',ifreq);fprintf(1,'\t uncorrected:\t');
  fprintf(1,'%7.4f',abs(z));
  fprintf(fid,'%2.0f',ifreq);
  fprintf(fid,'%7.4f  %7.4f',abs(z),p);
  if abs(z)<1.96
    fprintf(1,'\t passed');
    fprintf(fid,'  run passed');
  else
    fprintf(1,'\t failed');
    fprintf(fid,'  run failed');
  end
  fprintf(1,'\n');

  if(length(res(idx,ifreq)) > 20)
    [pass,value,bereik,opp,dplot] = kolmogorov(res(idx,ifreq)./sd(ifreq),k);
    if(value < foo(end,3))
        idx2 = find(value<foo(:,3));
        idx2 = idx2(1);
    else
        idx2 = length(foo);
    end
    p1 = (100-foo(idx2,2))/100;
  else
    bereik = zeros(size(idx));
    opp    = zeros(size(idx));
    dplot  = zeros(size(idx));
    value = 4.;
    p1 = 1.;
  end
  fprintf(fid,'   ks:');
  fprintf(fid,' %7.4f ',value,p1);
  if p1 >= 0.05
    fprintf(1,'\t ks passed');
    fprintf(fid,' ks passed');
  else
    fprintf(1,'\t ks failed');
    fprintf(fid,'  ks failed');
  end
  fprintf(1,'%8.4f',value,p1);
  fprintf(1,'\n');
  fprintf(fid,'\n');

  if(length(reshat(idx,ifreq))>20)
    z = runtest(reshat(idx,ifreq));
    p = normcdf(z,0,1);
  else
    z = 4.;
    p = 0;
  end;
  fprintf(1,'%2.0f',ifreq);fprintf(1,'\t corrected:\t');
  fprintf(1,'%6.4f',abs(z));
  fprintf(fid,'%2.0f',ifreq);
  fprintf(fid,'%7.4f  %7.4f',abs(z),p);
  if abs(z)<1.96
    fprintf(1,'\t passed');
    fprintf(fid,'  run passed');
  else
    fprintf(1,'\t failed');
    fprintf(fid,'  run failed');
  end
  fprintf(1,'\n');

  if(length(reshat(idx,ifreq)) > 20)
    [pass2,value2,bereik2,opp2,dplot2] = kolmogorov(reshat(idx,ifreq),k+1);
    if(value2 < foo(end,3))
        idx2 = find(value2<foo(:,3));
        idx2 = idx2(1);
    else
        idx2 = length(foo);
    end
    p2 = (100-foo(idx2,2))/100;
  else
    bereik2 = zeros(size(idx));
    opp2   = zeros(size(idx));
    dplot2 = zeros(size(idx));
    value2 = 4.;
    p2 = 1.;
  end;
  fprintf(fid,'  ks:');
  fprintf(fid,' %7.4f ',value2,p2);
  if p2 >= 0.05
    fprintf(1,'\t ks passed');
    fprintf(fid,'  ks passed');
  else
    fprintf(1,'\t ks failed');
    fprintf(fid,'  ks failed');
  end
  fprintf(1,'%8.4f',value2,p2);
  fprintf(1,'\n');
  fprintf(1,'\n');
  fprintf(fid,'\n');

%
% Plot uncorrected!
%
  gc3 = figure(4);
  hold on;
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);
  hold on; box on;
  plot(bereik,opp(1:length(bereik)),'--k');plot(bereik,dplot(1:length(bereik)),'k');
  if ((k==1) | (k==nx+1) | (k==(2*nx)+1)); 
     ylabel('cum. prob.','FontSize',14);
  end;
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k > (2*length(freq))-nx); 
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  set(gca,'YTickLabel',[]);
  if (k==1 | k==nx+1 | k==2*nx+1 | k==3*nx+1);set(gca,'YTickLabel',{'0';'';'1'});end;
  if(k > (2*length(freq))-nx);set(gca,'XTickLabel',{'-2';'0';'2'});end;
  text(-3.5,1.,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'FontSize',14);

%
% Plot corrected!
%
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);
  hold on; box on;
  plot(bereik2(1:size(bereik2,1)),opp2(1:size(bereik2,1)),'--k');
  plot(bereik2,dplot2,'k');
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k > (2*length(freq))-nx); 
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
  if(k > (2*length(freq))-nx);set(gca,'XTickLabel',{'-2';'0';'2'});end;
  text(-3.5,1.,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'FontSize',14);

  k = k+2;

end

save tmp reshat res F;

if(isavep == 1)
%    saveas(gc1,plotfile1,'epsc2');
%    saveas(gc2,plotfile2,'epsc2');
%    saveas(gc3,plotfile3,'epsc2');
    saveas(gc1,plotfile1,'png');
    saveas(gc2,plotfile2,'png');
    saveas(gc3,plotfile3,'png');
end

%
%  Saves both, inverse and orgiginal Cov mat.
%
if(isave == 1)
    save(cov_mat_file2,'F');

    covmat = inv(F(1).csave);
    save(cov_mat_file,'covmat','-ascii');
    for ifreq = 2:nfreq

        covmat = inv(F(ifreq).csave);
        save(cov_mat_file,'covmat','-ascii','-append');

    end;

    for ifreq = 1:nfreq

        covmat = F(ifreq).csave;
        save(cov_mat_file,'covmat','-ascii','-append');

    end;

end;

return;
