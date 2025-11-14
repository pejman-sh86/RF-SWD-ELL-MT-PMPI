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

data          = '2_b_80m_y_m0.txt';
cov_mat_file  = '2_b_80m_y_m0_2lay_cov.txt';
cov_mat_file2 = '2_b_80m_y_m0_2lay_cov.mat';
logfile       = '2_b_80m_y_m0_2lay_cov.log';
sdfile        = '2_b_80m_y_m0_2lay_sdev.txt';
residualfile  = '2_b_80m_y_m0_2lay_res.txt';
plotfile1     = '2_b_80m_y_m0_2lay_axx.png';
plotfile2     = '2_b_80m_y_m0_2lay_ks.png';
plotfile3     = '2_b_80m_y_m0_2lay_cum.png';
ibl = 1;

fid = fopen(logfile,'w');

res = load(residualfile);
load('critical_values.mat');
sd2 = load(sdfile);
dat = load(data);

freq = dat(2:end,:);
M = length(F(1).minlim);
widths = 4.;
widthe = 20.;

%
% loop over freqs 315-1600 Hz
%
k = 1;

nx = 6;
ny = 4;
xim = 0.01;
yim = 0.1/ny;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
nx2 = 4;
ny2 = 2;
[loc2,spw2,sph2] = get_loc(nx2,ny2,xim,yim,xymarg);

for ifreq = 1:NBAND

  N = length(res(ifreq).dat);
  if(ifreq == 1)
    disp('    M,  N,   NBAND, M/NBAND: ');
    fprintf(1,'%5i, %f10.4',M,N,NBAND,M/NBAND);
    fprintf(1,'\n');
  end;
  sd(ifreq) = mean(sd2(ifreq).dat);
  res(ifreq).dat = res(ifreq).dat';
  sd2(ifreq).dat = sd2(ifreq).dat';
  res(ifreq).dat = res(ifreq).dat - mean(res(ifreq).dat);

  if(istat == 1)
    res2(ifreq).dat = res(ifreq).dat./sd2(ifreq).dat;
  else
    res2 = res;
  end

  N2 = N-(M/NBAND);
%
% Zero padding:
%
%  npad = floor(N/6)
  npad = 1;
%
% Calc autocovariance:
%
  [scnew,dscnew,sdang,dC,sc,rscnew,rdscnew,sd3,x,C,fact] = ...
  acov_non(res(ifreq).dat,res(ifreq).ang,widths,widthe,ifreq);

  F(ifreq).csave = dC;
  L = chol(dC);
  Aresres = dscnew;
  Chat = dC;
  if(istat == 1)
      for i = 1:length(Chat)
      for j = 1:length(Chat)
          Chat(i,j) = Chat(i,j)*sd2(ifreq).dat(i)*sd2(ifreq).dat(j);
      end;end;
  end;
  F(ifreq).csave = Chat;
%
% Decompose that thing, L3 is UPPER triangular matrix:
%
  L3 = chol(Chat);

%%
%% "Uncorrelate" the residuals:
%%
  F(ifreq).reshat = inv(L3')*res(ifreq).dat;
  F(ifreq).res = res(ifreq).dat;

  [scnew2,dscnew2,sdang2,dC2,sc2,rscnew2,rdscnew2,sd4,x2,C2,fact2] = ...
  acov_non(F(ifreq).reshat,res(ifreq).ang,widths,widthe,ifreq);

  Aresres2 = dscnew2;

%%
%%  RESIDUALS 
%%
  gc0 = figure(1);hold on;
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);hold on;box on;
  plot(res(ifreq).ang,res(ifreq).dat,'-k.')
  plot([0 90],[0 0],'--k')
  text(40,2.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-3 -2 -1. 0 1. 2. 3.],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==nx+1) | (k==(2*nx)+1)); 
      ylabel('Data Res.','FontSize',14);
      set(gca,'YTickLabel',[-3. -2. -1 0 1. 2. 3.],'FontSize',14);
  end;
  set(gca,'XTick',[20 40 60],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('Angle (deg.)','FontSize',14);
     set(gca,'XTick',[20 40 60],'FontSize',14);
  end;
  set(gca,'XLim',[10 60]);
  set(gca,'YLim',[-4 4]);
  if(ibl == 0)
     set(gca,'YLim',[-1 1]);
  end
  set(gca,'FontSize',14);

  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);hold on;box on;
  plot(res(ifreq).ang,F(ifreq).reshat,'-k.')
  plot([0 90],[0 0],'--k')
  text(40,2.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-3 -2 -1. 0 1. 2. 3.],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  set(gca,'XTick',[20 40 60],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('Angle (deg.)','FontSize',14);
     set(gca,'XTick',[20 40 60],'FontSize',14);
  end;
  set(gca,'XLim',[10 60]);
  set(gca,'YLim',[-4 4]);
  if(ibl == 0)
     set(gca,'YLim',[-1 1]);
  end;
  set(gca,'FontSize',14);

%%
%%  AUTOCOVARIANCE
%%
  gc1 = figure(2);
  hold on;
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);hold on;box on;

  plot([-1.*fliplr(sdang),sdang],[fliplr(Aresres),Aresres]/...
       (max(Aresres)),'-k.')
  plot([-60 60],[0 0],'--k')
  text(20,.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-0.5 0 0.5 1],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==nx+1) | (k==(2*nx)+1)); 
      ylabel('A_{xx}','FontSize',14);
      set(gca,'YTickLabel',[-0.5 0 0.5 1],'FontSize',14);
  end;
  set(gca,'XTick',[-40 -20 0 20 40],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-40 -20 0 20 40],'FontSize',14);
  end;
  set(gca,'XLim',[-60 60]);
  set(gca,'YLim',[-0.5 1.2]);
  set(gca,'FontSize',14);

  
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);hold on;box on;
  plot([-1.*fliplr(sdang2),sdang2],[fliplr(Aresres2),Aresres2]/...
       (max(Aresres2)),'-k.')
  plot([-60 60],[0 0],'--k')
  text(20,0.8,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'YTick',[-0.5 0 0.5 1],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  set(gca,'XTick',[-40 -20 0 20 40],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>(2*length(freq))-nx);
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-40 -20 0 20 40],'FontSize',14);
  end;
  set(gca,'XLim',[-60 60]);
  set(gca,'YLim',[-0.5 1.2]);
  set(gca,'FontSize',14);

  x = -9.375:.75:9.375;
  xx = -9.25:.01:9.25;

%%
%%  HISTROGRAMS
%%
  gc2 = figure(3);
  hold on;
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);
  hold on;box on;
  [n1,xout] = hist(res(ifreq).dat./sd2(ifreq).dat,x);
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  plot(xx,nd,'--k')
  text(2.5,0.55,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
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
  [n1,xout] = hist(F(ifreq).reshat,x);
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
%  n1 = n1/(N2*(xout(2)-xout(1)));
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(2.5,0.55,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  axis([-6 6 0 0.6]);
  set(gca,'FontSize',14);

  z = runtest(res(ifreq).dat);
  p = normcdf(z,0,1);
%  p = 0;
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

  [pass,value,bereik,opp,dplot] = kolmogorov(res(ifreq).dat./sd2(ifreq).dat,k);
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
    fprintf(1,'\t ks passed');
    fprintf(fid,' ks passed');
  else
    fprintf(1,'\t ks failed');
    fprintf(fid,'  ks failed');
  end
  fprintf(1,'%8.4f',value,p1);
  fprintf(1,'\n');
  fprintf(fid,'\n');

  z = runtest(F(ifreq).reshat);
  p = normcdf(z,0,1);
%  p = 0;
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

  [pass2,value2,bereik2,opp2,dplot2] = kolmogorov(F(ifreq).reshat,k+1);
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
    fprintf(fid,'  ks passed');
  else
    fprintf(1,'\t ks failed');
    fprintf(fid,'  ks failed');
  end
  fprintf(1,'%8.4f',value2,p2);
  fprintf(1,'\n');
  fprintf(1,'\n');
  fprintf(fid,'\n');

%%
%%  CUMULATIVE DISTRIBUTIONS
%%
  gc3 = figure(4);
  hold on;
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);
  hold on; box on;
  plot(bereik,opp(1:length(bereik)),'--k');
  plot(bereik,dplot(1:length(bereik)),'k');
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
  if (k==1 | k==nx+1 | k==2*nx+1 | k==3*nx+1);
     set(gca,'YTickLabel',{'0';'';'1'});
  end;
  if(k > (2*length(freq))-nx);set(gca,'XTickLabel',{'-2';'0';'2'});end;
  text(-3.5,1.,[num2str(freq(ifreq)) ' Hz'],'FontSize',12)
  set(gca,'FontSize',14);

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


%%
%%  COVARIANCE MATRICES
%%
  gc4 = figure(5);
  hold on;
  subplot('Position',[loc2(1,ifreq) loc2(2,ifreq) spw2 sph2]);
  pcolor(res(ifreq).ang,res(ifreq).ang,Chat);shading flat;grid off;
  set(gca,'layer','top')
  subplot('Position',[loc2(1,ifreq) loc2(2,ifreq) spw2 sph2]);
  hold on; box on;
  set(gca,'YDir','reverse');
  set(gca,'XTick',[20 40 60],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  set(gca,'YTick',[20 40 60],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  set(gca,'XLim',[res(ifreq).ang(1) res(ifreq).ang(end)]);
  set(gca,'YLim',[res(ifreq).ang(1) res(ifreq).ang(end)]);
  if ((ifreq==1) | (ifreq==nx2+1) | (ifreq==(2*nx2)+1)); 
     ylabel('Angle (deg.)','FontSize',14);
     set(gca,'YTickLabel',[20 40 60],'FontSize',14);
  end;
  if (ifreq > NBAND-nx2); 
     xlabel('Angle (deg.)','FontSize',14);
     set(gca,'XTickLabel',[20 40 60],'FontSize',14);
  end;

  k = k+2;
end

save tmp res F;

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
    for ifreq = 2:NBAND

        covmat = inv(F(ifreq).csave);
        save(cov_mat_file,'covmat','-ascii','-append');

    end;

    for ifreq = 1:NBAND

        covmat = F(ifreq).csave;
        save(cov_mat_file,'covmat','-ascii','-append');

    end;

end;

return;

