function [] = est_cov_mat();
set(0, 'DefaultFigurePaperPosition', [0 0 8 6]);
%
% Test real data:
%
% UNCOMMENT SAVE STATEMENT!!!!
%
dex = 0;
M = 10;
N = 83;
nfreq = 4;
fact = 1;
istat = 1;

cov_mat_file = 'cov_offstat.txt';
cov_mat_file2 = 'cov_offstat.mat';
residualfile = 'residuals_map.mat';
sdfile = 'sigma_map.txt';
logfile = 'cov_offstat.log';
fid = fopen(logfile,'w');

load(residualfile);
load('critical_values.mat');
res = res';
sd = std(res)
sd2 = load(sdfile);
sd2 = sd2';

%
% loop over freqs 315-1600 Hz
%
ttl = [{'R1'},{'R2'},{'R3'},...
       {'R4'},{'R5'},{'R6'},{'R7'}];
k = 1;
if(istat == 1)
    res2 = res./sd2;
else
    res2 = res;
end

nx = 4
ny = 4
xim = 0.02;
yim = 0.02;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

for ifreq = 1:1:nfreq

  N2 = N-M;
%
% Zero padding:
%
%  npad = floor(N/6)
  npad = 1;
%
% Calc autocovariance:
%
  [Aresres,lag1] = xcov(res2(:,ifreq));
  Aresres = Aresres/(N2);
%------------------------------------------------
%  Damping
%------------------------------------------------
  x = cos(pi*[-N:N]'./(2*N)).^fact;
  figure(1);
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
  reshat = inv(L3')*res(:,ifreq);
  [Aresres2,lag2] = xcov(reshat);
  Aresres2 = Aresres2/N2;
  
  gc1 = figure(2);
  hold on;
%  subplot(4,4,k)
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);

  plot(lag1,Aresres/(max(Aresres)),'-k.')
  text(50,1,ttl(ifreq))
%  if (k>12 | k == 11 | k == 12); xlabel('lag','FontSize',14);end;
%  if ((k==1) | (k==5) | (k==9) | (k==13)); ylabel('Axx','FontSize',14);end;
  if ((k==5)); ylabel('Axx','FontSize',14);end;
  set(gca,'XTick',[-50 0 50],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>10); 
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-50 0 50],'FontSize',14);
  end;
  set(gca,'YLim',[-0.5 1.2]);
  set(gca,'XLim',[-N N]);
  set(gca,'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==5) | (k==9) | (k==13)); 
     set(gca,'YTickLabel',[-0.5 0.0 0.5 1.0],'FontSize',14);
  end;
%  subplot(4,4,k+1)
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);
  plot(lag2,Aresres2/(max(Aresres2)),'-k.')
  text(50,1,ttl(ifreq))
%  if (k>12 | k == 11 | k == 12); xlabel('lag','FontSize',14);end;
  set(gca,'XTick',[-50 0 50],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>10); 
     xlabel('lag','FontSize',14);
     set(gca,'XTickLabel',[-50 0 50],'FontSize',14);
  end;
  set(gca,'YLim',[-0.5 1.2]);
  set(gca,'XLim',[-N N]);
  set(gca,'YTickLabel',[],'FontSize',14);
  set(gca,'FontSize',14);

  x = -5:.5:5;
  xx = -5:.01:5;
  gc2 = figure(3);
  hold on;
%  subplot(4,4,k);
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);
  hold on;box on;
  [n1,xout] = hist(res(:,ifreq)./sd2(:,ifreq),x);
%  n1 = n1/(N2*(xout(2)-xout(1)));
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
  stairs(xout,n1,'k');
  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  plot(xx,nd,'--k')
  text(3.5,.6,ttl(ifreq));
%  if ((k==1) | (k==5) | (k==9) | (k==13)); ylabel('','FontSize',14);end;
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>10); 
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  ylabel('','FontSize',14);
  set(gca,'YTick',[0 0.5],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  if ((k==1) | (k==5) | (k==9) | (k==13)); 
     set(gca,'YTickLabel',[0 0.5],'FontSize',14);
  end;
  axis([-5 5 0 0.8]);
  set(gca,'FontSize',14);

%  subplot(4,4,k+1);
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);
%  if (k>12 | k == 11 | k == 12); xlabel('Res. (std. dev.)','FontSize',14);end;
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  if (k>10); 
     xlabel('Res. (std. dev.)','FontSize',14);
     set(gca,'XTickLabel',[-4 0 4],'FontSize',14);
  end;
  set(gca,'YTick',[0 0.5],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  hold on;box on;
  [n1,xout] = hist(reshat,x);
  area = sum(n1) * (xout(2)-xout(1));
  n1 = n1/area;
%  n1 = n1/(N2*(xout(2)-xout(1)));
  stairs(xout,n1,'k');
  plot(xx,nd,'--k')
  text(3.5,0.6,ttl(ifreq));
  axis([-5 5 0 0.8]);
  set(gca,'FontSize',14);

  z = runtest(res(:,ifreq));
  p = normcdf(z,0,1);
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

  save bla;
  [pass,value,bereik,opp,dplot] = kolmogorov(res(:,ifreq)./sd2(:,ifreq),k);
  if(value < foo(end,2))
      idx = find(value<foo(:,2));
      idx = idx(1);
  else
      idx = length(foo);
  end
  p1 = (100-foo(idx,1))/100;
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

  z = runtest(reshat(:));
  p = normcdf(z,0,1);
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

  [pass2,value2,bereik2,opp2,dplot2] = kolmogorov(reshat(:),k+1);
  if(value2 < foo(end,2))
      idx = find(value2<foo(:,2));
      idx = idx(1);
  else
      idx = length(foo);
  end
  p2 = (100-foo(idx,1))/100;
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
% Plot uncorrected!
%
  gc3 = figure(4);
  hold on;
%  subplot(4,4,k);
  subplot('Position',[loc(1,k) loc(2,k) spw sph]);
  hold on; box on;
  plot(bereik,opp,'--k');plot(bereik,dplot,'k');
%  if (k > 12 | k==11 | k == 12); xlabel('Res. (std. dev.)','FontSize',14);end;
%  if(k==1 | k==5 |k==9|k==13);ylabel('cum. prob.','FontSize',14);end;
  if(k==5);ylabel('cum. prob.','FontSize',14);end;
  if (k > 10); xlabel('Res. (std. dev.)','FontSize',14);end;
%  ylabel('cum. prob.','FontSize',14);
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
%  if(k > 12 | k==11 | k == 12);set(gca,'XTickLabel',{'-4';'-2';'0';'2';'4'});end;
  if (k==1 | k==5 |k==9|k==13);set(gca,'YTickLabel',{'0';'';'1'});end;
  if(k > 10);set(gca,'XTickLabel',{'-2';'0';'2'});end;
%  set(gca,'YTickLabel',{'0';'';'1'});
  text(-3,1.,ttl(ifreq),'FontSize',14)
  set(gca,'FontSize',14);

%
% Plot corrected!
%
%  subplot(4,4,k+1);
  subplot('Position',[loc(1,k+1) loc(2,k+1) spw sph]);
  hold on; box on;
  plot(bereik2,opp2,'--k');plot(bereik2,dplot2,'k');
% if (k > 12 | k == 11 | k == 12); xlabel('Res. (std. dev.)','FontSize',14);end;
  if (k > 10); xlabel('Res. (std. dev.)','FontSize',14);end;
  set(gca,'XLim',[-4 4],'YLim',[0 1.2]);
  set(gca,'XTick',[-2 0 2],'XTickLabel',{'';'';'';'';''});
  set(gca,'YTick',[0 0.5 1],'YTickLabel',{'';'';''});
%  if(k > 12 | k==11 | k == 12);set(gca,'XTickLabel',{'-4';'-2';'0';'2';'4'});end;
  if(k > 10);set(gca,'XTickLabel',{'-2';'0';'2'});end;
  text(-3,1.,ttl(ifreq),'FontSize',14)
  set(gca,'FontSize',14);




  k = k+2;

end;

saveas(gc1,'plot1.eps','epsc2');
saveas(gc2,'plot2.eps','epsc2');
saveas(gc3,'plot3.eps','epsc2');

save(cov_mat_file2,'F');

covmat = inv(F(1).csave);
save(cov_mat_file,'covmat','-ascii');
for ifreq = 2:nfreq

    covmat = inv(F(ifreq).csave);
    save(cov_mat_file,'covmat','-ascii','-append');

end;

return;

