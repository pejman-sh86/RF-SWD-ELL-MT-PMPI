function plot_hist9(nrbin);
             %plot_hist9(msyn,60);
%
% Plot histograms for noise FGS.
% 1) Saves MAP into map.mat
%             
set(0, 'DefaultFigurePaperPosition', [0 0 8 7]);
npar = 7;

sample1 = 'sample1.dat';
sample2 = 'sample2.dat';
plotfile1 = strrep(sample1,'sample1.dat','marginals.eps');
A = load(sample1);
B = load(sample2);

[cost_min i_min] = min([A(:,1) ; B(:,1)])

dat1 = [A(:,2:npar+1) ; B(:,2:npar+1)];
clear A;
clear B;

xmap = dat1(i_min,:);
xtrue = [2.0930214      1.3508944      1.4841878     0.73839097 ... 
         1472.8123 1467.5688     0.25152184]';
save map.mat xmap;
xmin = [min(dat1(:,1)) min(dat1(:,2)) min(dat1(:,3)) min(dat1(:,4)) ...
       min(dat1(:,5)) min(dat1(:,6)) min(dat1(:,7))];
xmax = [max(dat1(:,1)) max(dat1(:,2)) max(dat1(:,3)) max(dat1(:,4)) ...
        max(dat1(:,5)) max(dat1(:,6)) max(dat1(:,7))];
%xmin = [min(dat1(:,1)) min(dat1(:,2)) min(dat1(:,3)) min(dat1(:,4)) ...
%        min(dat1(:,5)) min(dat1(:,6)) min(dat1(:,7)) min(dat1(:,8)) ... 
%        min(dat1(:,9)) ];
%xmax = [max(dat1(:,1)) max(dat1(:,2)) max(dat1(:,3)) max(dat1(:,4)) ...
%        max(dat1(:,5)) max(dat1(:,6)) max(dat1(:,7)) max(dat1(:,8)) ... 
%        max(dat1(:,9)) ];
%xmin = [min(dat1(:,1)) min(dat1(:,2)) min(dat1(:,3)) min(dat1(:,4)) ...
%        min(dat1(:,5)) min(dat1(:,6)) min(dat1(:,7)) min(dat1(:,8)) ...
%        min(dat1(:,9)) ];
%xmax = [min(dat1(:,1))+0.05 min(dat1(:,2))+0.1 min(dat1(:,3))+0.1 ...
%        min(dat1(:,4))+0.6 min(dat1(:,5))+10 min(dat1(:,6))+10 ...
%        min(dat1(:,7))+0.1 min(dat1(:,8))+10 min(dat1(:,9))+0.1];
%xmin = [1.4 1.2 1.4 0.3 1460 1440 0.1];
%xmax = [3.5 1.5 1.7 1.3 1485 1465 0.5];
xdiff = abs(xmax - xmin);
delta = xdiff/nrbin;
N = length(dat1);
for i =0:nrbin
  edges(i+1,:) = xmin + xdiff/nrbin * i;
end
col = ['r' 'b' 'g'];

ttl=[{'h [m]'}, {'\rho_t [g/ccm]'},{'\rho_b [g/ccm]'},...
     {'\nu'}, {'c_t [m/s]'}, {'c_b [m/s]'},{'alpha'},...
     {'c_3 [m/s]'},{'\rho_3 [g/ccm]'}];
%     {'sd 315 Hz'},{'sd 400 Hz'}];

rank = [1 2 5 4 3 6 7];
    
figure(1);
llw = [0.03 0.36 0.69 0.03 0.36 0.69 0.03 0.36 0.69];
llh = [0.73 0.4 0.07];
width = 0.3;
height = 0.26
ii = 1;
for k=1:npar
  if k <= 3
    subplot('Position',[llw(ii) llh(1) width height]);
    ii = ii + 1;
  elseif k <= 6
    subplot('Position',[llw(ii) llh(2) width height]);
    ii = ii + 1;
  else
    subplot('Position',[llw(ii) llh(3) width height]);
    ii = ii + 1;
  end
  hold on;box on
  kk = rank(k);
  n1 = histc(dat1(:,kk),edges(:,kk));
  n1 = n1/(N * delta(kk))*xdiff(kk);
  stairs(edges(:,kk),n1,col(2));
  hold on;
  stdev = std(dat1(:,kk));
  exp_val = mean(dat1(:,kk))
  gamma = abs(exp_val - xmap(kk));
  
%  t = xmin(kk):delta(kk):xmax(kk);
%  realgauss = 1/(sqrt(2*pi)*stdev)*exp(-(t-exp_val).^2/(2*stdev*stdev));
%  realgauss = realgauss*xdiff(kk);
%  plot(t,realgauss,'--k');

%  hist(m(:,p),nrbin);hold on;
%  [n,x]=hist(m(:,p),nrbin); [ii jj]=  max(n);mx(k)=ii;
%  stat_fgs(p,:)=[x(jj) median(m(:,p)) mean(m(:,p)) 2*std(m(:,p))];
end

%
% 95% HPD's:
%
% Use finer bins for the HPDs!
%
nrbin = 29*nrbin;
for i =0:nrbin 
  edges(i+1,:) = xmin + xdiff/nrbin * i;
end
%outfile = 'hpds.dat';
%fid = fopen(outfile,'w');
ii = 1;
for k = 1:npar
  kk = rank(k);
  n1 = histc(dat1(:,kk),edges(:,kk));
  nintyfive = sum(n1)-5*sum(n1)/100;
  lb = 0;
  rb = 0;
  
  p = 1;
  for i = 1:length(n1)
    s = 0;
    stopj = 0;
    for j = i:length(n1)
      if(stopj == 0)
        s = s + n1(j);
        if(s >= nintyfive)
          lb(p) = edges(i,kk);
          rb(p) = edges(j,kk);
          p = p+1;
          stopj = 1;
        end
      end
    end
  end
  [nf_int,nf_index] = min(abs(lb-rb));
%  fprintf(fid,'%10.4f',lb(nf_index),rb(nf_index));fprintf(fid,'\n');

%
% Plot True and 95 HPD
% 
  if k <= 3
    subplot('Position',[llw(ii) llh(1) width height]);
    ii = ii + 1;
  elseif k <= 6
    subplot('Position',[llw(ii) llh(2) width height]);
    ii = ii + 1;
  else
    subplot('Position',[llw(ii) llh(3) width height]);
    ii = ii + 1;
  end
  hold on;box on
  set(gca,'XLim',[xmin(kk) xmax(kk)],'YLim',[0 10],'YTick',[],'FontSize',12);
%  set(gca,'XLim',[xmin(k) xmax(k)]);
  xlabel(ttl(kk),'FontSize',14);
  ylims = get(gca,'Ylim');
  exp_val = mean(dat1(:,kk));
  min_val = min(dat1(:,kk));
  max_val = max(dat1(:,kk));
%  plot([exp_val exp_val],[0 ylims(2)],'--k');
  plot([xmap(kk) xmap(kk)],[0 ylims(2)],'--k');
  plot([xtrue(kk) xtrue(kk)],[0 ylims(2)],'--r');
  plot([lb(nf_index) lb(nf_index)],[0 ylims(2)],':k','LineWidth',2);
  plot([rb(nf_index) rb(nf_index)],[0 ylims(2)],':k','LineWidth',2);
end
%fclose(fid);
saveas(gca,plotfile1,'epsc2');

return;
