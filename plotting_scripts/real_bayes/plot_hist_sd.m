function [stat_fgs]=plot_hist9(nrbin);
             %plot_hist9(msyn,60);
%
% Plot histograms for noise FGS.
%             
% msyn = [1.0 1.3 1.5 0.5 1465.0 1475.0 0.07]';

A = load('test1_sample1.dat');
B = load('test1_sample2.dat');
%m = [A(:,2:end-1) ; B(:,2:end-1)];
dat1(:,:) = [A(:,9:17) ; B(:,9:17)];

A = load('../../real_bayesian/tests/test1_sample1.dat');
B = load('../../real_bayesian/tests/test1_sample2.dat');
dat2(:,:) = [A(:,9:17) ; B(:,9:17)];

ml = [2.96 2.72 2.52 2.56 2.31 2.04 1.36 1.49 1.42];

col = ['r' 'b'];

xmin = [1 1 1 1 1 1 1 1 1];
xmax = [4 4 4 4 4 4 4 4 4];
xdiff = xmax - xmin;
for i =1:nrbin
  edges(i,:) = xmin + xdiff/nrbin * i;
end
knum=[1 2 3 4 5 6 7 8 9]; %plot order

%stat_fgs=zeros(size(A,2)-2,4);

ttl=[{'sd 315 Hz'}, {'sd 400 Hz'},{'sd 500 Hz'},...
     {'sd 630 Hz'}, {'sd 800 Hz'}, {'sd 1000 Hz'},{'sd 1250 Hz'},...
     {'sd 1600 Hz'},{'sd 2000 Hz'}];
%     {'\rho_3 [g/cm^3]'},{'c_3 [m/s]'}];
    

figure
for k=1:size(dat1,2)
  subplot(3,3,k)
  p=knum(k);
  n1 = histc(dat1(:,p),edges(:,p));
  n1 = n1/max(n1);
  n2 = histc(dat2(:,p),edges(:,p));
  n2 = n2/max(n2);
  stairs(edges(:,p),n1,col(1));
  hold on;
  stairs(edges(:,p),n2,col(2));
  hold on;
  x = [ml(k) ml(k)];
  y = [0 1];
  plot(x,y,'--y','LineWidth',2);
%    hist(m(:,p),nrbin);hold on;
%    [n,x]=hist(m(:,p),nrbin); [ii jj]=  max(n);mx(k)=ii;
%    stat_fgs(p,:)=[x(jj) median(m(:,p)) mean(m(:,p)) 2*std(m(:,p))];
end

for k=1:size(dat1,2)
  subplot(3,3,k);
  p=knum(k);
  set(gca,'Xlim',[xmin(p) xmax(p)]);title(ttl(p));
end

