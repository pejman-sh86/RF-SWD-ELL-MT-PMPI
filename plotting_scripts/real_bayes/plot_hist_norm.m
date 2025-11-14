function [stat_fgs]=plot_hist9(nrbin);
             %plot_hist9(msyn,60);
%
% Plot histograms for noise FGS.
%             
% msyn = [1.0 1.3 1.5 0.5 1465.0 1475.0 0.07]';

A = load('test3_sample1.dat');
B = load('test3_sample2.dat');
%m = [A(:,2:end-1) ; B(:,2:end-1)];
dat1(:,:) = [A(:,2:8) ; B(:,2:8)];

%A = load('../../real_bayesian/tests/test1_sample1.dat');
%B = load('../../real_bayesian/tests/test1_sample2.dat');
%dat2(:,:) = [A(:,2:8) ; B(:,2:8)];

%A = load('../../mlestimate/tests/test4_sample1.dat');
%B = load('../../mlestimate/tests/test4_sample2.dat');
%dat3(:,:) = [A(:,2:8) ; B(:,2:8)];

col = ['r' 'b' 'g'];

xmin = [0.8 1.3 1.3 0.0 1465 1460 0.15];
xmax = [2.5 1.4 1.6 1.5 1475 1480 0.35];
xdiff = xmax - xmin;
for i =1:nrbin
  edges(i,:) = xmin + xdiff/nrbin * i;
end
knum=[1 2 5 4 3 6 7]; %plot order

%stat_fgs=zeros(size(A,2)-2,4);

ttl=[{'h [m]'}, {'\rho_t [g/cm^3]'},{'\rho_b [g/cm^3]'},...
     {'\nu'}, {'c_t [m/s]'}, {'c_b [m/s]'},{'alpha'},...
     {'sd 315 Hz'},{'sd 400 Hz'}];
%     {'\rho_3 [g/cm^3]'},{'c_3 [m/s]'}];
    

figure
for k=1:size(dat1,2)
  subplot(4,2,k)
  p=knum(k);
  n1 = histc(dat1(:,p),edges(:,p));
  n1 = n1/max(n1);
%  n2 = histc(dat2(:,p),edges(:,p));
%  n2 = n2/max(n2);
%  n3 = histc(dat3(:,p),edges(:,p));
%  n3 = n3/max(n3);
  stairs(edges(:,p),n1,col(1));
  hold on;
%  stairs(edges(:,p),n2,col(2));
%  hold on;
%  stairs(edges(:,p),n3,col(3));
%    hist(m(:,p),nrbin);hold on;
%    [n,x]=hist(m(:,p),nrbin); [ii jj]=  max(n);mx(k)=ii;
%    stat_fgs(p,:)=[x(jj) median(m(:,p)) mean(m(:,p)) 2*std(m(:,p))];
end

for k=1:size(dat1,2)
  subplot(4,2,k);
  p=knum(k);
  set(gca,'Xlim',[xmin(p) xmax(p)]);title(ttl(p));
end

