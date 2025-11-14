%
% Compute and plot correlation matrix
%

function corr_mat;

set(0, 'DefaultFigurePaperPosition', [0 0 3 3]);

sample1 = 'odd_sample1.dat';
sample2 = 'odd_sample2.dat';
plotfile1 = strrep(sample1,'sample1.dat','corr_mat.eps');
A = load(sample1);
B = load(sample2);

m = [A(:,2:17) ; B(:,2:17)];
size(m)

sdlow = [1 1 1 1 1 1 1 1 1];
sdup =  [4 4 4 4 4 4 4 4 4];
minlim = [00.10  1.0  1.3  0.0  1450  1450  0.005 sdlow]';
maxlim = [2.50  2.2  2.2  1.50  1550  1550  0.5 sdup]';

npar = 7;
ndat = length(m);

mm = zeros(size(m));
for i = 1:npar
  m(:,i) = (m(:,i) - minlim(i))/maxlim(i);
  mm(:,i) = m(:,i) - sum(m(:,i))/ndat;
end 

mcov = zeros(npar,npar);
mcor = zeros(size(mcov));
mcorp = zeros(size(mcov));
for i = 1:npar
  for j = 1:npar
    mcov(i,j) = sum(mm(:,i).*mm(:,j));
  end
end
mcov = mcov/ndat;

for i = 1:npar
  for j = 1:npar
    mcor(i,j) = mcov(i,j)/sqrt(mcov(i,i)*mcov(j,j));
  end
end

%k = 1.;
%for i = 1:npar
%  mcorp(i,:) = mcor(i,:) + k;
%  k = k + 2.;
%end

mcorp = mcor+repmat(((2*npar)-1:-2:0)',1,npar);

size(mcorp)
figure(1);
for i = 1:npar
  grid on
  plot(mcorp(i,:));
  set(gca,'YTick',[1,3,5,7,9,11,13]);
  set(gca,'XTick',[1,2,3,4,5,6,7]);
  set(gca,'YTickLabel',{'\alpha','c_b','c_t','\nu','\rho_b','\rho_t','h'});
  set(gca,'XTickLabel',{'h','\rho_t','\rho_b','\nu','c_t','c_b',...
          '\alpha'});
  hold on;
end
saveas(gca,plotfile1,'epsc2');

mcor
return

