%
% Compute and plot correlation matrix
%

function [] = plot_cor();

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',8,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

sample    = 'tmp_sample.mat';
plotfile1 = 'x_s01_1_3_20_3lay_corr_mat.eps';
corrfile  = 'x_s01_1_3_20_3lay_corr_mat.mat';
covfile   = 'x_s01_1_3_20_3lay_cov_mat.mat';

load(sample);

m = A(:,5:8);

minlim = [ 1   -1  -1  1];
maxlim = [50  -30 -60 20];

dlim = maxlim-minlim;
npar = length(m(1,:));
ndat = length(m);

for i = 1:npar
  m(:,i) = (m(:,i) - minlim(i))/dlim(i);
end 

mcov = zeros(npar,npar);
mcor = zeros(size(mcov));
mcorp = zeros(size(mcov));

mcov = cov(m);

k = 1;
for i = 1:npar
  for j = 1:npar
    mcor(i,j) = mcov(i,j)/sqrt(mcov(i,i)*mcov(j,j));
  end
end
for i = 1:npar
  for j = i:npar
    if(abs(mcor(i,j)) > 0.55 & abs(mcor(i,j)) < 1.0)
      a(k,:) = [i j];
      mcorplot(k)=mcor(i,j);
      k = k + 1;
    end
  end
end
mcorp = flipud(mcor);

disp('mcor');
disp(mcor);
disp('mcov');
disp(mcov);
save('cov0.txt','mcov','-ascii');
%ylabels = fliplr({'h1','c1','r1t','r1b','a1','h2','c2','r2','a2',...
%                  'h3','c3','r3','a3','h4','c4','r4','a4',...
%                       'c5','r5','a5'});
%                  'h5','c5','r5','a5','h6','c6','r6','a6',...
%xlabels = {'h1','c1','r1t','r1b','a1','h2','c2','r2','a2',...
%                  'h3','c3','r3','a3','h4','c4','r4','a4',...
%                       'c5','r5','a5'};
%                  'h5','c5','r5','a5','h6','c6','r6','a6',...

ylabels = fliplr({'h1','c1','r1','a1','h2','c2','r2','a2',...
                  'h3','c3','r3','a3','h4','c4','r4','a4','h5','c5','r5','a5','c6','r6','a6'});
xlabels = {'h1','c1','r1','a1','h2','c2','r2','a2',...
                  'h3','c3','r3','a3','h4','c4','r4','a4','h5','c5','r5','a5','c6','r6','a6'};
%ylabels = fliplr({'h1','c1','h2','c2','h3','c3','h4','c4',...
%                  'h5','c5','h6','c6','a1','a2','a3','a4','a5','a6','a7'});
%xlabels = {'h1','c1','h2','c2','h3','c3','h4','c4',...
%           'h5','c5','h6','c6','a1','a2','a3','a4','a5','a6','a7'};
gc1=figure(1);
hold on;box on;grid on;
for i = 1:npar
%  plot(mcorp(i,:));
  offset = ((i-1)*2)+1;
  area([(npar-i)+1:npar],mcorp(i,(npar-i)+1:end)+offset,offset,'FaceColor',[0 0 0]);
  set(gca,'YLim',[0 2*npar+1],'YTick',[1:2:npar*2]);
  set(gca,'XLim',[1 npar],'XTick',[1:1:npar]);
  set(gca,'YTickLabel',ylabels);
  set(gca,'XTickLabel',xlabels);
%  set(gca,'YTickLabel',{'\alpha','c_b','c_t','\nu','\rho_b','\rho_t','h'});
%  set(gca,'XTickLabel',{'h','\rho_t','\rho_b','\nu','c_t','c_b',...
%          '\alpha'});

end
%saveas(gca,plotfile1,'epsc2');
exportfig(gc1,plotfile1,opts);

save(corrfile,'mcor');
save(covfile,'mcov');
a'
mcorplot
svd(mcor)
svd(mcov)

return;

