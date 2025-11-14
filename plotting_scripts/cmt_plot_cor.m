function [] = plot_cor(filename,NCMT,NLOC);
%
% Compute and plot correlation matrix
%

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
% Output files
corrfile        = strcat(filebase,'_corr.dat');
covfile        = strcat(filebase,'_cov.dat');

plotfilecor    = strcat(filebase,'_cormat.');

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

[IMAP,ICOV,NVMX,NPV,ILOC,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay]=cmt_read_parfile(parfile);
[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

load(filename);

m = A(:,5:4+NVMX*NPV);

minlim(1:NVMX,1:5) = -6.;
maxlim(1:NVMX,1:5) =  6.;
minlim(1:NVMX,NCMT+1:NPV) = [minlat',minlon',mindepth',mindelay'];
maxlim(1:NVMX,NCMT+1:NPV) = [maxlat',maxlon',maxdepth',maxdelay'];

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

xlabels = ({'M_1','M_2','M_3','M_4','M_5',...
            'lat','lon','depth','t_d'});
ylabels = fliplr({'M_1','M_2','M_3','M_4','M_5',...
            'lat','lon','depth','t_d'});

xlabels2 = xlabels;
ylabels2 = ylabels;
for ivo = 2:NVMX;
  xlabels2 = [xlabels2,xlabels];
  ylabels2 = [ylabels2,ylabels];
end;
figcormat=figure();
hold on;box on;grid on;
for i = 1:npar
  offset = ((i-1)*2)+1;
  area([(npar-i)+1:npar],mcorp(i,(npar-i)+1:end)+offset,offset,'FaceColor',[0 0 0]);
  set(gca,'YLim',[0 2*npar+1],'YTick',[1:2:npar*2]);
  set(gca,'XLim',[1 npar],'XTick',[1:1:npar]);
  set(gca,'YTickLabel',ylabels2);
  set(gca,'XTickLabel',xlabels2);
end
print(figcormat,'-painters','-r250',strcat(plotfilecor,plotext2),'-dpng');

save(corrfile,'mcor');
save(covfile,'mcov');
a'
mcorplot
svd(mcor)
svd(mcov)

return;
