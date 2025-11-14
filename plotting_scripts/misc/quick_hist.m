function plot_hist(nrbin);
set(0, 'DefaultFigurePaperPosition', [0 0 8.5 7]);

sample = 'sample_tmp.mat';
load(sample);

npar = size(A,2)-2;

dat1 = A(:,2:npar+1);
clear A;

xmin = min(dat1)
xmax = max(dat1)
xdiff = abs(xmax - xmin);
delta = xdiff/nrbin;
N = length(dat1);
for i =0:nrbin
  edges(i+1,:) = xmin + xdiff/nrbin * i;
end

nx = 4
ny = 3
xim = 0.03;
yim = 0.1;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);
nsubpfig = nx*ny;
nfig = ceil(npar/nsubpfig);

kdat = 1;
for ifig = 1:nfig

  figure(ifig);hold on;
  for k=1:nsubpfig

    subplot('Position',[loc(1,k) loc(2,k) spw sph]);
    hold on;box on
    n1 = histc(dat1(:,kdat),edges(:,kdat));
    n1 = n1/(N * delta(kdat))*xdiff(kdat);
    stairs(edges(:,kdat),n1,'k');
    hold on;
    stdev = std(dat1(:,kdat));
    exp_val = mean(dat1(:,kdat))
  
    kdat = kdat+1;
    if kdat>npar;break;end; 
  end
end
return;
