%
% Plot 2-D Histograms from FGS samples
%
function plot_2dhist;

set(gcf, 'renderer', 'painters')

sample    = 'noisy_data_sample.mat';
plotfile1 = 'x_jan_1_sph_ci2_joints_rot.eps';

load(sample);
titlearray = {'a' 'b' 'c' 'd' 'e' 'f'};
%
% h = 1; r_t = 2; r_b = 3; nu = 4; c_t = 5; c_b = 6; alpha = 7
%
ip1 = [1  1  1  2  2  3];
ip2 = [2  3  4  3  4  4];

minlim = [ 0, -30, -60, 0.0];
maxlim = [50,  -1,  -1, 0.1];

nplot = length(ip1);
nx = 3;
ny = 2;
xim = 0.28/nx;
yim = 0.28/ny;
xymarg = [0.14 0.04 0.04 0.14];
opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',14,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

m = A(:,[5:7,13]);

ndat = length(m);

xlabels = {'h','dvp','dvs','drho'};

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc1=figure(1);
for i = 1:nplot
  %[n,x,nbins] = histmulti5([m(:,ip2(i)) m(:,ip1(i))],[100 100]);

  subplot('Position',[loc(1,i) loc(2,i) spw sph]);
  set(gcf, 'renderer', 'painters');
  %colormap(1-gray);
%  cloudPlot(m(:,ip1(i)),m(:,ip2(i)),[minlim(ip1(i)) maxlim(ip1(i)) minlim(ip2(i)) maxlim(ip2(i))],false,[200 200]);
  cloudPlot(m(:,ip1(i)),m(:,ip2(i)),[],false,[2000 2000]);
  set(gca,'FontSize',14,'layer','top');
  xlabel(xlabels(ip1(i)),'FontSize',14);
  ylabel(xlabels(ip2(i)),'FontSize',14);
  hold on;box on;
end

%saveas(gca,plotfile1,'epsc2');
%exportfig(gc1,plotfile1,opts);


return
