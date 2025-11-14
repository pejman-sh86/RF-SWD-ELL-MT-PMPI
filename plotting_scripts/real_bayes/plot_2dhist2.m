%
% Plot 2-D Histograms from FGS samples
%
function plot_2dhist;

set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);

sample1 = 'third_second_sample1.dat';
sample2 = 'third_second_sample2.dat';
plotfile1 = strrep(sample1,'sample1.dat','joints.eps');
A = load(sample1);
B = load(sample2);


%
% h = 1; r_t = 2; r_b = 3; nu = 4; c_t = 5; c_b = 6; alpha = 7
%
% Test 6,8
%ip1 = [1 2 2 4];
%ip2 = [3 4 7 7];
ip1 = [1 2 4 6];
ip2 = [3 4 7 7];

% Test 3,5
%ip1 = [1 1 1 2 2 4];
%ip2 = [2 4 6 4 6 6];
%xlim1 = [0.98 0.98 0.98 1.298 1.298 0.482];
%xlim2 = [1.02 1.02 1.02 1.302 1.302 0.52];
%ylim1 = [1.298 0.482 1474.8 0.482 1474.8 1474.8];
%ylim2 = [1.302 0.52 1475.2 0.52 1475.2 1475.2];
nplot = length(ip1);
nrow = ceil(nplot/2)

m = [A(:,2:8) ; B(:,2:8)];
msyn = [1.7936    1.3626    1.4684    0.8392 1471.8052 1466.0385    0.2505];

npar = 7;
ndat = length(m);

figure(1);
for i = 1:nplot
  [n,x,nbins] = histmulti5([m(:,ip2(i)) m(:,ip1(i))],[50 50]);

%set(gca,'ydir','reverse');
  subplot(nrow,2,i);
  pcolor(x(:,2),x(:,1),n),colorbar,shading interp;
%  pcolor(x(:,2),x(:,1),n),colorbar,shading flat;

  if ip1(i) == 1
    xlabel('h [m]');
  elseif ip1(i) == 2
    xlabel('\rho_t [g/ccm]');
  elseif ip1(i) == 3
    xlabel('\rho_b [g/ccm]');
  elseif ip1(i) == 4
    xlabel('\nu');
  elseif ip1(i) == 5
    xlabel('c_t [m/s]');
  elseif ip1(i) == 6
    xlabel('c_b [m/s]');
  elseif ip1(i) == 7
    xlabel('\alpha [dB/\lambda]');
  end
  
  if ip2(i) == 1
    ylabel('h [m]');
  elseif ip2(i) == 2
    ylabel('\rho_t [g/ccm]');
  elseif ip2(i) == 3
    ylabel('\rho_b [g/ccm]');
  elseif ip2(i) == 4
    ylabel('\nu');
  elseif ip2(i) == 5
    ylabel('c_t [m/s]');
  elseif ip2(i) == 6
    ylabel('c_b [m/s]');
  elseif ip2(i) == 7
    ylabel('\alpha [dB/\lambda]');
  end
  hold on;
  plot(msyn(ip1(i)),msyn(ip2(i)),'w+','MarkerSize',10);
  hold off;

%  set(gca,'XLim',[xlim1(i) xlim2(i)]);
%  set(gca,'YLim',[ylim1(i) ylim2(i)]);
end

saveas(gca,plotfile1,'epsc2');


return
