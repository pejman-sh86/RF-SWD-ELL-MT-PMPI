function plot_models;
%
% Loads: map.mat (from plot_hist)
%

%set(0, 'DefaultFigurePaperPosition', [0 0 7 4]); 
dex = 0;

i_rc = 0;
paper = 1;
errbar = 1;	%  0 == print cores only
		%  1 == print errbars on every 10th datum
		%  2 == print averaged cores and bounds
barstep = 10;
nlayer = 20;
laythick = 0.1;
nmod = 10000;
nrbin = 100;

sample = 'real2_sample.dat';
mapfile = strrep(sample,'sample.dat','map.mat');
plotfile2 = strrep(sample,'sample.dat','core_inv.eps');

load(mapfile);

%
% Use ML model
%
xc5 = [1.03 1.3091698  1.4918757  0.49265393  1468.6338  1469.0314]';
xc6 = [1.12 1.1989922  1.4990106  0.21987834  1471.061   1460.0173]';

sc_5 = 1.003551;
sc_6 = 1.010277;
sr_5 = 1.006025;
sr_6 = 0.996599;
C = load('core5.dat');
D = load('core6.dat');
F(1).depth = C(:,1);
F(1).depth_mx = F(1).depth(end);
F(2).depth = D(:,1);
F(2).depth_mx = F(2).depth(end);
if(i_rc == 0)
  F(1).dat = C(:,3)/sr_5;
  F(2).dat = D(:,3)/sr_6;
  F(1).datc = C(:,2)/sc_5;
  F(2).datc = D(:,2)/sc_6;
end

F(1).lay_thick = 1/length(F(1).depth);
F(1).znorm=F(1).lay_thick/2:F(1).lay_thick:1-F(1).lay_thick/2;
F(2).lay_thick = 1/length(F(2).depth);
F(2).znorm=F(2).lay_thick/2:F(2).lay_thick:1-F(2).lay_thick/2;
size(F(1).znorm)

znorm = zeros(10,1);
lay_thick = 0.01;
znorm=lay_thick/2:lay_thick:1-lay_thick/2;

c = zeros(1,length(znorm)+1);
rho = zeros(1,length(znorm)+1);

%
% ML model:
%
h_xc5 = xc5(1);
rhot_xc5 = xc5(2);
rhob_xc5 = xc5(3);
nu_xc5 = xc5(4);
ct_xc5 = xc5(5);
cb_xc5 = xc5(6);

h_xc6 = xc6(1);
rhot_xc6 = xc6(2);
rhob_xc6 = xc6(3);
nu_xc6 = xc6(4);
ct_xc6 = xc6(5);
cb_xc6 = xc6(6);
%
% MAP model:
%
h_map = xmap(1);
rhot_map = xmap(2);
rhob_map = xmap(3);
nu_map = xmap(4);
ct_map = xmap(5);
cb_map = xmap(6);
alpha_map = xmap(7);

z_xc5 = [0 F(1).znorm] *h_xc5;
rho_xc5 = rhot_xc5 + (sin([0 F(1).znorm] .* pi/2)).^nu_xc5 .* (rhob_xc5 - rhot_xc5);
c_xc5 = ct_xc5 + (cb_xc5 - ct_xc5) .* [0 F(1).znorm];
z_xc6 = [0 F(2).znorm] *h_xc6;
rho_xc6 = rhot_xc6 + (sin([0 F(2).znorm] .* pi/2)).^nu_xc6 .* (rhob_xc6 - rhot_xc6);
c_xc6 = ct_xc6 + (cb_xc6 - ct_xc6) .* [0 F(2).znorm];

z_map = [0 znorm] *h_map;
rho_map = rhot_map + (sin([0 znorm] .* pi/2)).^nu_map .* (rhob_map - rhot_map);
c_map = ct_map + (cb_map - ct_map) .* [0 znorm];

%saveas(gca,plotfile1,'epsc2');

  figure(2);
  hold on;
  subplot('Position',[.12 .15 .4 .8]);
  box on;
  hold on;
  plot('v6',rho_xc5,z_xc5,'-g','MarkerSize',8);
  plot('v6',rho_xc6,z_xc6,'-r','MarkerSize',8);
  plot('v6',rho_map,z_map,'-k.','MarkerSize',8);
  if(paper == 1)
    set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
        'YTickLabel',[0 0.5 1.0 1.5],'FontSize',16);
  else
    set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
        'YTickLabel',[0 0.5 1.0 1.5],'FontSize',16);
  end
  set(gca,'XLim',[1.25 1.615]);
  set(gca,'XTickLabel',[1.3 1.4 1.5 1.6],'XTick',[1.3 1.4 1.5 1.6]);
  xlabel('density (g/cm^3)');
  ylabel('depth (m)');
%return;
  subplot('Position',[.56 .15 .4 .8]);
  box on;
  hold on;
  plot('v6',c_xc5,z_xc5,'-g','MarkerSize',8);
  plot('v6',c_xc6,z_xc6,'-r','MarkerSize',8);
  plot('v6',c_map,z_map,'-k','MarkerSize',8);
  set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
      'YTickLabel',[],'FontSize',16);
  set(gca,'XLim',[1445 1512]);
  set(gca,'XTickLabel',[1450 1470 1490 1510],'XTick',[1450 1470 1490 1510]);
  xlabel('velocity (m/s)');
%  ylabel('depth [m]');
  saveas(gca,plotfile2,'epsc2');
  hgsave('fig17.fig','-v6');
%  saveas(gca,plotfile2,'png');

%save bla rho_xc5 z_xc5 rho_xc6 z_xc6 rho_map z_map c_xc5 c_xc6 c_map;
%
% Calc core misfit
%
mis1 = 1/length(F(1).dat)*sum(abs(F(1).dat - rho_xc5(2:end)')./F(1).dat)
mis2 = 1/length(F(2).dat)*sum(abs(F(2).dat - rho_xc6(2:end)')./F(2).dat)

mis3 = 1/length(F(1).datc)*sum(abs(F(1).datc - c_xc5(2:end)')./F(1).datc)
mis4 = 1/length(F(2).datc)*sum(abs(F(2).datc - c_xc6(2:end)')./F(2).datc)

abs(F(1).datc - c_xc5(2:end)')

return;

