
set(0, 'DefaultFigurePaperPosition', [0 0 4 8]);
load syn_cov_mat_est2;
%G = F;
%load blde4_syn.mat

width = 0.4
height = 0.8
llw = [0.1 0.54];
llh = [0.1];

gca1 = figure(1);
%subplot('Position',[llw(1) llh(1) width height])
pcolor(F(8).csave);
shading interp;
%colormap('gray');
set(gca,'YDir','reverse','FontSize',14);
xlabel('Data points','FontSize',14);
ylabel('Data points','FontSize',14);

return;
subplot('Position',[llw(2) llh(1) width height])
pcolor(G(8).csave);
shading interp;
set(gca,'YTick',[]);
%colormap('gray');
set(gca,'YDir','reverse','FontSize',14);
xlabel('Data points','FontSize',14);

%saveas(gca1,'CAA_plot1.tif','tiffn')
return;
