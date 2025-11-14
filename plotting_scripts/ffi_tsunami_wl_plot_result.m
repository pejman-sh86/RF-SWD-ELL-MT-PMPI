load tohoku_data_sample.mat;

figure();

pmean=load('mean.txt');
figure();
imagesc(pmean);
colormap(darkb2r(-1,6.5));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
axis equal;
set(gca,'XLim',[0 9],'YLim',[0 17]);

tru=load('sl1_tot_sim.dat');
figure();
imagesc(tru);
colormap(darkb2r(-1,6.5));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
axis equal;
set(gca,'XLim',[0 9],'YLim',[0 17]);

figure();
imagesc(tru-pmean);
colormap(darkb2r(-2,2));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
axis equal;
set(gca,'XLim',[0 9],'YLim',[0 17]);

tru=load('var.txt');
figure();
imagesc(tru);
colormap(darkb2r(-2,2));
cb = colorbar('peer',gca,'FontSize',14,'location','NorthOutside');
axis equal;
set(gca,'XLim',[0 9],'YLim',[0 17]);


