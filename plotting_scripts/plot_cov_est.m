
dex = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 8 8]);

sample = 'real2_sample.dat';

if dex == 0
  cov_mat_file = strrep(sample,'sample.dat','cov_mat_est.mat');
  plotfile = strrep(sample,'sample.dat','cov_mat_est.eps');
else
  cov_mat_file = strrep(sample,'sample.dat','cov_mat_est_dexp.mat');
  plotfile = strrep(sample,'sample.dat','cov_mat_est_dexp.eps');
end

load(cov_mat_file);
nfreq = 8;
freq = [315 400 500 630 800 1000 1250 1600];

width = 0.38
height = 0.38
llw = [0.1 0.6 0.1 0.6];
llh = [0.6 0.1];

gca1 = figure(1);
set(hf2, 'renderer', 'painters')

k = 1;
for ifreq = 2:2:nfreq
  if k <= 2
    subplot('Position',[llw(k) llh(1) width height]);
  else
    subplot('Position',[llw(k) llh(2) width height]);
  end
  hold on;
  pcolor(F(ifreq).csave);
  shading interp;
  colormap('gray');
  set(h2,'layer','top')
  box on;
  set(gca,'YDir','reverse','FontSize',12);
  if k <= 2
    set(gca,'XTick',[0 20 40  60],'XTickLabel',{'0'; '20'; '40'; '60'})
    set(gca,'YTick',[0 20 40  60],'YTickLabel',{'0'; '20'; '40'; '60'})
  else
    set(gca,'XTick',[0 40 80 120],'XTickLabel',{'0'; '40'; '80'; '120'})
    set(gca,'YTick',[0 40 80 120],'YTickLabel',{'0'; '40'; '80'; '120'})
  end
  if k >= 3; xlabel('Data points','FontSize',12);end;
  if k == 1;ylabel('Data points','FontSize',12);end;
  if k == 3;ylabel('Data points','FontSize',12);end;
  set(gca,'XLim',[1 length(F(ifreq).csave)],'YLim',[1 length(F(ifreq).csave)]);
  text(4.5*length(F(ifreq).csave)/6,length(F(ifreq).csave)/6,...
      [num2str(freq(ifreq)) ' Hz'],'FontSize',10,'BackgroundColor',[1 1 1]);

  k = k +1;

end

saveas(gca1,plotfile,'eps2')
%saveas(gca1,'CAA_plot1.tif','tiffn')
return;
