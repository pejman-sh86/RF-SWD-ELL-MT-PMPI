function [] = plot_pvalues();

set(0, 'DefaultFigurePaperPosition', [0 0 7 3]);
pfile = 'x_jan_2_sph_ciw_pvalues.txt';
plotfile1 = 'x_jan_2_sph_ciw_pvalues_runs.eps';
plotfile2 = 'x_jan_2_sph_ciw_pvalues_ks.eps';
%pfile = 'cov_offstat_pvalues.txt';
%plotfile1 = 'cov_offstat_pvalues_runs.eps';
%plotfile2 = 'cov_offstat_pvalues_ks.eps';

%freq = [315 400 500 630 800 1000 1250 1600];
nfreq = 5;
freq = [400 600 900 1250 1600];
frtck = [1 2 3 4 5];
%nfreq = 7;
%freq  = [1 2 3 4 5 6 7];
%frtck = [1 2 3 4 5 6 7];
p = load(pfile)

figure(1);box on;
semilogy([0 nfreq+1],[0.05 0.05],'--k');
hold on;
semilogy(p(1:nfreq,1),'x--k');
semilogy(p(nfreq+2:end,1),'o:k')
xlabel('freq','FontSize',14);
ylabel('p','FontSize',14);
set(gca,'FontSize',14,'XLim',[0 9]);
set(gca,'XTick',frtck,'XTickLabel',freq);
set(gca,'XLim',[0 frtck(end)+1]);
saveas(gca,plotfile1,'epsc2');

figure(2);box on;
semilogy([0 nfreq+1],[0.05 0.05],'--k')
hold on;
semilogy(p(1:nfreq,2),'x--k')
semilogy(p(nfreq+2:end,2),'o:k')
xlabel('freq','FontSize',14);
ylabel('p','FontSize',14);
set(gca,'FontSize',14,'XLim',[0 9]);
set(gca,'XTick',frtck,'XTickLabel',freq);
set(gca,'XLim',[0 frtck(end)+1]);
saveas(gca,plotfile2,'epsc2');

return;
