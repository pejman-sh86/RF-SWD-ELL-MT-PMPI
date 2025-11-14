
set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);
load cov_mat_est;
G = F;
load blde4_syn.mat
%load bla.mat
load residuals.mat
load ml_std.mat

ml_sd = stdv_dB(:,2)

L = chol(G(8).csave);
Axx1 = xcov(B(8).diff,'coeff');
Axx2 = xcov(inv(L')*B(8).diff,'coeff');

width = 0.35
height = 0.35
llw = [0.1 0.55 0.1 0.55];
llh = [0.55 0.1];

figure(1);
set(0, 'DefaultFigurePaperPosition', [0 0 6 8]);
subplot('Position',[llw(1) llh(1) width height])
hold on;box on
Err = sqrt(G(8).csave(1,1))*ones(length(F(8).dat(1:3:end)),1);
plot(F(8).ang(1:3:end), F(8).dat(1:3:end),'k.');
errorbar(F(8).ang(1:3:end),F(8).dat(1:3:end),Err,'k.');
plot(F(8).ang,Frep(8).dat,'-k')
xlabel('Angle [deg]','FontSize',14);
ylabel('BL [dB]','FontSize',14);
axis([0 90 0 40]); set(gca,'Xtick',[0:30:90],'FontSize',14);

subplot('Position',[llw(2) llh(1) width height])
hold on;box on
plot(F(8).ang(1:3:end), B(8).diff(1:3:end)/ml_sd(8),'kx');
axis([0 90 -5 5]); set(gca,'Xtick',[0:30:90],'FontSize',14);
xlabel('Angle [deg]','FontSize',14);
ylabel('Residual [\sigma]','FontSize',14);
set(gca,'YGrid','on');

subplot('Position',[llw(1) llh(2) width height])
hold on; box on
plot([-130:1:130],Axx1,'-k.')
axis([-130 130 -.5 1]);set(gca,'YTick',[0 1],'FontSize',14);
xlabel('Lag','FontSize',14);
ylabel('c','FontSize',14);
subplot('Position',[llw(2) llh(2) width height])
hold on;box on
plot([-130:1:130],Axx2,'-k.')
axis([-130 130 -.5 1]);set(gca,'YTick',[0 1],'FontSize',14);
xlabel('Lag','FontSize',14);
ylabel('c','FontSize',14);

saveas(gca,'CAA_plot2.eps','epsc2')

return;
