function []=plot_evidence_ais();

data = load('evidence_with_power.dat');
N = 20;
npar = data(1,:);
logZ = data(2,:);
logL = -data(3,:);
[npar,idx] = sort(npar)
logZ = logZ(idx);
logL = logL(idx);
para = [{'1 Layer'},{'2 Layer'},{'3 Layer'},{'Lin. Gradient'},{'Power Law'}];
%para = [{'1 Layer'},{'2 Layer'},{'3 Layer'},{'Lin. Gradient'}];
para = para(idx);

BIC = -2.*logL+npar*log(N);
BIC
BIC = BIC - min(BIC);

figure(1);
subplot(4,1,1);hold on;box on;
plot(logZ,'--xk','MarkerSize',10)
plot(3,logZ(3),'or','MarkerSize',12,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('log(Z)');
set(gca,'YLim',[-54 -26],'XLim',[0.8 5.2])
set(gca,'XTick',[1,2,3,4,5])
%set(gca,'XTickLabel',para)
set(gca,'XTickLabel',[])

subplot(4,1,2);hold on;box on;
plot(logL,'--xk','MarkerSize',10);
%plot(2,logL(2),'ob','MarkerSize',12,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('log(L)');
set(gca,'YLim',[-13 16],'XLim',[0.8 5.2])
set(gca,'XTick',[1,2,3,4,5])
%set(gca,'XTickLabel',para)
set(gca,'XTickLabel',[])

subplot(4,1,3);hold on;box on;
plot(BIC,'--xk','MarkerSize',10)
%plot(2,BIC(2),'ob','MarkerSize',12,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('BIC');
set(gca,'YLim',[0 51],'XLim',[0.8 5.2])
set(gca,'XTick',[1,2,3,4,5])
%set(gca,'XTickLabel',para)
set(gca,'XTickLabel',[])

subplot(4,1,4);hold on;box on;
plot(npar,'--xk','MarkerSize',10);
set(gca,'FontSize',14)
ylabel('No. Parameters');
set(gca,'YLim',[7 17],'XLim',[0.8 5.2])
set(gca,'XTick',[1,2,3,4,5])
set(gca,'XTickLabel',para)

return;
