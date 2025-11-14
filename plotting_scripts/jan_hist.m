function h = jan_hist(x,bins);

[n1,lim]=hist(x,bins);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
hold on; box on;set(gca,'FontSize',14);
n1 = n1/sum(n1);
lim = lim - (lim(3)-lim(2))/2;
[xx,yy]=stairs(lim,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n1,'k');
return;