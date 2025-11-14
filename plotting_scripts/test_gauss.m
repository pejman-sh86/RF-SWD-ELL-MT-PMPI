%% log Normal 

C=load('lin_post_cov.txt');

xx=[0:.001:1];
mux=-2;
sigx=1;
lng = 1./(xx.*sqrt(2.*pi).*sigx).*exp(-1./(2.*sigx^2).*(log(xx)-mux).^2);

figure();
plot(xx,lng)


stop









%% Normal distribution marginal test 
mu = [0,0];
sig = [1.,.9;.9,1.];
rng default
r = mvnrnd(mu,sig,10000);
figure();
plot(r(:,1),r(:,2),'.');

x = [-10:.1:10];
figure();hold on;
[n1,out]=hist(r(:,1),x);
dx = mean(diff(out));
n1 = n1/(sum(n1)*dx);
n1 = [0, n1, 0];out = [out(1) out out(end)];
[xx,yy]=stairs(out,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(out,n1,'k');
axis tight;

xx=[-10:.01:10];
mux=0;
sigx=sqrt(sig(1,1));
g = 1./(sqrt(2.*pi)*sigx)*exp(-1./(2.*sigx^2).*(xx-mux).^2);

g2 = normpdf(xx,mux,sigx);

plot(xx,g);
plot(xx,g2,'--');
