function []=plot_n();

load ../band1/x_sim_1_1_2_biot_sample.mat;
A1 = A;
afdep1 = A1(:,49);
load ../band2/x_sim_1_3_8_biot_sample.mat;
A2 = A;
afdep2 = A2(:,49);
load ../band3/x_sim_1_10_25_biot_sample.mat;
A3 = A;
afdep3 = A3(:,49);
load ../band4/x_sim_1_31_100_biot_sample.mat;
A4 = A;
afdep4 = A4(:,49);

fig1=figure;
nx = 2;
ny = 2;
xim = 0.01;
yim = 0.01;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

subplot('Position',[loc(1,1) loc(2,1) spw sph]);
set(gca,'FontSize',16);
hold on; box on;
[n,lim]=hist(afdep1,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
ylabel('Probability');
set(gca,'Layer','top');
set(gca,'XLim',[0.0 2.0]);
set(gca,'YLim',[0 0.06],'YTick',[0,0.02,0.04],'XTick',[0,0.5,1.0,1.5,2.0]);
set(gca,'XTickLabel',[]);
subplot('Position',[loc(1,1) loc(2,1) spw sph]);
text(0.1,0.05,'100-250 Hz','FontSize',14)
box on;

subplot('Position',[loc(1,2) loc(2,2) spw sph]);
set(gca,'FontSize',16);
hold on; box on;
[n,lim]=hist(afdep2,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
set(gca,'Layer','top');
set(gca,'XLim',[0.0 2.0]);
set(gca,'YLim',[0 0.06],'YTick',[0,0.02,0.04],'XTick',[0,0.5,1.0,1.5,2.0]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot('Position',[loc(1,2) loc(2,2) spw sph]);
text(0.1,0.05,'300-800 Hz','FontSize',14)
box on;

subplot('Position',[loc(1,3) loc(2,3) spw sph]);
set(gca,'FontSize',16);
hold on; box on;
[n,lim]=hist(afdep3,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('Attenuation exponent n');
ylabel('Probability');
set(gca,'Layer','top');
set(gca,'XLim',[0.0 2.0]);
set(gca,'YLim',[0 0.06],'YTick',[0,0.02,0.04],'XTick',[0,0.5,1.0,1.5,2.0]);
subplot('Position',[loc(1,3) loc(2,3) spw sph]);
text(0.1,0.05,'1000-2,500 Hz','FontSize',14)
box on;

subplot('Position',[loc(1,4) loc(2,4) spw sph]);
set(gca,'FontSize',16);
hold on; box on;
[n,lim]=hist(afdep4,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('Attenuation exponent n');
set(gca,'Layer','top');
set(gca,'XLim',[0.0 2.0]);
set(gca,'YLim',[0 0.06],'YTick',[0,0.02,0.04],'XTick',[0,0.5,1.0,1.5,2.0]);
set(gca,'YTickLabel',[]);
text(0.1,0.05,'3150-10,000 Hz','FontSize',14)
subplot('Position',[loc(1,4) loc(2,4) spw sph]);
box on;

return;
