function plot_meanshade;

set(0, 'DefaultFigurePaperPosition', [0 0 11 5.5]);
set(gcf, 'renderer', 'painters')
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');


%%
%% Site 01
%%
%profilefile = 'x_s01_1_3_20_5lay_profile.mat';
%plotfile1   = 'x_s01_1_3_20_5lay_cmeanprofile.eps';
%plotfile2   = 'x_s01_1_3_20_5lay_rmeanprofile.eps';
%%
%% Site 04
%%
%profilefile = 'x_s04_1_3_16_0lay_profile.mat';
%plotfile1   = 'x_s04_1_3_16_0lay_cmeanprofile.eps';
%plotfile2   = 'x_s04_1_3_16_0lay_rmeanprofile.eps';
%%
%% Site 05
%%
%profilefile = 'x_s05_1_8_32_0lay_profile.mat';
%plotfile1   = 'x_s05_1_8_32_0lay_cmeanprofile.eps';
%plotfile2   = 'x_s05_1_8_32_0lay_rmeanprofile.eps';
%%
%% Site 19
%%
%profilefile = 'x_s19_1_8_25_0lay_profile.mat';
%plotfile1   = 'x_s19_1_8_25_0lay_cmeanprofile.eps';
%plotfile2   = 'x_s19_1_8_25_0lay_rmeanprofile.eps';
%%
%% Site 20
%%
profilefile = 'x_s20_1_8_40_0lay_profile.mat';
plotfile1   = 'x_s20_1_8_40_0lay_cmeanprofile.eps';
plotfile2   = 'x_s20_1_8_40_0lay_rmeanprofile.eps';
%%
%% Site 21
%%
%profilefile = 'x_s21_1_5_40_6lay_profile.mat';
%plotfile1   = 'x_s21_1_5_40_6lay_cmeanprofile.eps';
%plotfile2   = 'x_s21_1_5_40_6lay_rmeanprofile.eps';

load(profilefile);

NBIN    = 80;
percent = 90;
[M,NZ] = size(c);
NC = 400;
NR = 400;
NA = 400;
pmin = [1450 1.3 0.0];
pmax = [1500 1.8 1.0];
clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
cplot = pmin(1)*ones(NC,NZ)-pmin(1)/400;
rplot = pmin(2)*ones(NR,NZ)-pmin(2)/400;
aplot = pmin(3)*ones(NA,NZ)-pmin(3)/400;

for i=1:NZ;
   c_hpd_int(i,:) = hpd(c(:,i),NBIN,percent);
   r_hpd_int(i,:) = hpd(r(:,i),NBIN,percent);
   a_hpd_int(i,:) = hpd(a(:,i),NBIN,percent);
end;


for i=1:NZ;
   idx = find(c_hpd_int(i,1)<clim & c_hpd_int(i,2)>clim);
   cplot(idx,i) = c_mean(i);
   idx = find(r_hpd_int(i,1)<rlim & r_hpd_int(i,2)>rlim);
   rplot(idx,i) = r_mean(i);
   idx = find(a_hpd_int(i,1)<alim & a_hpd_int(i,2)>alim);
   aplot(idx,i) = a_mean(i);
end;

h1 = figure(1);
pcolor(clim,z,cplot');
set(gca,'Fontsize',14);
set(gca,'layer','top')
set(gca,'YDir','reverse','CLim',[pmin(1) pmax(1)]);
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
hold on;box on;
shading flat;
cmap = colormap(jet);
%cmap = colormap(flipud(bone));
cmap(1,:) = [1 1 1];
colormap(cmap);
colorbar;

h2 = figure(2);
pcolor(rlim,z,rplot');
set(gca,'Fontsize',14);
xlabel('Density (g/ccm)');
ylabel('Depth (m)');
set(gca,'layer','top')
set(gca,'YDir','reverse','CLim',[pmin(2) pmax(2)]);
hold on;box on;
shading flat;
cmap = colormap(jet);
%cmap = colormap(flipud(bone));
cmap(1,:) = [1 1 1];
colormap(cmap);
colorbar;

exportfig(h1,plotfile1,opts);
exportfig(h2,plotfile2,opts);
%saveas(h1,plotfile1,'png');
%saveas(h2,plotfile2,'png');

return;
