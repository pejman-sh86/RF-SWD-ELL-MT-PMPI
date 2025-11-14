%function [] = auv_plot_bathy_profile(filename);

filename = 'AUV_May22_NAV.mat';
filebase = strrep(filename,'.mat','');

i_inter = 0;  %% 1 = resample onto even grid
i_frave = 1;  %% 1 = average frequencies 
i_pave  = 1;  %% 1 = apply ping averaging
nskip = 2;  %% stride between pings
npave = 5;  %% stride between pings

pinglist = load('ping_list.txt');
NPING = pinglist(1);
dping = 300;      %% Interval for ping axis label
pinglist(1)=[];
pstart = pinglist(1);
pend = pinglist(end);

nang_re = 31;
load(filename);
load track_cra.mat

%% Start end end of track estimated from Charles' email.
ist = 1746;
iend = ist + 100;
%iend = ist + 107;
%iend = 12044;

%% Bathy is samples every second; source fires every 3 seconds.
bathy=Abathy(ist:iend);
bathy2=bathy(1:2:end);

lat=Alat(ist:iend);
lat2=lat(1:2:end);

time = auv_t(ist:iend);
time2= time(1:2:end);

%% Frequencies:
freq = [988., 1113., 1288., 1913., 2263., 2513.];

idat = 1;
for iping = pstart:nskip:pend;

   bathy_ave(idat) = mean(bathy2(iping-npave:iping+npave));
   lat_ave(idat) = mean(lat2(iping-npave:iping+npave));

   time_ave(idat)  = mean(time2(iping-npave:iping+npave));

   idx(idat) = iping;
   idat = idat + 1;
end

bathy_rel = bathy_ave-min(bathy_ave);
bathy_mag = max(bathy_rel)-min(bathy_rel);

NZI = 500;
NZ  = 300;
hmax = 7.5;
dz  = hmax/(NZ-1);
z   = cumsum(dz*ones(1,NZ))-dz;
dzi = hmax/(NZI-1);
zi  = cumsum(dzi*ones(1,NZI))-dzi;

znew = [0:dz:max(z)+bathy_mag]+min(bathy_ave);
c_mean_new = zeros(length(znew),1468)+1511;
c_med_new = zeros(length(znew),1468)+1511;

r_mean_new = zeros(length(znew),339)+1511;
r_med_new = zeros(length(znew),339)+1511;

a_mean_new = zeros(length(znew),339)+1511;
a_med_new = zeros(length(znew),339)+1511;

%hgl = hgl';
znewi = [0:dzi:max(zi)+bathy_mag]+min(bathy_ave);
hgl_new = zeros(length(znewi),339);

save tmp.mat bathy_ave bathy_rel time_ave;

idx = floor(bathy_rel/dz);
idxi = floor(bathy_rel/dzi);

for iping = 1:5;
    c_mean_new(1:idx(iping),iping) = min(min(c_mean))-10;
    c_mean_new([1:250]+idx(iping),iping) = c_mean(:,iping);
    c_mean_new(250+idx(iping):end,iping) = min(min(c_mean))-10;

    c_med_new(1:idx(iping),iping) = min(min(c_mead))-10;
    c_med_new([1:250]+idx(iping),iping) = c_mead(:,iping);
    c_med_new(250+idx(iping):end,iping) = min(min(c_mead))-10;

    r_mean_new(1:idx(iping),iping) = min(min(r_mean))-.1;
    r_mean_new([1:250]+idx(iping),iping) = r_mean(:,iping);
    r_mean_new(250+idx(iping):end,iping) = min(min(r_mean))-.1;

    r_med_new(1:idx(iping),iping) = min(min(r_mead))-.1;
    r_med_new([1:250]+idx(iping),iping) = r_mead(:,iping);
    r_med_new(250+idx(iping):end,iping) = min(min(r_mead))-.1;

    a_mean_new(1:idx(iping),iping) = min(min(a_mean))-.1;
    a_mean_new([1:250]+idx(iping),iping) = a_mean(:,iping);
    a_mean_new(250+idx(iping):end,iping) = min(min(a_mean))-.1;

    a_med_new(1:idx(iping),iping) = min(min(a_mead))-.1;
    a_med_new([1:250]+idx(iping),iping) = a_mead(:,iping);
    a_med_new(250+idx(iping):end,iping) = min(min(a_mead))-.1;

    hgl_new([1:400]+idxi(iping),iping) = hgl(:,iping);
end

%save tmp.mat;

figw = 18;
figh = 8;
nx = 1;
ny = 3;
ML = .05;
MR = .03;
MB = .08;
MT = .02;
SP = .02;
PAD = 0;
FNT = 12;
inter = 20;

fig_mean=figure();hold on;
set(fig_mean,'PaperUnits','inches','PaperPosition',[0 0 18 8]);
subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,c_mean_new);
set(gca, 'Layer', 'top')
mycolormap = [ones(1,3); jet(128)];
colormap(mycolormap);
%xlabel('Latitude');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Velocity (m/s)');
box on;

subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,r_mean_new);
set(gca, 'Layer', 'top')
mycolormap = [ones(1,3); jet(128)];
colormap(mycolormap);
xlabel('Latitude N');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
%set(gca,'YTickLabel',[]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Density (g/ccm)');
box on;

subaxis(ny,nx,3,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,a_mean_new);
set(gca, 'Layer', 'top')
mycolormap = [ones(1,3); jet(128)];
colormap(mycolormap);
xlabel('Latitude N');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
%set(gca,'YTickLabel',[]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Density (g/ccm)');
box on;

fig_med=figure();hold on;
set(fig_med,'PaperUnits','inches','PaperPosition',[0 0 18 8]);
subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,c_med_new);
set(gca, 'Layer', 'top')
mycolormap = [ ones(1,3); jet(128)];
colormap(mycolormap);
%xlabel('Latitude');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Velocity (m/s)');
box on;

subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,r_med_new);
mycolormap = [ ones(1,3); jet(128)];
set(gca, 'Layer', 'top')
colormap(mycolormap);
xlabel('Latitude N');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Density (g/ccm)');
box on;

subaxis(ny,nx,3,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znew,a_med_new);
mycolormap = [ ones(1,3); jet(128)];
set(gca, 'Layer', 'top')
colormap(mycolormap);
xlabel('Latitude N');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
set(gca,'XTick',[36:.01:37]);
c1=colorbar;ylabel(c1,'Density (g/ccm)');
box on;

fig_int=figure();hold on;
set(fig_int,'PaperUnits','inches','PaperPosition',[0 0 18 8]);
subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);hold on;
imagesc(lat_ave,znewi,hgl_new);
set(gca, 'Layer', 'top')
mycolormap = [ones(1,3); jet(128)];
colormap(flipud(gray));
xlabel('Latitude N');
ylabel('Depth (m)');
plot(lat_ave,bathy_ave,'-k','LineWidth',2);
plot(lat_ave,bathy_ave,'--w','LineWidth',2);
set(gca,'XLim',[lat_ave(end),lat_ave(1)],'FontSize',FNT);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse','TickDir','out');
set(gca,'YTick',[0:5:500],'Ylim',[142,160]);
set(gca,'XTick',[36:.01:37]);
set(gca,'CLim',[0 0.1]);
c1=colorbar;ylabel(c1,'Interface probability');
box on;

print(fig_mean,'-r300','bathy_results_mean.png','-dpng');
print(fig_med,'-r300','bathy_results_med.png','-dpng');
print(fig_int,'-r300','bathy_results_interface.png','-dpng');

saveas(fig_mean,'bathy_results_mean.fig','fig');
saveas(fig_med,'bathy_results_med.fig','fig');
saveas(fig_int,'bathy_results_interface.fig','fig');

%return;