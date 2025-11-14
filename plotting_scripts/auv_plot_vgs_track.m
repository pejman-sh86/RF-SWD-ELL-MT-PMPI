load track.mat
z_por = z;
load track_cra.mat

bad_pings = [5,9];

for iping=1:length(bad_pings);
    meact(bad_pings(iping),:,:) = NaN;
    meart(bad_pings(iping),:,:) = NaN;
    meaat(bad_pings(iping),:,:) = NaN;
    nfct(bad_pings(iping),:,:,:) = NaN;
    nfrt(bad_pings(iping),:,:,:) = NaN;
    nfat(bad_pings(iping),:,:,:) = NaN;
end;

x = 4.*pinglist;
ndep = 200;
ndep2 = 200;
NP = 106;
FNT = 22;

fig1=figure(1);
nx = 1;
ny = 3;
xim = 0.01;
yim = 0.02;
xymarg = [0.06 0.04 0.012 0.085];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
figw = 20;
figh = 12;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 figw figh];
set(fig1,'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
%subplot(3,1,1)
pcolor(x(1:NP),z(1:ndep),squeeze(nfct(1:NP,1,1:ndep,2))'-squeeze(nfct(1:NP,1,1:ndep,1))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[0,150],'TickDir','out','FontSize',FNT);
shading flat;box on;
%ylabel('Depth below seafloor (m)');
set(gca,'XTickLabel',[]);
h = colorbar;
ylabel(h, '95% CI Velocity (m/s)','FontSize',FNT-2)

h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
%subplot(3,1,2)
pcolor(x(1:NP),z(1:ndep),squeeze(nfrt(1:NP,1,1:ndep,2))'-squeeze(nfrt(1:NP,1,1:ndep,1))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[0,1],'TickDir','out','FontSize',FNT);
shading flat;box on;
ylabel('Depth below seafloor (m)');
set(gca,'XTickLabel',[]);
h = colorbar;
ylabel(h, '95% CI Density (g/ccm)','FontSize',FNT-2)

h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
pcolor(x(1:NP),z(1:ndep),squeeze(nfat(1:NP,1,1:ndep,2))'-squeeze(nfat(1:NP,1,1:ndep,1))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[0,2],'TickDir','out','FontSize',FNT);
shading flat;box on;
%ylabel('Depth below seafloor (m)');
xlabel('Distance (m)')
h = colorbar;
ylabel(h, '95% CI Atten. log_{10}(dB/m/kHz)','FontSize',FNT-2)

fig2=figure(2);
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 figw figh];
set(fig2, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
colormap jet;
pcolor(x(1:NP),z(1:ndep),squeeze(meact(1:NP,1,1:ndep))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[1450,1600],'TickDir','out','FontSize',FNT);
shading flat;box on;
%ylabel('Depth below seafloor (m)');
%xlabel('Distance (m)')
h = colorbar;
ylabel(h, 'Velocity (m/s)','FontSize',FNT-2)

h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
colormap jet;
pcolor(x(1:NP),z(1:ndep),squeeze(meart(1:NP,1,1:ndep))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[1.20,2.00],'TickDir','out','FontSize',FNT);
shading flat;box on;
ylabel('Depth below seafloor (m)');
xlabel('Distance (m)')
h = colorbar;
ylabel(h, 'Density (g/ccm)','FontSize',FNT-2)

h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
colormap jet;
pcolor(x(1:NP),z(1:ndep),squeeze(meaat(1:NP,1,1:ndep))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[-2.35,0.8],'TickDir','out','FontSize',FNT);
shading flat;box on;
%ylabel('Depth below seafloor (m)');
xlabel('Distance (m)')
h = colorbar;
ylabel(h, 'Atten. log_{10}(dB/m/kHz)','FontSize',FNT-2)

fig3=figure(3);
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 figw figh];
set(fig3, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
colormap jet;
pcolor(x(1:NP),z_por(1:ndep2),c_mean(1:ndep2,1:NP));
colormap(jet);
set(gca,'YDir','reverse','CLim',[0.2,1.],'TickDir','out','FontSize',FNT);
shading flat;box on;
ylabel('Depth below seafloor (m)');
xlabel('Distance (m)')
h = colorbar;
ylabel(h, 'Porosity','FontSize',FNT-2)

fig4=figure(4);
fig4.PaperUnits = 'inches';
fig4.PaperPosition = [0 0 figw figh];
set(fig4, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
colormap jet;
pcolor(x(1:NP),z(1:ndep2),squeeze(meact(1:NP,3,1:ndep))');
colormap(jet);
set(gca,'YDir','reverse','CLim',[1450,1600],'TickDir','out','FontSize',FNT-4);
shading flat;box on;
ylabel('Depth below seafloor (m)');
xlabel('Distance (m)')
h = colorbar;
ylabel(h, 'Velocity (m/s)','FontSize',FNT-6)

saveas(fig1,'track_cra_meadian.png','png');
saveas(fig2,'track_cra_CI.png','png');
saveas(fig3,'track_por.png','png');
saveas(fig4,'track_vel.png','png');