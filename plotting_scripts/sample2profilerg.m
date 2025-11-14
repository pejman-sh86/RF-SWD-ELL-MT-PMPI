%
% Plot 2-D profile density
%
function plot_profile_shade;

set(0, 'DefaultFigurePaperPosition', [0 0 10 5]);
set(gcf, 'renderer', 'painters')

sample    = 'x_s13_2_10_50_5layrg_sample.mat';
mapfile   = 'x_s13_2_10_50_5layrg_map.txt';
plotfile1 = 'x_s13_2_10_50_5layrg_profileshade.eps';
corefile  = 'core13.mat';
%sample    = 'x_s21_1_10-50_4layrg_b_ci_samplec.mat';
%mapfile   = 'x_s21_1_10-50_4layrg_b_ci_map.txt';
%plotfile1 = 'x_s21_1_10-50_4layrg_b_ci_profileshade.eps';
%corefile  = 'core03.mat';
%sample    = 'x_s16_1_15_40_3layrg_sample.mat';
%mapfile   = 'x_s16_1_15_40_3layrg_map.txt';
%plotfile1 = 'x_s16_1_15_40_3layrg_profileshade.eps';
%corefile  = 'core15.mat';
%sample    = 'tmp_sample.mat';
%mapfile   = 'tmp_ci_map.txt';
%plotfile1 = 'tmp_ci_profileshade.eps';
%corefile  = 'core15.mat';

icore = 1;
iplot = 1;
ifull = 0;
Tstar = 1;

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',6,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

NZ = 100;
NC = 60;
NR = 60;
NA = 60;
NLAY  = 5;
NPARL = 4;
hmax = 8;

pmin = [1450 1.2 0.0];
pmax = [1750 2.4 1.0];

clim = pmin(1)+cumsum((pmax(1)-pmin(1))/(NC-1)*ones(1,NC))-((pmax(1)-pmin(1))/(NC-1));
rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/(NR-1)*ones(1,NR))-((pmax(2)-pmin(2))/(NR-1));
alim = pmin(3)+cumsum((pmax(3)-pmin(3))/(NA-1)*ones(1,NA))-((pmax(3)-pmin(3))/(NA-1));

dz = hmax/(NZ-1);
z = cumsum(dz*ones(1,NZ))-dz;

load(sample);
map = load(mapfile);
maptmp = [map(1:3) map(5:end)];
if(icore == 1)
    load(corefile);
end;

m = A(:,2:end);
E = A(:,1);
NSAMP = length(m);
if(ifull == 0)
    NPROF = 160000;
    m = m(round(NSAMP*rand(NPROF,1)),:);
else
    NPROF = NSAMP;
end

idxh = (([1:NLAY-1])*NPARL)+2;
idxc = (([1:NLAY-1])*NPARL)+3;
idxr = (([1:NLAY-1])*NPARL)+4;
idxa = (([1:NLAY-1])*NPARL)+5;
idxh = [1 idxh idxh(end)];
idxc = [2 idxc idxc(end)+3];
idxr = [3 idxr idxr(end)+3];
idxa = [5 idxa idxa(end)+3];

disp('Starting profiles...');

prof(1:NPROF,:,1) = cumsum(m(1:NPROF,idxh),2);
prof(1:NPROF,:,2) = m(1:NPROF,idxc);
prof(1:NPROF,:,3) = m(1:NPROF,idxr);
prof(1:NPROF,:,4) = m(1:NPROF,idxa);

disp('Sample size: '),disp(size(m))
idx=1;
for iprof = 1:NPROF

    if(iprof==idx*1000)
        fprintf(1,'%8i',iprof)
        idx=idx+1;
    end
    izold = 1;
    for ilay=1:NLAY
        idxz = find(z <= prof(iprof,ilay,1));
        c(iprof,idxz(izold:end)) = prof(iprof,ilay,2);
        a(iprof,idxz(izold:end)) = prof(iprof,ilay,4);
        if(ilay>1)
            r(iprof,idxz(izold:end)) = prof(iprof,ilay,3);
        else
            for iz=1:idxz(end)
                r(iprof,iz) = prof(iprof,1,3)+(m(iprof,4)-prof(iprof,1,3))/...
                              prof(iprof,1,1)*z(iz);
            end
        end
        izold = idxz(end)+1;
    end;
    c(iprof,izold:NZ) = prof(iprof,end,2);
    r(iprof,izold:NZ) = prof(iprof,end,3);
    a(iprof,izold:NZ) = prof(iprof,end,4);
end;
fprintf(1,'\n')
disp('Done with profiles.');

disp('Starting histograms...');
for iz=1:NZ
    Nc(iz,:) = histc_tstar(c(:,iz),E,clim,Tstar);
    Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
    Nr(iz,:) = histc_tstar(r(:,iz),E,rlim,Tstar);
    Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
    Na(iz,:) = histc_tstar(a(:,iz),E,alim,Tstar);
    Na(iz,:) = Na(iz,:)/max(Na(iz,:));
%    disp(iz);
end;
disp('Done histograms.');

nx = 3;
ny = 1;
xim = 0.06;
yim = 0.06;
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure(1);hold on; box on;
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box on;
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.8 2.2]);
set(gca,'XTick',[1.4 1.8 2.2]);
hold on; box on;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[0.0 0.5 1.0]);
set(gca,'XTick',[0.0 0.5 1.0]);
hold on; box on;

subplot(h1)
Nc = [Nc Nc(:,end)];
pcolor(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
plprof(maptmp,hmax,'k',1);
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
box on;

subplot(h2)
Nr = [Nr Nr(:,end)];
pcolor(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
plprofg(map,hmax,'k',2);
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
box on;

subplot(h3)
Na = [Na Na(:,end)];
pcolor(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)');
plprof(maptmp,hmax,'k',3);
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;


if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'k','Linewidth',1);
%    plot(c2(:,2),c2(:,1),'r','Linewidth',2);

    subplot(h2);
    plot(r1(:,2),r1(:,1),'k','Linewidth',1);
%    plot(r2(:,2),r2(:,1),'r','Linewidth',2);
end;

if(iplot == 1)
%    saveas(hf2,plotfile1,'epsc2');
    exportfig(hf2,plotfile1,opts);
end;

return;
