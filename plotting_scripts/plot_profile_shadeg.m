%
% Plot 2-D profile density
%
function plot_profile_shade;

set(0, 'DefaultFigurePaperPosition', [0 0 10 5]);
set(gcf, 'renderer', 'painters')

%sample    = 'sim_C_1_7_40_6layrg_sample.mat';
%mapfile   = 'sim_C_1_7_40_6layrg_fgsmap.txt';
%plotfile1 = 'sim_C_1_7_40_6layrg_profileshade.eps';

%sample    = 'x_s01_1_3_20_2layrg_sample.mat';
%mapfile   = 'x_s01_1_3_20_2layrg_fgsmap.txt';
%plotfile1 = 'x_s01_1_3_20_2layrg_profileshade.eps';
%corefile  = 'core.mat';

%sample    = 'x_s02_1_3_16_2layrg_sample.mat';
%mapfile   = 'x_s02_1_3_16_2layrg_fgsmap.txt';
%plotfile1 = 'x_s02_1_3_16_2layrg_profileshade.eps';
%corefile  = 'core.mat';

%sample    = 'x_s13_2_10_50_5layrgb_sample.mat';
%mapfile   = 'x_s13_2_10_50_5layrgb_fgsmap.txt';
%plotfile1 = 'x_s13_2_10_50_5layrgb_profileshade.eps';
%corefile  = 'core13.mat';

%sample    = 'x_s16_1_5_40_6layrgcov_sample.mat';
%mapfile   = 'x_s16_1_5_40_6layrgcov_fgsmap.txt';
%plotfile1 = 'x_s16_1_5_40_6layrgcov_profileshade.eps';
%corefile  = 'core15.mat';

sample    = 'x_s21_1_5_25_4layrg_sample.mat';
mapfile   = 'x_s21_1_5_25_4layrg_fgsmap.txt';
plotfile1 = 'x_s21_1_5_25_4layrg_profileshade.eps';
corefile  = 'core.mat';

icore = 1;
isyn = 0;
iplot = 1;
ifull = 0;
Tstar = 1;
if(isyn == 1)
   mtru = [0.31, 1480, 1.30, 1.6, 0.01, ...
            0.35, 1510, 1.65, 0.30, ...
            1.20, 1500, 1.55, 0.30, ...
            0.40, 1700, 1.90, 0.30, ...
                  1560, 1.70, 0.01 ];
   mtrutmp = [mtru(1:3) mtru(5:end)];
end;
nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',6,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

NZ = 140;
NC = 120;
NR = 100;
NA = 80;
NPARL = 4;
hmax = 4;

pmin = [1450 1.2 0.0];
pmax = [1850 2.2 1.0];

clim = pmin(1)+cumsum((pmax(1)-pmin(1))/(NC-1)*ones(1,NC))-((pmax(1)-pmin(1))/(NC-1));
rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/(NR-1)*ones(1,NR))-((pmax(2)-pmin(2))/(NR-1));
alim = pmin(3)+cumsum((pmax(3)-pmin(3))/(NA-1)*ones(1,NA))-((pmax(3)-pmin(3))/(NA-1));

dz = hmax/(NZ-1);
z = cumsum(dz*ones(1,NZ))-dz;

load(sample);
NLAY  = (size(A,2)-5)/4
map = load(mapfile);
maptmp = [map(1:3) map(5:end)];
if(icore == 1)
    load(corefile);
end;

m = A(:,2:end);
E = A(:,1);
NSAMP = length(m);
if(ifull == 0)
    NPROF = 3000;
    m = m(round(NSAMP*rand(NPROF,1)),:);
else
    NPROF = NSAMP;
end

if(NLAY>1)
   idxh = (([1:NLAY-1])*NPARL)+2;
   idxc = (([1:NLAY-1])*NPARL)+3;
   idxr = (([1:NLAY-1])*NPARL)+4;
   idxa = (([1:NLAY-1])*NPARL)+5;
   idxh = [1 idxh idxh(end)];
   idxc = [2 idxc idxc(end)+3];
   idxr = [3 idxr idxr(end)+3];
   idxa = [5 idxa idxa(end)+3];
else
   idxh = [1];
   idxc = [2 6];
   idxr = [3 7];
   idxa = [5 8];
end;

disp('Starting profiles...');
disp('Sample size: '),disp(size(m))

prof(1:NPROF,:,1) = cumsum(m(1:NPROF,idxh),2);
prof(1:NPROF,:,2) = m(1:NPROF,idxc);
prof(1:NPROF,:,3) = m(1:NPROF,idxr);
prof(1:NPROF,:,4) = m(1:NPROF,idxa);

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
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
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
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',1);
end;
plprof(maptmp,hmax,'k',1);
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
box on;

subplot(h2)
Nr = [Nr Nr(:,end)];
pcolor(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
   plprofg(mtru,hmax,'--k',2);
end;
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
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',3);
end;
plprof(maptmp,hmax,'k',3);
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;

barstep1 = 5;
barstep1r = 1;
barstep2 = 2;
barstep2r = 1;
barstep3 = 10;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'k','Linewidth',1);
%    plot(c2(:,2),c2(:,1),'r','Linewidth',2);

    subplot(h2);
    plot(r1(:,2),r1(:,1),'k','Linewidth',1);
%    plot(r2(:,2),r2(:,1),'r','Linewidth',2);

%    subplot(h1);
%    plot(c1(:,2),c1(:,1),'g','Linewidth',1);
%    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
%               10.*ones(size(c1(1:barstep1:end,2))),...
%               zeros(size(c1(1:barstep1:end,2))),'g','g');
%    plot(c2(:,2),c2(:,1),'b','Linewidth',1);
%    errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%               2*10.*ones(size(c2(1:barstep2:end,2))),...
%               zeros(size(c2(1:barstep2:end,2))),'b','b');
%    plot(c3(:,2),c3(:,1),'r','Linewidth',1);
%    errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%               10.*ones(size(c3(1:barstep3:end,2))),...
%               zeros(size(c3(1:barstep3:end,2))),'r','r');
%
%
%    subplot(h2);
%    plot(r1(:,2),r1(:,1),'g','Linewidth',1);
%    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
%               2./100.*r1(1:barstep1r:end,2),...
%               zeros(size(r1(1:barstep1r:end,2))),'g','g');
%    plot(r2(:,2),r2(:,1),'b','Linewidth',1);
%    errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%               4./100.*r2(1:barstep2r:end,2),...
%               zeros(size(r2(1:barstep2r:end,2))),'b','b');
%    plot(r3(:,2),r3(:,1),'r','Linewidth',1);
%    errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%               2./100.*r3(1:barstep3:end,2),...
%               zeros(size(r3(1:barstep3:end,2))),'r','r');


end;

if(iplot == 1)
%    saveas(hf2,plotfile1,'epsc2');
    exportfig(hf2,plotfile1,opts);
end;

return;
