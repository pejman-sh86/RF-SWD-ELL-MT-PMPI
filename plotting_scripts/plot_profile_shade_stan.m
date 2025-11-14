%
% Plot 2-D profile density
%
function plot_profile_shade;

set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);
set(gcf, 'renderer', 'painters')

icore = 1;
imap = 0;
imarg = 1;
isyn = 0;
iplot = 1;
ifull = 0;

sample      = 'x_s01_1_3_20_5lay_sample.mat';
mapfile     = 'x_s01_1_3_20_5lay_fgsmap.txt';
profilefile = 'x_s01_1_5_20_5lay_profile.mat';
if(imarg == 1)
   plotfile1 = 'x_s01_1_3_20_5lay_profileshade.eps';
else
   plotfile1 = 'x_s01_1_3_20_5lay_profilemap.eps';
end
corefile  = 'core.mat';

if(isyn == 1)
%
% Sim B
%
   mtru = [0.10, 1490, 1.45, 0.30, ...
           0.30, 1550, 1.75, 0.25, ...
           0.35, 1580, 1.95, 0.30, ...
           3.50, 1540, 1.70, 0.10, ...
                 1600, 1.75, 0.50 ];
   mtrutmp = mtru;
end;

%%
%% Arrange panels for subplots:
%%
nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');


NPARL = 4;	% # par per layer
hmax = 4.0;	% max depth plotted

%%
%% plot limits
%%
pmin = [1450 1.2 0.0];
pmax = [1900 2.4 1.0];

%
% # bins in z and parameter direction 
%
if(imarg == 1)
   NZ = 400;
   NC = 160;
   NR = 160;
   NA = 160;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

   dz = hmax/(NZ-1);
   z = cumsum(dz*ones(1,NZ))-dz;
end

if(imap == 1)
   map = load(mapfile);
end;
if(icore == 1)
    load(corefile);
end;

if(imarg == 1)
   load(sample);
   NLAY  = (size(A,2)-4)/4
   m = A(:,2:end);
   NSAMP = length(m);

%
% Set # profiles NPROF (full sample might be too much)
%
   if(ifull == 0)
       NPROF = 8000;
       m = m(round(NSAMP*rand(NPROF,1)),:);
   else
       NPROF = NSAMP;
   end;

idxh = (([1:NLAY]-1)*NPARL)+1;
idxc = (([1:NLAY]-1)*NPARL)+2;
idxr = (([1:NLAY]-1)*NPARL)+3;
idxa = (([1:NLAY]-1)*NPARL)+4;
idxh = [idxh idxh(end)];
idxc = [idxc idxc(end)+3];
idxr = [idxr idxr(end)+3];
idxa = [idxa idxa(end)+3];

disp('Starting profiles...');

%
% Sort model paramerters into profiles of absolute depth
%
prof(1:NPROF,:,1) = cumsum(m(1:NPROF,idxh),2);
prof(1:NPROF,:,2) = m(1:NPROF,idxc);
prof(1:NPROF,:,3) = m(1:NPROF,idxr);
prof(1:NPROF,:,4) = m(1:NPROF,idxa);

%
% Compute discretized form of profile
%
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
        r(iprof,idxz(izold:end)) = prof(iprof,ilay,3);
        a(iprof,idxz(izold:end)) = prof(iprof,ilay,4);
        izold = idxz(end)+1;
    end;
    c(iprof,izold:NZ) = prof(iprof,end,2);
    r(iprof,izold:NZ) = prof(iprof,end,3);
    a(iprof,izold:NZ) = prof(iprof,end,4);
end;
fprintf(1,'\n')
disp('Done with profiles.');

%
% Compute histograms for each depth
%
disp('Starting histograms...');
for iz=1:NZ
    Nc(iz,:) = histc(c(:,iz),clim);
    Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
    Nr(iz,:) = histc(r(:,iz),rlim);
    Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
    Na(iz,:) = histc(a(:,iz),alim);
    Na(iz,:) = Na(iz,:)/max(Na(iz,:));
end;
disp('Done histograms.');
end; % end imarg

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure(1);hold on; box on;
set(hf2, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box off;

subplot(h1)
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',1);
end;
if(imap == 1)
   plprof(map,hmax,'k',1);
end;
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[1500 1600 1700 1800]);
set(gca,'XTick',[1500 1600 1700 1800]);
box on;

subplot(h2)
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
   plprof(mtru,hmax,'--k',2);
end
if(imap == 1)
   plprof(map,hmax,'k',2);
end;
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',[1.4 1.6 1.8 2.0]);
box on;

subplot(h3)
if(imarg == 1)
   pcolor(alim,z,Na);shading flat;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)');
if(isyn == 1)
   plprof(mtrutmp,hmax,'--k',3);
end;
if(imap == 1)
   plprof(map,hmax,'k',3);
end;
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;

if(isrc == 1)
   subplot(h1);
   plot((source(1:100,2)*100)+1680,(1500*source(1:100,1)-.1),'k','Linewidth',1);
end
barstep1 = 10;
barstep1r = 1;
barstep2 = 10;
barstep2r = 1;
barstep3 = 10;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'g','Linewidth',1);
    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
               10.*ones(size(c1(1:barstep1:end,2))),...
               zeros(size(c1(1:barstep1:end,2))),'g','g');

    subplot(h2);
    plot(r1(:,2),r1(:,1),'g','Linewidth',1);
    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
               2./100.*r1(1:barstep1r:end,2),...
               zeros(size(r1(1:barstep1r:end,2))),'g','g');

end;

if(iplot == 1)
%    saveas(hf2,plotfile1,'epsc2');
    exportfig(hf2,plotfile1,opts);
end;

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');

return;
