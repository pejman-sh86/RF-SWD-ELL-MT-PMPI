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
ifull = 1;
inorm = 1;  % 1=line by line; 2=normalize by global max
isrc = 0;
Tstar = 1.;

%sample       = 'x_s20_1_8_80_5lay_sample.mat';
%mapfile      = 'x_s20_1_8_80_5lay_fgsmap.txt';
%profilefile  = 'x_s20_1_8_80_5lay_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s20_1_8_80_5lay_profileshade.png';
%else
%   plotfile1 = 'x_s20_1_8_80_5lay_profilemap.png';
%end
%corefile  = 'core.mat';

%sample      = 'x_s16_1_5_40_7lay_sample.mat';
%mapfile     = 'x_s16_1_5_40_7lay_fgsmap.txt';
%profilefile = 'x_s16_1_5_40_7lay_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s16_1_5_40_7lay_profileshade.eps';
%else
%   plotfile1 = 'x_s16_1_5_40_7lay_profilemap.eps';
%end
%corefile  = 'core15.mat';

sample      = 'x_s13_2_10_50_5layb_sample.mat';
mapfile     = 'x_s13_2_10_50_5layb_fgsmap.txt';
profilefile = 'x_s13_2_10_50_5layb_profile.mat';
if(imarg == 1)
   plotfile1 = 'x_s13_2_10_50_5layb_profileshade.eps';
else
   plotfile1 = 'x_s13_2_10_50_5layb_profilemap.eps';
end
corefile  = 'core13.mat';

%sample      = 'x_s07_1_1_100_2lay_T2_sample.mat';
%mapfile     = 'x_s07_1_1_100_2lay_T2_fgsmap.txt';
%profilefile = 'x_s07_1_1_100_2lay_T2_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s07_1_1_100_2lay_T2_profileshade.png';
%else
%   plotfile1 = 'x_s07_1_1_100_2lay_T2_profilemap.eps';
%end
%corefile  = 'core.mat';

%sample      = 'x_s05_2_8_25_2lay_sample.mat';
%mapfile     = 'x_s05_2_8_25_2lay_fgsmap.txt';
%profilefile = 'x_s05_2_8_25_2lay_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s05_2_8_25_2lay_profileshade.png';
%else
%   plotfile1 = 'x_s05_2_8_25_2lay_profilemap.eps';
%end
%corefile  = 'core.mat';
%nffile = 'nf_band3.mat';

%sample      = 'x_s02_1_3_25_6lay_sample.mat';
%mapfile     = 'x_s02_1_3_25_6lay_fgsmap.txt';
%profilefile = 'x_s02_1_3_25_6lay_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s02_1_3_25_6lay_profileshade.eps';
%else
%   plotfile1 = 'x_s02_1_3_25_6lay_profilemap.eps';
%end
%corefile  = 'core.mat';

%sample      = 'x_s01_1_3_20_5lay_sample.mat';
%mapfile     = 'x_s01_1_3_20_5lay_fgsmap.txt';
%profilefile = 'x_s01_1_5_20_5lay_profile.mat';
%if(imarg == 1)
%   plotfile1 = 'x_s01_1_3_20_5lay_profileshade.eps';
%else
%   plotfile1 = 'x_s01_1_3_20_5lay_profilemap.eps';
%end
%corefile  = 'core.mat';

if(isrc == 1)
   sourcefile = '../source_13_low.dat';
   source = load(sourcefile);
   source(:,1) = source(:,1) - source(1,1);
end;

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
%
% Sim C
%
%   mtru = [0.31, 1480, 1.30, 1.6, 0.01, ...
%           0.35, 1510, 1.65, 0.30, ...
%           1.20, 1500, 1.55, 0.30, ...
%           0.40, 1700, 1.90, 0.30, ...
%                 1560, 1.70, 0.01 ];
%   mtrutmp = [mtru(1:3) mtru(5:end)];
end;
nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

NPARL = 4;
hmax = 5.2;

pmin = [1450 1.35 0.0]';
pmax = [1650 2.10 1.0]';

%
% Site 01
%
%pmin = [1450 1.2 0.0];
%pmax = [1900 2.4 1.0];

%
% # bins in z and parameter direction 
%
if(imarg == 1)
   NZ = 500;
   NC = 500;
   NR = 400;
   NA = 300;
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
   E = A(:,1);
   NSAMP = length(m);
%
% Set # profiles NPROF 
%
   if(ifull == 0)
       NPROF = 60000;
       ranuni = rand(NPROF,1);
       m = m(round(NSAMP*ranuni),:);
       E = E(round(NSAMP*ranuni));
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
c = zeros(NPROF,NZ);
r = zeros(NPROF,NZ);
a = zeros(NPROF,NZ);
for iprof = 1:NPROF
    
    if(rem(iprof,5000)==0)
        fprintf(1,'%8i',iprof)
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
Nc = zeros(NZ,NC);
Nr = zeros(NZ,NR);
Na = zeros(NZ,NA);
%if(Tstar == 1.)
%   for iz=1:NZ
%      Nc(iz,:) = histc(c(:,iz),clim);
%      Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
%      Nr(iz,:) = histc(r(:,iz),rlim);
%      Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
%      Na(iz,:) = histc(a(:,iz),alim);
%      Na(iz,:) = Na(iz,:)/max(Na(iz,:));
%   end;
%else
for iz=1:NZ
%   Nc(iz,:) = histc_tstar(c(:,iz),E,clim,Tstar);
%   Nr(iz,:) = histc_tstar(r(:,iz),E,rlim,Tstar);
%   Na(iz,:) = histc_tstar(a(:,iz),E,alim,Tstar);
   Nc(iz,:) = histc(c(:,iz),clim);
   Nr(iz,:) = histc(r(:,iz),rlim);
   Na(iz,:) = histc(a(:,iz),alim);
end;

%clim = [pmin(1)*ones(NZ,1),clim,pmax(1)*ones(NZ,1)];
%rlim = [pmin(2)*ones(NZ,1),rlim,pmax(2)*ones(NZ,1)];
%alim = [pmin(3)*ones(NZ,1),alim,pmax(3)*ones(NZ,1)];
%Nc   = [zeros(NZ,1),Nc,zeros(NZ,1)];
%Nr   = [zeros(NZ,1),Nr,zeros(NZ,1)];
%Na   = [zeros(NZ,1),Na,zeros(NZ,1)];

%
% Normalize Histograms
%
disp('Normalize histograms...');
if(inorm == 1)
   for iz=1:NZ
      Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
      Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
      Na(iz,:) = Na(iz,:)/max(Na(iz,:));
      [nfc(iz,:)] = hpd(c(:,iz),100,95);
      [nfr(iz,:)] = hpd(r(:,iz),100,95);
      [nfa(iz,:)] = hpd(a(:,iz),100,95);
      meac(iz) = median(c(:,iz));
      mear(iz) = median(r(:,iz));
      meaa(iz) = median(a(:,iz));
   end;
elseif(inorm == 2)
   Nc = Nc/max(max(Nc));
   Nr = Nr/max(max(Nr));
   Na = Na/max(max(Na));
end;
save tmp.mat z nfc nfa nfr meac mear meaa;
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
   size(clim)
   size(z)
   size(Nc)
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
    plot(c1(:,2),c1(:,1),'k','Linewidth',1);
    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
               10.*ones(size(c1(1:barstep1:end,2))),...
               zeros(size(c1(1:barstep1:end,2))),'k','k');
%   plot(c2(:,2),c2(:,1),'g','Linewidth',1);
%    errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%               10.*ones(size(c2(1:barstep2:end,2))),...
%               zeros(size(c2(1:barstep2:end,2))),'g','g');
%    plot(c3(:,2),c3(:,1),'r','Linewidth',1);
%    errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%               10.*ones(size(c3(1:barstep3:end,2))),...
%               zeros(size(c3(1:barstep3:end,2))),'r','r');


    subplot(h2);
    plot(r1(:,2),r1(:,1),'k','Linewidth',1);
    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
               2./100.*r1(1:barstep1r:end,2),...
               zeros(size(r1(1:barstep1r:end,2))),'k','k');
%   plot(r2(:,2),r2(:,1),'g','Linewidth',1);
%   errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%              2./100.*r2(1:barstep2r:end,2),...
%              zeros(size(r2(1:barstep2r:end,2))),'g','g');
%    plot(r3(:,2),r3(:,1),'r','Linewidth',1);
%    errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%               2./100.*r3(1:barstep3:end,2),...
%               zeros(size(r3(1:barstep3:end,2))),'r','r');


end;

if(iplot == 1)
    saveas(hf2,plotfile1,'png');
%    exportfig(hf2,plotfile1,opts);
end;

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');

return;
