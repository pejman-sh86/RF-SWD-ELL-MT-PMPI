function [] = plot_sinrg_core(samplefile);

set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);
set(gcf, 'renderer', 'painters')

icore = 1;
ifull = 1;
inorm = 1;
imap  = 0;
iplot = 1;
iscratch = 1;
IBFIX = 0;
NLGRAD = 100;
hmax = 2.0;
pmin = [1450 1.1  0.0];
pmax = [1500 1.75 1.0];

filebase = strrep(samplefile,'sample.mat','');
%%
%%  Files:
%%
depthmarg   = strcat(filebase,'depthmarg.mat');
mapfile     = strcat(filebase,'map.dat');
plotfile1   = strcat(filebase,'profileshade.png');
profilefile = strcat(filebase,'profile.mat');
corefile    = 'core.mat';

NPARL = 4;
NZ = 300;
NC = 300;
NR = 300;
NA = 100;

NPAR = NPARL * NLGRAD + 3;
load(samplefile);
if(imap == 1)
   map = load(mapfile);
end;
if(icore == 1)
   load(corefile);
end;

if(imap == 1)
   for i = 1:NLGRAD
      znorm(i) = ((1./NLGRAD)-(1./NLGRAD/2.))+((1./NLGRAD)*(i-1));
      m_rg(((i-1)*4)+1) = map(1)/NLGRAD;
      m_rg(((i-1)*4)+2) = map(2)+(map(3)-map(2))/(NLGRAD-1)*(i-1);
      m_rg(((i-1)*4)+3) = map(4)+sin(znorm(i)*pi/2.)^map(6)...
                          *(map(5)-map(4));
      m_rg(((i-1)*4)+4) = map(end);
   end
   if(IBFIX == 0)
      m_rg(end-2:end) = [map(end-4),map(end-2),map(end)];
   else
      m_rg(end-2:end) = [map(end-2:end)];
   end
end;

m = A(:,5:11);
NSAMP = length(m);

clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

dz = hmax/(NZ-1);
z = cumsum(dz*ones(1,NZ))-dz;

if(ifull == 0)
    NPROF = 1000;
    m = m(round(NSAMP*rand(NPROF,1)),:);
else
    NPROF = NSAMP;
end
disp('Make gradient sample...')
m2 = zeros(NPROF,NLGRAD*4);
for j = 1:NPROF
   for i = 1:NLGRAD
      znorm(i) = ((1./NLGRAD)-(1./NLGRAD/2.))+((1./NLGRAD)*(i-1));
      m2(j,((i-1)*4)+1) = m(j,1)/NLGRAD;
      m2(j,((i-1)*4)+2) = m(j,2)+(m(j,3)-m(j,2))/(NLGRAD-1)*(i-1);
      m2(j,((i-1)*4)+3) = m(j,4)+sin(znorm(i)*pi/2.)^m(j,6)...
                          *(m(j,5)-m(j,4));
      m2(j,((i-1)*4)+4) = m(j,end);
   end
   if(IBFIX == 0)
      m2(j,NPAR-2:NPAR) = [m(j,end-4),m(j,end-2),m(j,end)];
   else
      m2(j,NPAR-2:NPAR) = [m(j,end-2:end)];
   end
end;
disp('...done gradient sample.')

NLAY = NLGRAD;

idxh = (([1:NLAY]-1)*NPARL)+1;
idxc = (([1:NLAY]-1)*NPARL)+2;
idxr = (([1:NLAY]-1)*NPARL)+3;
idxa = (([1:NLAY]-1)*NPARL)+4;
idxh = [idxh idxh(end)];
idxc = [idxc idxc(end)+3];
idxr = [idxr idxr(end)+3];
idxa = [idxa idxa(end)+3];

disp('Starting profiles...');

prof(1:NPROF,:,1) = cumsum(m2(1:NPROF,idxh),2);
prof(1:NPROF,:,2) = m2(1:NPROF,idxc);
prof(1:NPROF,:,3) = m2(1:NPROF,idxr);
prof(1:NPROF,:,4) = m2(1:NPROF,idxa);

disp('Sample size: '),disp(size(m2))
idx=1;
c = zeros(NPROF,NZ);
r = zeros(NPROF,NZ);
a = zeros(NPROF,NZ);
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

save(depthmarg,'c','r','a');

disp('Starting histograms...');
for iz=1:NZ
  [Nc(iz,:),binsc] = hist(c(:,iz),clim);
  [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
  [Na(iz,:),binsa] = hist(a(:,iz),alim);
  if(inorm == 0);
    Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
    Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
    Na(iz,:) = Na(iz,:)/max(Na(iz,:));
  else;
    Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
    Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
    Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
  end;
end;
disp('Done histograms.');

save tmp.mat

NLAY  = NLGRAD;

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure(1);hold on; box on;
set(hf2, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box off;

subplot(h1);
pcolor(clim,z,Nc);shading flat;
set(h1,'layer','top')
if(imap == 1)
   plprof(m_rg,hmax,'k',1);
end;
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Velocity (m/s)','Fontsize',18);
ylabel('Depth (m)','Fontsize',18);
set(gca,'Fontsize',16,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[1450 1500 1550]);
set(gca,'XTick',[1450 1500 1550]);
box on;

subplot(h2);
pcolor(rlim,z,Nr);shading flat;
set(h1,'layer','top')
if(imap == 1)
   plprof(m_rg,hmax,'k',2);
end;
set(h2,'layer','top')
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)','Fontsize',18);
set(gca,'Fontsize',16,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6]);
set(gca,'XTick',[1.4 1.6]);
box on;

subplot(h3);
pcolor(alim,z,Na);shading flat;
set(h1,'layer','top')
if(imap == 1)
   plprof(m_rg,hmax,'k',3);
end;
set(h3,'layer','top')
set(gca,'YDir','reverse');
xlabel('Attenuation (dB/L)','Fontsize',18);
set(gca,'Fontsize',16,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YTickLabel',[]);
cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;

barstep1 = 10;
barstep1r = 10;
barstep2 = 10;
barstep2r = 10
barstep3 = 10;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'w','Linewidth',3);
    plot(c1(:,2),c1(:,1),'k','Linewidth',1);
    errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
               10.*ones(size(c1(1:barstep1:end,2))),...
               zeros(size(c1(1:barstep1:end,2))),'k','k');
%    plot(c2(:,2),c2(:,1),'k','Linewidth',1);
%    errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%               2*10.*ones(size(c2(1:barstep2:end,2))),...
%               zeros(size(c2(1:barstep2:end,2))),'k','k');
%    plot(c3(:,2),c3(:,1),'r','Linewidth',1);
%    errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%               10.*ones(size(c3(1:barstep3:end,2))),...
%               zeros(size(c3(1:barstep3:end,2))),'r','r');


    subplot(h2);
    plot(r1(:,2),r1(:,1),'w','Linewidth',3);
    plot(r1(:,2),r1(:,1),'k','Linewidth',1);
    errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
               2./100.*r1(1:barstep1r:end,2),...
               zeros(size(r1(1:barstep1r:end,2))),'k','k');
%    plot(r2(:,2),r2(:,1),'k','Linewidth',1);
%    errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%               2./100.*r2(1:barstep2r:end,2),...
%               zeros(size(r2(1:barstep2r:end,2))),'k','k');
%    plot(r3(:,2),r3(:,1),'r','Linewidth',1);
%    errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%               2./100.*r3(1:barstep3:end,2),...
%               zeros(size(r3(1:barstep3:end,2))),'r','r');
%

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
