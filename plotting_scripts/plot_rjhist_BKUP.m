function [] = plot_rjhist(filename);

imarg = 1; %% Plot depth-marginal distributions?
inorm = 1; %% Normalize profile marginals line by line
hmax = 4.0;
Tstar = 1.;
isyn = 1;
imap = 0;
imead = 1;
imean = 0;
isave = 1;

filebase    = strrep(filename,'sample.mat','')
mapfile     = strcat(filebase,'fgsmap.txt')
profilefile = strcat(filebase,'profile.mat')
plotfile1   = strcat(filebase,'profileshade.png')
corefile    = 'core.mat';

load(filename);

NPROF = 1e3;
%NPROF = length(A);
idxm = randperm(length(A));

if(isyn == 1)
   mtru = [0.60, 1480, 1.30, 0.10, ...
           1.00, 1550, 1.60, 0.30, ...
           0.60, 1585, 1.78, 0.30, ...
                 1600, 1.80, 0.40 ];
   mtrutmp = mtru;
end;

k = A(:,3);
NFP = (k*4)+3;
m = A(:,4:end-5);
logL = A(:,1);

for i = 1:size(A,1);
   sd(i) = m(i,NFP(i)+1);
   idxh = [0:k(i)-1]*4+1;
   h(i) = sum(m(i,idxh));
   clear idxh;
end

nx = 3;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

figure;
subplot(2,1,1);
plot(h)
subplot(2,1,2);
plot(sd);

figure;
subplot(2,1,1);
plot(NFP)
subplot(2,1,2);
hist(k,[1:10]);
set(gca,'XLim',[0.5 10.5]);

figure;
hist(logL,100);

%%
%% PLOT PROFILE MARGINALS
%%
NPARL = 4;
pmin = [1450 1.2 0.0]';
pmax = [1700 2.0 1.5]';
if(imarg == 1)
   NZ = 200;
   NC = 200;
   NR = 200;
   NA = 200;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

   dz = hmax/(NZ-1);
   z = cumsum(dz*ones(1,NZ))-dz;

   disp('Sample size: '),disp(size(m))
   m = m(idxm(1:NPROF),:);
   disp('Plotting size: '),disp(size(m))
   for iprof = 1:NPROF

      if(rem(iprof,1000)==0)
          fprintf(1,'%8i',iprof)
     end
     clear idxh idxc idxr idxa prof;
      %% Find index for current model
      idxh = (([1:k(iprof)]-1)*NPARL)+1;
      idxc = (([1:k(iprof)]-1)*NPARL)+2;
      idxr = (([1:k(iprof)]-1)*NPARL)+3;
      idxa = (([1:k(iprof)]-1)*NPARL)+4;
      idxh = [idxh idxh(end)];
      idxc = [idxc idxc(end)+3];
      idxr = [idxr idxr(end)+3];
      idxa = [idxa idxa(end)+3];

      %% Compute the profile for current model
      prof(:,1) = cumsum(m(iprof,idxh),2);
      prof(:,2) = m(iprof,idxc);
      prof(:,3) = m(iprof,idxr);
      prof(:,4) = m(iprof,idxa);

      izold = 1;
      for ilay=1:k(iprof)  %% k is # layers of current model
          idxz = find(z <= prof(ilay,1));
          c(iprof,idxz(izold:end)) = prof(ilay,2);
          r(iprof,idxz(izold:end)) = prof(ilay,3);
          a(iprof,idxz(izold:end)) = prof(ilay,4);
          izold = idxz(end)+1;
      end;
      c(iprof,izold:NZ) = prof(end,2);
      r(iprof,izold:NZ) = prof(end,3);
      a(iprof,izold:NZ) = prof(end,4);
   end;
   fprintf(1,'\n')
   disp('Done with profiles.');

   %
   % Compute histograms for each depth
   %
   disp('Starting histograms...');
   for iz=1:NZ
      Nc(iz,:) = histc_tstar(c(:,iz),logL,clim,Tstar);
      Nr(iz,:) = histc_tstar(r(:,iz),logL,rlim,Tstar);
      Na(iz,:) = histc_tstar(a(:,iz),logL,alim,Tstar);
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         Na(iz,:) = Na(iz,:)/max(Na(iz,:));
      end;
   elseif(inorm == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      Na = Na/max(max(Na));
   end;
   disp('Done histograms.');
end; % end imarg

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
c_mead = median(c);
r_mead = median(r);
a_mead = median(a);

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure;hold on; box on;
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
   if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
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
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
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
   if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
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

if(isave == 1)
   saveas(hf2,plotfile1,'png');
%   exportfig(hf2,plotfile1,opts);
   save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');
end;

return;
