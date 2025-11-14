%
% Plot 2-D profile density
%
function plot_profile_disp;

set(0, 'DefaultFigurePaperPosition', [0 0 5 8]);
set(gcf, 'renderer', 'painters')

i_smp   = 0;   % 0=MH; 1=nested; 2=ais
i_grad  = 0;   % 0=; 1=linear; 2=power
ilpl    = 1;   % Layer that has power law/linear gradient
isd     = 1;   % are std devs sampled over?
NSD     = 1;   % how many std devs were sampled over?
iosler= 1;
imarg = 1;
iplot = 1;
ifull = 1;
inorm = 1;  % 1=line by line; 2=normalize by global max
Tstar = 1.;
NLGRAD = 20;
NPARL = 4;
hmax = 60.0;
NPROF = 100000;

sample      = 'mh_3lay_sample.mat';
profilefile = 'mh_3lay_profile.mat';
if(imarg == 1)
%   plotfile1 = 'mh_lg_profileshade.png';
   plotfile1 = 'mh_3lay_profileshade.eps';
else
   plotfile1 = 'mh_3lay_profilemap.eps';
end
oslerfile1  = 'eb_osler_curve.dat';

nx = 1;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',4.6,'Height',7,'Color','rgb',...
              'Renderer','painters','Format','eps',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

pmin = [0   1450 1.3]';
pmax = [210 2200 2.2]';

%
% # bins in z and parameter direction 
%
if(imarg == 1)
   NZ = 1000;
   NC = 800;
   NR = 50;
   NA = 50;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));

   dz = hmax/(NZ-1);
   z = cumsum(dz*ones(1,NZ))-dz;
end

if(iosler == 1)
    osler1 = load(oslerfile1);
end;

if(imarg == 1)
   load(sample);
   if(i_smp == 0)
      istart = 2;
      E = A(:,1);
   elseif(i_smp == 1)
      istart = 6;
      logL     = A(:,1);
      logwidth = A(:,2);
      logwt    = A(:,3);
      H1        = A(:,4);
      logZ1     = A(:,5);
      w = exp(logwt-logZ1(end));
      fprintf(1,'Nested Sampling Evidence estimate: Znst = %12.4f \n',logZ1(end));
   elseif(i_smp == 2)
      istart = 3;
      logwt3 = A(:,1);
      logL3  = A(:,2);
   end;
   if(isd == 0)
      m = A(:,istart:end);
   else
      m = A(:,istart:end-NSD);
   end;
   NFP = length(m(1,:));
   NSAMP = length(m);
%
% Set # profiles NPROF 
%
   if(ifull == 0)
      ranuni = rand(NPROF,1);
      m = m(round(NSAMP*ranuni),:);
      if(i_smp == 0)
         E = E(round(NSAMP*ranuni));
      elseif(i_smp == 1)
         w = w(round(NSAMP*ranuni));
      elseif(i_smp == 2)
         logwt3 = logwt3(round(NSAMP*ranuni));
      end;
   else
       NPROF = NSAMP;
   end;

if(i_grad == 0)
   NLAY  = (size(m,2)-3)/4
   NLAY2 = NLAY;
   idxh  = (([1:NLAY]-1)*NPARL)+1;
   idxvs = (([1:NLAY]-1)*NPARL)+2;
   idxvp = (([1:NLAY]-1)*NPARL)+3;
   idxr  = (([1:NLAY]-1)*NPARL)+4;
elseif(i_grad == 1)
   %% LINEAR GRADIENT
   NLAY  = (size(m,2)-4)/4
   NLAY2 = NLAY - 1 + NLGRAD
   j = 1;
   for i=1:NLAY
      idxh(i)  = j;j=j+1;
      if(i == ilpl);
         idxvs1 = j;j=j+1;
%         idxvs2 = j;j=j+1;
      end;
      idxvs(i) = j;j=j+1;
      idxvp(i) = j;j=j+1;
      idxr(i)  = j;j=j+1;
   end
elseif(i_grad == 2)
   %% POWER LAW
   NLAY  = (size(m,2)-5)/4
   NLAY2 = NLAY - 1 + NLGRAD
   j = 1;
   for i=1:NLAY
      idxh(i)  = j;j=j+1;
      if(i == ilpl);
         idxvs1 = j;j=j+1;
         idxvs2 = j;j=j+1;
      end;
      idxvs(i) = j;j=j+1;
      idxvp(i) = j;j=j+1;
      idxr(i)  = j;j=j+1;
   end
end
idxh  = [idxh  idxh(end)];
idxvs = [idxvs idxvs(end)+3];
idxvp = [idxvp idxvp(end)+3];
idxr  = [idxr  idxr(end)+3];

disp('Starting profiles...');

%
% Sort model paramerters into profiles of absolute depth
%
%mtmp = repmat([120., 100., 200., 2000., 1.5, 250., 2250., 2.5],length(m),1);
%m = mtmp+0.1*mtmp.*randn(size(mtmp));
%mtmp = repmat([120., 100.,110., 200., 2000., 1.5, 250., 2250., 2.5],length(m),1);
%m = mtmp;

if(i_grad == 0)
   prof(1:NPROF,:,1) = cumsum(m(1:NPROF,idxh),2);
   prof(1:NPROF,:,2) = m(1:NPROF,idxvs);
%   prof(1:NPROF,:,3) = m(1:NPROF,idxvp);
%   prof(1:NPROF,:,4) = m(1:NPROF,idxr);
elseif(i_grad == 1)
   %% LINEAR
   il2 = 1;
   for il = 1:NLAY
      if(il == ilpl)
         zsum = 0.;
         for i = 1:NLGRAD+1
            h = m(1:NPROF,idxh(il))/NLGRAD;
            if(i == 1)
               if(ilpl == 1)
%                  prof(1:NPROF,il2,1) = h;
                  prof(1:NPROF,il2,1) = 0;
               else
                  prof(1:NPROF,il2,1) = prof(1:NPROF,il2-1,1);
               end;
               prof(1:NPROF,il2,2) = m(1:NPROF,idxvs1);
            else
               prof(1:NPROF,il2,1) = prof(1:NPROF,il2-1,1)+h;
               zsum = prof(1:NPROF,il2,1);
               prof(1:NPROF,il2,2) = m(1:NPROF,idxvs1)+...
                                    (m(1:NPROF,idxvs(il))-m(1:NPROF,idxvs1))./...
                                     m(1:NPROF,idxh(il)).*zsum;
            end
%            prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(il));
%            prof(1:NPROF,il2,4) = m(1:NPROF,idxr(il));
            il2 = il2+1;
         end;
      else
         prof(1:NPROF,il2,1) = prof(1:NPROF,il2-1,1)+m(1:NPROF,idxh(il));
         prof(1:NPROF,il2,2) = m(1:NPROF,idxvs(il));
%         prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(il));
%         prof(1:NPROF,il2,4) = m(1:NPROF,idxr(il));
         il2 = il2+1;
      end;
   end;
   prof(1:NPROF,il2,1) = hmax;
   prof(1:NPROF,il2,2) = m(1:NPROF,idxvs(end));
%   prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(end));
%   prof(1:NPROF,il2,4) = m(1:NPROF,idxr(end));
elseif(i_grad == 2)
   %% POWER
   il2 = 1;
   for il = 1:NLAY
      if(il == ilpl)
         v0   = m(1:NPROF,idxvs1);
         b    = m(1:NPROF,idxvs2)-v0;
         a    = log((m(1:NPROF,idxvs(il))-v0)./(m(1:NPROF,idxvs2)-v0))./log(m(1:NPROF,idxh(il)));
         for i = 1:NLGRAD+1
            if(i == 1)
               if(ilpl == 1)
                  prof(1:NPROF,il2,1) = 0.;
               else
                  prof(1:NPROF,il2,1) = m(1:NPROF,idxh(il-1));
               end;
                  prof(1:NPROF,il2,2) = m(1:NPROF,idxvs1);
            else
               h = m(1:NPROF,idxh(il))/NLGRAD;
               prof(1:NPROF,il2,1) = prof(1:NPROF,il2-1,1)+h;
               zsum = prof(1:NPROF,il2,1);
               prof(1:NPROF,il2,2)  = (v0.*h+b./(a+1.).*(zsum.^(a+1.)- ...
                                      (zsum-h).^(a+1.)))./h;
            end;
%            prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(il));
%            prof(1:NPROF,il2,4) = m(1:NPROF,idxr(il));
            il2 = il2+1;
         end;
      else
         prof(1:NPROF,il2,1) = prof(1:NPROF,il2-1,1)+m(1:NPROF,idxh(il));
         prof(1:NPROF,il2,2) = m(1:NPROF,idxvs(il));
%         prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(il));
%         prof(1:NPROF,il2,4) = m(1:NPROF,idxr(il));
         il2 = il2+1;
      end;
   end;
   prof(1:NPROF,il2,1) = hmax;
   prof(1:NPROF,il2,2) = m(1:NPROF,idxvs(end));
%   prof(1:NPROF,il2,3) = m(1:NPROF,idxvp(end));
%   prof(1:NPROF,il2,4) = m(1:NPROF,idxr(end));
   save tmp prof;
end;

%
% Compute discretized form of profile
%
disp('Sample size: '),disp(size(m))
for iprof = 1:NPROF
    
    if(rem(iprof,1000)==0);fprintf(1,'%8i',iprof);end;
    izold = 1;
    for ilay=1:NLAY2
        idxz = find(z <= prof(iprof,ilay,1));
        vs(iprof,idxz(izold:end)) = prof(iprof,ilay,2);
%        vp(iprof,idxz(izold:end)) = prof(iprof,ilay,3);
%        r(iprof,idxz(izold:end)) = prof(iprof,ilay,4);
        izold = idxz(end)+1;
    end;
    vs(iprof,izold:NZ) = prof(iprof,end,2);
%    vp(iprof,izold:NZ) = prof(iprof,end,3);
%    r(iprof,izold:NZ) = prof(iprof,end,4);
end;
fprintf(1,'\n')
disp('Done with profiles.');

%
% Compute histograms for each depth
%
disp('Starting histograms...');
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
   if(i_smp == 0)
      Nc(iz,:) = histc_tstar(vs(:,iz),E,clim,Tstar);
   elseif(i_smp == 1)
      Nc(iz,:) = histc_wt(vs(:,iz),w,clim);
   elseif(i_smp == 2)
      Nc(iz,:) = histc_ais(vs(:,iz),logwt3,clim);
   end
%   Nr(iz,:) = histc_tstar(vp(:,iz),E,rlim,Tstar);
%   Na(iz,:) = histc_tstar(r(:,iz),E,alim,Tstar);
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
if(inorm == 1)
   for iz=1:NZ
      Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
%      Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
%      Na(iz,:) = Na(iz,:)/max(Na(iz,:));
   end;
elseif(inorm == 2)
   Nc = Nc/max(max(Nc));
%   Nr = Nr/max(max(Nr));
%   Na = Na/max(max(Na));
end;
disp('Done histograms.');
end; % end imarg

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf2 = figure(1);hold on; box on;
set(hf2, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
%h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
%hold on; box off;
%h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
%hold on; box off;

subplot(h1)
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'YDir','reverse');
xlabel('Shear Vel. (m/s)');
ylabel('Depth (m)');
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'XTickLabel',[100 200 300 400]);
set(gca,'XTick',[100 200 300 400]);
box on;

%subplot(h2)
%if(imarg == 1)
%   pcolor(rlim,z,Nr);shading flat;
%end;
%%surf(rlim,z,Nr);shading flat;
%set(h2,'layer','top')
%set(gca,'YDir','reverse');
%xlabel('Comp. Vel. (m/s)');
%set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
%set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[1400 1500 1600 1700]);
%set(gca,'XTick',[1400 1500 1600 1700]);
%box on;
%
%subplot(h3)
%if(imarg == 1)
%   pcolor(alim,z,Na);shading flat;
%end;
%%surf(alim,z,Na);shading flat;
%set(h3,'layer','top')
%set(gca,'YDir','reverse');
%xlabel('Density (g/ccm)');
%set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
%set(gca,'YTickLabel',[]);

cmap = colormap(flipud(gray));
cmap = colormap(jet);
cmap(1,:) = [1 1 1];
colormap(cmap);
box on;

barstep1 = 10;
barstep1r = 1;
barstep2 = 10;
barstep2r = 1;
barstep3 = 10;
if(iosler == 1)
    subplot(h1);
    plot(osler1(:,1),osler1(:,2),'--k','Linewidth',2);
end;

if(iplot == 1)
%    saveas(hf2,plotfile1,'png');
    exportfig(hf2,plotfile1,opts);
end;

vs_mean = mean(vs);
save(profilefile,'z', 'vs', 'vs_mean');
%vp_mean = mean(vp);
%r_mean = mean(r);
%save(profilefile,'z', 'vs', 'vp', 'r','vs_mean','vp_mean','r_mean');

return;
