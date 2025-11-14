function [] = sim_longtrack();

%fr=[1000,1200,2000,2400,2800,3200];
fr=[1000,1200,2000,2400];
angobs = [2.91363560E+01  3.03432690E+01  3.15501830E+01  3.27570960E+01  3.39640100E+01  3.51709230E+01  3.63778370E+01  3.75847510E+01  3.87916640E+01  3.99985780E+01  4.12054910E+01  4.24124050E+01  4.36193180E+01  4.48262320E+01  4.60331450E+01  4.72400590E+01  4.84469720E+01  4.96538860E+01  5.08608000E+01  5.20677130E+01  5.32746270E+01  5.44815400E+01  5.56884540E+01  5.68953670E+01  5.81022810E+01  5.93091940E+01  6.05161080E+01  6.17230220E+01  6.29299350E+01  6.41368490E+01  6.53437620E+01  6.65506760E+01];

track = dlmread('track_environment.dat');
zt = 2.3930223e+01;
cw = 1.5120000e+03;
rw = 1.0290000e+00;
zmx= 4.0000000e+00;
NZ = 400;
dz = zmx/(NZ-1);
z = cumsum(dz*ones(1,NZ))-dz;
h = zeros(1,NZ);
NPROF = size(track,1);
NPARL = 4;
c = zeros(NPROF,NZ);
r = zeros(NPROF,NZ);
a = zeros(NPROF,NZ);

fig1=figure();hold on;box on;
fig2=figure();hold on;box on;
set(gca,'YDir','reverse');
fig3=figure();hold on;box on;
set(gca,'YDir','reverse');
offset1 = 0;
offset2 = 0;
for iprof = 1:size(track,1);
   clear m;
   filename = strcat('p',num2str(iprof,'%03i'),'_1000_2400.txt');
   plotname = strcat('p',num2str(iprof,'%03i'),'_1000_2400.png');
   k = track(iprof,1);
   NFP = k*4+3;
   sd = track(iprof,NFP+2)
   m = track(iprof,2:NFP+1);
   [ref]=forward_ref_nlay3(m,angobs,fr);
   ref = ref+randn(size(ref))*sd;
   ref = ref';
   save(filename,'zt','-ascii');
   save(filename,'cw','-ascii','-append');
   save(filename,'rw','-ascii','-append');
   save(filename,'zmx','-ascii','-append');
   save(filename,'ref','-ascii','-append');
   save(filename,'angobs','-ascii','-append');
   figure(fig1);
   for j = 1:length(fr);
%      subplot(2,3,j);hold on;box on;
      subplot(2,3,j);box on;
%      plot(angobs,ref(j,:),'-b');
%      set(gca,'XLim',[angobs(1) angobs(end)],'YLim',[0 .65]);
%      if(j==1 | j==4);ylabel('Refl. coeff.');end;
%      if(j > 3);xlabel('Angle (deg.)');end;
   end;
%   saveas(fig1,plotname,'png');
%   figure(fig2);
%   offset1 = offset1 + 120;
%   plprof(m,4,'-k',1,offset1);
%   figure(fig3);
%   offset2 = offset2 + .25;
%   plprof(m,4,'-k',2,offset2);

   %%
   %% Plot pcolor track profile
   %%
   clear idxh idxc idxr idxa prof;
   %% Find index for current model
   if(k > 0)
      idxh = (([1:k]-1)*NPARL)+1;
      idxc = (([1:k]-1)*NPARL)+2;
      idxr = (([1:k]-1)*NPARL)+3;
      idxa = (([1:k]-1)*NPARL)+4;
      idxh = [idxh idxh(end)];
      idxc = [idxc idxc(end)+3];
      idxr = [idxr idxr(end)+3];
      idxa = [idxa idxa(end)+3];
   else
      idxh = [];
      idxc = [1];
      idxr = [2];
      idxa = [3];
   end

   %% Compute the profile for current model
   if(k > 0)
      prof(1:k,1) = cumsum(m(idxh(1:end-1)),2);
      prof(k+1,1) = prof(k,1)+m(idxh(end));
   else
      prof(1,1) = zmx;
   end
   prof(:,2) = m(idxc);
   prof(:,3) = m(idxr);
   prof(:,4) = m(idxa);

   c(iprof,:) = prof(1,2);
   r(iprof,:) = prof(1,3);
   a(iprof,:) = prof(1,4);
   for ilay=2:k+1  %% k is # layers of current model
      idxz = round(prof(ilay-1,1)/dz);
%      h(idxz)     = h(idxz) + 1;
      c(iprof,idxz:end) = prof(ilay,2);
      r(iprof,idxz:end) = prof(ilay,3);
      a(iprof,idxz:end) = prof(ilay,4);
   end;
end;

figw = 10;
figh = 6.6667;
fig4=figure();
set(fig4,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
nx = 1;
ny = 3;
xim = 0.01;
yim = 0.02;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
pcolor([1:NPROF],z,c');shading flat;
set(gca,'XTickLabel',[])
ylabel('Depth (m)');
set(gca,'Layer','top','YDir','reverse');box on;
set(gca,'XLim',[1 NPROF],'YLim',[0 zmx],'CLim',[1450 1750]);box on;
c1=colorbar;ylabel(c1,'Velocity (m/s)')

subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
pcolor([1:NPROF],z,r');shading flat;
set(gca,'XTickLabel',[])
ylabel('Depth (m)');
set(gca,'Layer','top','YDir','reverse');box on;
set(gca,'XLim',[1 NPROF],'YLim',[0 zmx],'CLim',[1.2 2.2]);box on;
c2=colorbar;ylabel(c2,'Density (g/ccm)')

subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
pcolor([1:NPROF],z,a');shading flat;
xlabel('Ping no.');
ylabel('Depth (m)');
set(gca,'Layer','top','YDir','reverse');
set(gca,'XLim',[1 NPROF],'YLim',[0 zmx],'CLim',[0 1]);
box off;
box on;
c3=colorbar;ylabel(c3,'Attenuation (dB/m/kHz)')
%print(fig4,'-painters','-r300','track_true_model.png','-dpng');
print(fig4,'-painters','-r300','track_true_model.eps','-depsc');
return;
