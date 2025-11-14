clear; close all;

%%
%% Set up grid
%%
filebase = 'tohoku_data';
datfile  = strcat(filebase,'.hdf5');
NRAN    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx     = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx     = h5readatt(datfile,'/Sensitivity_kernel','max_y');
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
deltrn = double(deltr)/double(Rmx);
deltzn = double(deltz)/double(Zmx);
x_evn = [deltrn/2.:deltrn:1.-deltrn/2.];
z_evn = [deltzn/2.:deltzn:1.-deltzn/2.];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';

kmax = 30;
NTW = 3;
load sample.mat

for i=1:length(buf);
  step(:,:,i)=buf(i).step_sz;
  accept(:,:,i)=buf(i).iaccept;
  BIC(i)=buf(i).BIC;
  k(i) = buf(i).k;
  [a,b] = max(BIC);
  sl1(:,:,:,i)=buf(i).sl1;
end;
ens_sl1=mean(sl1,4);

mapobj = buf(b);
kmap = buf(b).k;
save('mapstart.mat','mapobj');

figure;
subplot(2,2,1);
plot(BIC);
subplot(2,2,3);
hist(BIC,100);
subplot(2,2,2);
plot(k);
subplot(2,2,4);
hist(k,[0:kmax]);

f1=figure;
f2=figure;
f3=figure;
for i=1:NTW;
  figure(f1);
  subplot(1,3,i);hold on;box on;
  imagesc(x_ev,z_ev,buf(b).sl1(:,:,i));
  axis equal;set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx])
  plot(buf(b).voro(1:kmap,1),buf(b).voro(1:kmap,2),'.k','Markersize',5);
  [vx,vz]=voronoi(buf(b).voro(1:kmap,1),buf(b).voro(1:kmap,2));plot(vx,vz,'k','LineWidth',2);
  colormap(darkb2r(-7.5,7.5));
  set(gca,'YDir','reverse');

  figure(f2);
  subplot(1,3,i);hold on;box on;
  imagesc(x_ev,z_ev,ens_sl1(:,:,i));
  axis equal;set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx])
  plot(buf(b).voro(1:kmap,1),buf(b).voro(1:kmap,2),'.k','Markersize',5);
  [vx,vz]=voronoi(buf(b).voro(1:kmap,1),buf(b).voro(1:kmap,2));plot(vx,vz,'--k','LineWidth',1);
  colormap(darkb2r(-7.5,7.5));
  set(gca,'YDir','reverse');

  figure(f3);
  subplot(1,3,i)
  imagesc(x_ev,z_ev,sl1tru2(:,:,i));
  colormap(darkb2r(-7.5,7.5));
  axis equal;set(gca,'XLim',[0 Rmx],'YLim',[0 Zmx])
end;

%% Excluded grid points:
mat_excl2 = repmat(mat_excl,[NTW,1]);
NDAT = sum(sum(mat_excl2));
rex = mat_excl2(:);

figure;hold on;
plot(buf(1).obs);
rep = buf(b).sl1(:);
plot(rep.*rex,'--r')

