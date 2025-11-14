load hists.mat

NX = size(sl1_ens,2);
NY = size(sl1_ens,1);
over = 5; % amount of overlap (grid points)
nx = 15;
ny = 15;
nx2 = nx + 2*over;
ny2 = ny + 2*over;
nx3 = NX*nx + 2*over;
ny3 = NY*ny + 2*over;

l  = 1.0; % fine grid spacing (km)
L  = 15; % coarse grid spacing between source patches (grid points)

ltape = over*2;     % taper length (should be two times the overlap)

X = zeros(ny3,nx3); % total source region
ss = ones(ny2,nx2); % source array

%% Now design a cos taper:
tx=[-pi/2:pi/(ltape-1):pi/2];
ty=[0:pi/(ltape-1):pi];
ta = (sin(tx)+1)/2;
tb = (cos(ty)+1)/2;
tape = ones(nx2,1);
tape(1:ltape) = ta;
tape(end-ltape+1:end) = tb;

%% Apply taper to source array ss:
for i=1:nx2;
  ss(:,i) = ss(:,i).*tape;
  ss(i,:) = ss(i,:).*tape';
end;
figure,surf(ss);

load hists;
sl1_ens = sum(sl1_ens,3);
sl1_ens_sm = zeros(ny3,nx3);
j1 = 1;
for j2=1:NY;
  i1 = 1;
  for i2=1:NX;
    sl1_ens_sm(j1:j1+ny2-1,i1:i1+nx2-1) = sl1_ens_sm(j1:j1+ny2-1,i1:i1+nx2-1) + sl1_ens(j2,i2)*ss;

    i1 = i1 + nx;
  end;
  j1 = j1 + ny;
end;

figure();
contourf(sl1_ens_sm,20)
colorbar;colormap(darkb2r(-3,10));
set(gca,'YDir','reverse');axis equal;
set(gca,'XLim',[0 nx3],'YLim',[0 ny3])

figure();
imagesc(sl1_ens_sm)
colorbar;colormap(darkb2r(-3,10));
set(gca,'YDir','reverse');axis equal;
set(gca,'XLim',[0 nx3],'YLim',[0 ny3])

nx = 8;
ny = 2;
xim = 0.005/nx;
yim = 0.01/ny;
xymarg = [0.02 0.02 0.02 0.02];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
figure();
for i=1:16;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  imagesc(bins(1,:),[1:32],squeeze(sl1_tot_hst2(:,i,:)));
  axis equal;set(gca,'YDir','reverse','YLim',[0 33],'XLim',[-3 10]);
end;

figure();
for i=1:16;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  imagesc(bins(1,:),[1:32],squeeze(sl1_hst2(:,i,:,1)));
  axis equal;set(gca,'YDir','reverse','YLim',[0 33],'XLim',[-3 10]);
end;

figure();
for i=1:16;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  imagesc(bins(1,:),[1:32],squeeze(sl1_hst2(:,i,:,2)));
  axis equal;set(gca,'YDir','reverse','YLim',[0 33],'XLim',[-3 10]);
end;

figure();
for i=1:16;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  imagesc(bins(1,:),[1:32],squeeze(sl1_hst2(:,i,:,3)));
  axis equal;set(gca,'YDir','reverse','YLim',[0 33],'XLim',[-3 10]);
end;

figure();
for i=1:16;
  subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
  imagesc(bins(1,:),[1:32],squeeze(sl1_hst2(:,i,:,4)));
  axis equal;set(gca,'YDir','reverse','YLim',[0 33],'XLim',[-3 10]);
end;


