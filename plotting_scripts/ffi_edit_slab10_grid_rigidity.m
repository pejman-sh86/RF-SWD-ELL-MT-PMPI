fileid  = '24km_2km_mat';
fileext1 = '.txt';
fileext2 = '_new.txt';

gfile      = strcat('grid_',fileid,fileext1);
gfilenew   = strcat('grid_',fileid,fileext2);
rigfile    = 'rig_mat.dat';
rigfilenew = strcat('rig_',fileid,fileext2);
edgefilenew = strcat('edge_',fileid,fileext2);

g=load(gfile);


nsx = 12;
nsy = 28;

%% Cut out desired part of grid:
for i=1:nsx;
    g2(i,:,:)=g((i-1)*nsy+1:i*nsy,:);
end;

nsx1 = 1;
nsx2 = 9;
nsy1 =  1;
nsy2 = 25;
NX = nsx2-nsx1+1;
NY = nsy2-nsy1+1;

N = nsx2-nsx1+1;
ii = 1;
for i=nsy1:nsy2;
    g3((ii-1)*N+1:ii*N,:)=squeeze(g2(nsx1:nsx2,i,:));
    ii = ii+1;
end;
plot(g(:,1),g(:,2));
hold on;
plot(g3(:,1),g3(:,2));
set(gca,'XLim',[140,145],'YLim',[35,41]);
save(gfilenew,'g3','-ascii');

re = [g3(1:NX:end,1),g3(1:NX:end,2)];
le = [g3(NX:NX:end,1),g3(NX:NX:end,2)];
plot(re(:,1),re(:,2),'ok')
plot(le(:,1),le(:,2),'ok')
edges = [le,re];
save(edgefilenew,'edges','-ascii');

%% Make rigidity file for each depth


z = g3(:,3);
II=[1:length(z)];
[zs,I]=sort(z);
[tmp,Iu]=sort(II(I));

rig=dlmread(rigfile);
rig(6,:)=[200,66.825];
rig2 = interp1(rig(:,1),rig(:,2),zs,'next');


rig3=rig2(Iu);
rig3=[z,rig2(Iu)];
figure();hold on;
stairs(rig(:,2),rig(:,1));
plot(rig2(Iu),z,'.');
save(rigfilenew,'rig3','-ascii');

%% After saving, paste first line to be:
% "#Rigidities (For tohoku eq interpolated to Slab1.0 grid):"

%clear;