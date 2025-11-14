function [] = compare_marg();

set(0, 'DefaultFigurePaperPosition', [0 0 8 4]);
set(gcf, 'renderer', 'painters')
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');


%%
%% Site 19
%%
profilefile19 = 'site19/packet_2m_new/exp_std/x_s19_1_8_25_0lay_profile.mat';
load(profilefile19);
c19 = c;
r19 = r;
a19 = a;

%%
%% Site 04
%%
profilefile4 = 'site04/holland_etal05/x_s04_1_3_16_0lay_profile.mat';
load(profilefile4);
c4 = c;
r4 = r;
a4 = a;

%%
%% Site 05
%%
profilefile5 = 'site05/packet_6m/x_s05_2_8_25_2lay_profile.mat';
load(profilefile5);
c5 = c;
r5 = r;
a5 = a;

%%
%% Site 01
%%
profilefile1 = 'site01/packet_3m/x_s01_1_3_20_5lay_profile.mat';
load(profilefile1);
c1 = c;
r1 = r;
a1 = a;

%%
%% Site 21
%%
profilefile21 = 'site21/packet_4m_new/tmp_6lay_20ave/x_s21_1_5_40_6lay_profile.mat';
load(profilefile21);
c21 = c;
r21 = r;
a21 = a;

[M19,NZ] = size(c19);
[M4,NZ]  = size(c4);
[M5,NZ]  = size(c5);
[M1,NZ]  = size(c1);
[M21,NZ] = size(c21);

NC1 = 800;
NC2 = 300;
NC3 = 200;
NC4 = 200;

NR1 = 300;
NR2 = 300;
NR3 = 300;
NR4 = 300;

pmin = [1450 1.2 0.0];
pmax = [1850 2.4 1.0];
clim1 = pmin(1)+cumsum((pmax(1)-pmin(1))/NC1*ones(1,NC1));
clim2 = pmin(1)+cumsum((pmax(1)-pmin(1))/NC2*ones(1,NC2));
clim3 = pmin(1)+cumsum((pmax(1)-pmin(1))/NC3*ones(1,NC3));
clim4 = pmin(1)+cumsum((pmax(1)-pmin(1))/NC4*ones(1,NC4));

rlim1 = pmin(2)+cumsum((pmax(2)-pmin(2))/NR1*ones(1,NR1));
rlim2 = pmin(2)+cumsum((pmax(2)-pmin(2))/NR2*ones(1,NR2));
rlim3 = pmin(2)+cumsum((pmax(2)-pmin(2))/NR3*ones(1,NR3));
rlim4 = pmin(2)+cumsum((pmax(2)-pmin(2))/NR4*ones(1,NR4));

%%
%% KS TEST
%%

%for i=1:NZ;
%   [hc(i),pc(i),kc(i)] = kstest2(c1(:,i),c2(:,i));
%   [hr(i),pr(i),kr(i)] = kstest2(r1(:,i),r2(:,i));
%   [ha(i),pz(i),ka(i)] = kstest2(a1(:,i),a2(:,i));
%end;

%%
%% Bhattacharyya distance
%%
for iz=1:NZ
   %%
   %% Site 19 & 4
   %%
   Nc19(iz,:) = histc(c19(:,iz),clim1);
   Nc19(iz,:) = Nc19(iz,:)/(sum(Nc19(iz,:)));
   Nc4_a(iz,:) = histc(c4(:,iz),clim1);
   Nc4_a(iz,:) = Nc4_a(iz,:)/(sum(Nc4_a(iz,:)));
   BHc(iz,1) = sum(sqrt(Nc19(iz,:).*Nc4_a(iz,:)));
   %%
   %%
   %% Site 4 & 5
   %%
   Nc4_b(iz,:) = histc(c4(:,iz),clim2);
   Nc4_b(iz,:) = Nc4_b(iz,:)/(sum(Nc4_b(iz,:)));
   Nc5_a(iz,:) = histc(c5(:,iz),clim2);
   Nc5_a(iz,:) = Nc5_a(iz,:)/(sum(Nc5_a(iz,:)));
   BHc(iz,2) = sum(sqrt(Nc4_b(iz,:).*Nc5_a(iz,:)));
   %%
   %% Site 5 & 1
   %%
   Nc5_b(iz,:) = histc(c5(:,iz),clim3);
   Nc5_b(iz,:) = Nc5_b(iz,:)/(sum(Nc5_b(iz,:)));
   Nc1_a(iz,:) = histc(c1(:,iz),clim3);
   Nc1_a(iz,:) = Nc1_a(iz,:)/(sum(Nc1_a(iz,:)));
   BHc(iz,3) = sum(sqrt(Nc5_b(iz,:).*Nc1_a(iz,:)));
   %%
   %% Site 1 & 21
   %%
   Nc1_b(iz,:) = histc(c1(:,iz),clim4);
   Nc1_b(iz,:) = Nc1_b(iz,:)/(sum(Nc1_b(iz,:)));
   Nc21(iz,:) = histc(c21(:,iz),clim4);
   Nc21(iz,:) = Nc21(iz,:)/(sum(Nc21(iz,:)));
   BHc(iz,4) = sum(sqrt(Nc1_b(iz,:).*Nc21(iz,:)));

   %%
   %% Site 19 & 4
   %%
   Nr19(iz,:) = histc(r19(:,iz),rlim1);
   Nr19(iz,:) = Nr19(iz,:)/(sum(Nr19(iz,:)));
   Nr4_a(iz,:) = histc(r4(:,iz),rlim1);
   Nr4_a(iz,:) = Nr4_a(iz,:)/(sum(Nr4_a(iz,:)));
   BHr(iz,1) = sum(sqrt(Nr19(iz,:).*Nr4_a(iz,:)));
   %%
   %%
   %% Site 4 & 5
   %%
   Nr4_b(iz,:) = histc(r4(:,iz),rlim2);
   Nr4_b(iz,:) = Nr4_b(iz,:)/(sum(Nr4_b(iz,:)));
   Nr5_a(iz,:) = histc(r5(:,iz),rlim2);
   Nr5_a(iz,:) = Nr5_a(iz,:)/(sum(Nr5_a(iz,:)));
   BHr(iz,2) = sum(sqrt(Nr4_b(iz,:).*Nr5_a(iz,:)));
   %%
   %% Site 5 & 1
   %%
   Nr5_b(iz,:) = histc(r5(:,iz),rlim3);
   Nr5_b(iz,:) = Nr5_b(iz,:)/(sum(Nr5_b(iz,:)));
   Nr1_a(iz,:) = histc(r1(:,iz),rlim3);
   Nr1_a(iz,:) = Nr1_a(iz,:)/(sum(Nr1_a(iz,:)));
   BHr(iz,3) = sum(sqrt(Nr5_b(iz,:).*Nr1_a(iz,:)));
   %%
   %% Site 1 & 21
   %%
   Nr1_b(iz,:) = histc(r1(:,iz),rlim4);
   Nr1_b(iz,:) = Nr1_b(iz,:)/(sum(Nr1_b(iz,:)));
   Nr21(iz,:) = histc(r21(:,iz),rlim4);
   Nr21(iz,:) = Nr21(iz,:)/(sum(Nr21(iz,:)));
   BHr(iz,4) = sum(sqrt(Nr1_b(iz,:).*Nr21(iz,:)));


end;

ii(1,:) = [11, 52, 112, 223];
ii(2,:) = [4,  27, 212, 335];
ii(3,:) = [9,  56, 181, 349];
ii(4,:) = [7,  38, 139, 304];

nx = 4;
ny = 1;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

figure(1);hold on;box on;
for i=1:4;
   subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   hold on; box on;
   plot(BHc(:,i),z,'k');
   for j=1:4;
      plot([0,1],[z(ii(i,j)),z(ii(i,j))],'--k');
   end;
   plot([0.3,0.3],[0,4],':k');
   set(gca,'YDir','reverse','XLim',[0,1]);
   set(gca,'XTick',[0.2 0.4 0.6 0.8]);
   if(i>1);
      set(gca,'YTickLabel',[]);
   end;
   if(i==1);
      ylabel('Depth (m)');
   end;
   xlabel('BC');
   if(i==1);title('Sites 19 and 4');end;
   if(i==2);title('Sites 4 and 5');end;
   if(i==3);title('Sites 5 and 1');end;
   if(i==4);title('Sites 1 and 21');end;
end

figure(2);hold on;box on;
for i=1:4;
   subplot('Position',[loc(1,i) loc(2,i) spw sph]);
   hold on; box on;
   plot(BHr(:,i),z,'k');
   for j=1:4;
      plot([0,1],[z(ii(i,j)),z(ii(i,j))],'--k');
   end;
   plot([0.3,0.3],[0,4],':k');
   set(gca,'YDir','reverse','XLim',[0,1]);
   set(gca,'XTick',[0.2 0.4 0.6 0.8]);
   if(i>1);
      set(gca,'YTickLabel',[]);
   end;
   if(i==1);
      ylabel('Depth (m)');
   end;
   xlabel('BC');
   if(i==1);title('Sites 19 and 4');end;
   if(i==2);title('Sites 4 and 5');end;
   if(i==3);title('Sites 5 and 1');end;
   if(i==4);title('Sites 1 and 21');end;
end


nx = 4;
ny = 4;
xim = 0.04;
yim = 0.04;
xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

figure(3);hold on;box on;
format short;
for i=1:4;
   subplot('Position',[loc(1,((i-1)*4)+1) loc(2,((i-1)*4)+1) spw sph]);hold on; box on;
   F1 = stairs(clim1,Nc19(ii(1,i),:),'k');
   hold on;
   F2 = stairs(clim1,Nc4_a(ii(1,i),:),':k');
   legend(['BC = ',num2str(BHc(ii(1,i),1),'%4.2f')]);
   legend('boxoff');
   lolim  = min([min(c19(:,ii(1,i))),min( c4(:,ii(1,i)))]);
   uplim  = max([max(c19(:,ii(1,i))),max( c4(:,ii(1,i)))]);
   set(gca,'XLim',[1450,1500]);
   set(gca,'XTick',[1460 1480 1500]);
   if(i==1);title('Sites 19 and 4');end;
   if(i==4);xlabel('Sound speed (m/s)');end;
   ylabel('Probability');
end;
for i=1:4;
   subplot('Position',[loc(1,((i-1)*4)+2) loc(2,((i-1)*4)+2) spw sph]);hold on; box on;
   F1 = stairs(clim2,Nc4_b(ii(2,i),:),'k');
   hold on;
   F2 = stairs(clim2,Nc5_a(ii(2,i),:),':k');
   legend(['BC = ',num2str(BHc(ii(2,i),2),'%4.2f')]);
   legend('boxoff');
   lolim  = min([min(c4(:,ii(2,i))),min( c5(:,ii(2,i)))]);
   uplim  = max([max(c4(:,ii(2,i))),max( c5(:,ii(2,i)))]);
   set(gca,'XLim',[lolim,uplim]);
   if(i==1);title('Sites 4 and 5');end;
   if(i==4);xlabel('Sound speed (m/s)');end;
end;
for i=1:4;
   subplot('Position',[loc(1,((i-1)*4)+3) loc(2,((i-1)*4)+3) spw sph]);hold on; box on;
   F1 = stairs(clim3,Nc5_b(ii(3,i),:),'k');
   hold on;
   F2 = stairs(clim3,Nc1_a(ii(3,i),:),':k');
   legend(['BC = ',num2str(BHc(ii(3,i),3),'%4.2f')]);
   legend('boxoff');
   lolim  = min([min(c5(:,ii(3,i))),min( c1(:,ii(3,i)))]);
   uplim  = max([max(c5(:,ii(3,i))),max( c1(:,ii(3,i)))]);
   set(gca,'XLim',[lolim,uplim]);
   if(i==1);title('Sites 5 and 1');end;
   if(i==4);xlabel('Sound speed (m/s)');end;
end;
for i=1:4;
   subplot('Position',[loc(1,((i-1)*4)+4) loc(2,((i-1)*4)+4) spw sph]);hold on; box on;
   F1 = stairs(clim4,Nc1_b(ii(4,i),:),'k');
   hold on;
   F2 = stairs(clim4,Nc21(ii(4,i),:),':k');
   legend(['BC = ',num2str(BHc(ii(4,i),4),'%4.2f')]);
   legend('boxoff');
   lolim  = min([min(c1(:,ii(4,i))),min( c21(:,ii(4,i)))]);
   uplim  = max([max(c1(:,ii(4,i))),max( c21(:,ii(4,i)))]);
   set(gca,'XLim',[lolim,uplim]);
   if(i==1);title('Sites 1 and 21');end;
   if(i==4);xlabel('Sound speed (m/s)');end;
end;

%%
%% Bhattacharyya distance (special case; track high speed layer)
%%
Nc5 = histc(c5(:,276),clim3);
Nc5 = Nc5/(sum(Nc5));

Nc1 = histc(c1(:,124),clim3);
Nc1 = Nc1/(sum(Nc1));

clear Nc21;
Nc21 = histc(c21(:,78),clim3);
Nc21 = Nc21/(sum(Nc21));

BHc2(1) = sum(sqrt(Nc5.*Nc1));
BHc2(2) = sum(sqrt(Nc1.*Nc21));
BHc2(3) = sum(sqrt(Nc5.*Nc21));

figure(4);hold on;box on;
F1 = stairs(clim3,Nc5,'k');
F2 = stairs(clim3,Nc1,'--k');
F3 = stairs(clim3,Nc21,':k');
legend(['BC = ',num2str(BHc2(1),'%4.2f')],...
       ['BC = ',num2str(BHc2(2),'%4.2f')],...
       ['BC = ',num2str(BHc2(3),'%4.2f')]);
legend('boxoff');
set(gca,'XLim',[1580,1850]);
xlabel('Sound speed (m/s)');
ylabel('Probability');

%%
%% Bhattacharyya distance (special case; track high speed layer)
%%
%Nr5 = histc(r5(:,276),rlim3);
%Nr5 = Nr5/(sum(Nr5));
%
%Nr1 = histc(r1(:,124),rlim3);
%Nr1 = Nr1/(sum(Nr1));
%
%clear Nr21;
%Nr21 = histc(r21(:,78),rlim3);
%Nr21 = Nr21/(sum(Nr21));
%
%BHr2(1) = sum(sqrt(Nr5.*Nr1));
%BHr2(2) = sum(sqrt(Nr1.*Nr21));
%BHr2(3) = sum(sqrt(Nr5.*Nr21));
%
%figure(5);hold on;box on;
%F1 = stairs(rlim3,Nr5,'k');
%F2 = stairs(rlim3,Nr1,'--k');
%F3 = stairs(rlim3,Nr21,':k');
%legend(['BC = ',num2str(BHr2(1),'%4.2f')],...
%       ['BC = ',num2str(BHr2(2),'%4.2f')],...
%       ['BC = ',num2str(BHr2(3),'%4.2f')]);
%legend('boxoff');
%set(gca,'XLim',[1.2,2.4]);
%xlabel('Sound speed (m/s)');
%ylabel('Probability');


return;
