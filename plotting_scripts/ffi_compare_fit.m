function []=ffi_compare_fit(filebase);

datfile1 = strcat(filebase,'.hdf5');
datfile2 = strcat('../maule_sub4_50x25/',filebase,'.hdf5');
repfile1 = strcat(filebase,'_replica.dat');
repfile2 = strcat('../maule_sub4_50x25/',filebase,'_replica.dat');

NSTN1  = h5readatt(datfile1,'/Observed_data','N_sta');
deltt1 = h5readatt(datfile1,'/Observed_data','Sample_rate');
dat1 = h5read(datfile1,'/Observed_data/displacements');
dat1 = dat1';
NTSMP1 = h5read(datfile1,'/Observed_data/Ntraces');

NSTN2  = h5readatt(datfile2,'/Observed_data','N_sta');
deltt2 = h5readatt(datfile2,'/Observed_data','Sample_rate');
dat2 = h5read(datfile2,'/Observed_data/displacements');
dat2 = dat2';
NTSMP2 = h5read(datfile2,'/Observed_data/Ntraces');

rep1=dlmread(repfile1);
rep2=dlmread(repfile2);

res1 = dat1(1,:)-rep1(3,:);
res2 = dat2(1,:)-rep2(3,:);
sum(abs(res1))/sum(NTSMP1)
sum(abs(res2))/sum(NTSMP2)

  %%
  %% Full Data plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  for i=1:NSTN1;
    subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    iend1 = sum(NTSMP1(1:i));
    istart1 = iend1-NTSMP1(i)+1;
    t1(i).t = deltt1*[0:NTSMP1(i)-1];
    plot(t1(i).t,dat1(1,istart1:iend1),'b','LineWidth',1.);
    plot(t1(i).t,rep1(3,istart1:iend1),'--r','LineWidth',1.);
    plot([t1(i).t(1) t1(i).t(end)],[0 0],'--k','LineWidth',1);

    iend2 = sum(NTSMP2(1:i));
    istart2 = iend2-NTSMP2(i)+1;
    t2 = deltt2*[0:NTSMP2(i)-1];
    %plot(t2,dat2(1,istart2:iend2),'.k','LineWidth',1.);
    plot(t2,rep2(3,istart2:iend2),'--k','LineWidth',1.);

    miny = min(dat1(1,istart1:iend1))-(-min(dat1(1,istart1:iend1))+max(dat1(1,istart1:iend1)))/10;
    maxy = max(dat1(1,istart1:iend1))+(-min(dat1(1,istart1:iend1))+max(dat1(1,istart1:iend1)))/10;
    set(gca,'YLim',[miny maxy],'XLim',[t1(i).t(1) t1(i).t(end)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  %%
  %% Full res plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  for i=1:NSTN1;
    subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    iend1 = sum(NTSMP1(1:i));
    istart1 = iend1-NTSMP1(i)+1;
    plot(t1(i).t,res1(1,istart1:iend1),'-r','LineWidth',1.);
    plot([t1(i).t(1) t1(i).t(end)],[0 0],'--k','LineWidth',1);

    miny(i) = min(res1(1,istart1:iend1))-(-min(res1(1,istart1:iend1))+max(res1(1,istart1:iend1)))/10;
    maxy(i) = max(res1(1,istart1:iend1))+(-min(res1(1,istart1:iend1))+max(res1(1,istart1:iend1)))/10;
    set(gca,'YLim',[miny(i) maxy(i)],'XLim',[t1(i).t(1) t1(i).t(end)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;
  %%
  %% Full res plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 6;
  ny = 5;
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  for i=1:NSTN1;
    subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
    set(gca,'FontSize',14,'layer','top','LineWidth',1)

    iend2 = sum(NTSMP2(1:i));
    istart2 = iend2-NTSMP2(i)+1;
    t2 = deltt2*[0:NTSMP2(i)-1];
    plot(t2,res2(1,istart2:iend2),'b','LineWidth',1.);
    plot([t2(1) t2(end)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[miny(i) maxy(i)],'XLim',[t1(i).t(1) t1(i).t(end)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:200:1500]);
    if(i == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
  end;





return;
