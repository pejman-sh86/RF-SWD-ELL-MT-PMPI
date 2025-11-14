%function []=cmt_plot_data(filebase);
%filebase = 'am16_data';
%filebase = 'bo10_data';
%filebase = 'ma10_data';
%filebase = 'mi10_data';
%filebase = 'ak17_6.3';
filebase = 'an17_data';

parfile = strcat(filebase,'_parameter.dat');
obsfile = strcat(filebase,'_observed.dat');
prefile = strcat(filebase,'_replica.dat');

FNT = 14;
deltt = 1;

[IMAP,ICOV,NVMX,ILOC,IAR,IEXCHANGE,IDBLCPL,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay,minMw,minstr,mindip,minrak,...
 maxMw,maxstr,maxdip,maxrak]=cmt_read_parfile(parfile);

[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

pre=dlmread(prefile);

  %%
  %% Full Data plots:
  %%
  figdat = figure();hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
  nx = 5;
  ny = ceil(NSTN/nx);
  xim = 0.01/nx;
  yim = 0.01/ny;
  xymarg = [0.07 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  jjstn = 1;
  for istn = 1:NSTN;
    if(istn == (nx*ny)+1);
      figdat2 = figure();hold on;box on;
      jjstn = 1;
    elseif(istn == (nx*ny)*2+1);
      figdat3 = figure();hold on;box on;
      jjstn = 1;
    end;
    subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
  
    iend = sum(NTSMP(1:istn));
    istart = iend-NTSMP(istn)+1;
    plot(deltt*[0:NTSMP(istn)-1],dobs(istart:iend),'b','LineWidth',1.5);
    plot(deltt*[0:NTSMP(istn)-1],pre(3,istart:iend),'--r','LineWidth',1.5);
    plot([0 deltt*NTSMP(istn)],[0 0],'--k','LineWidth',1);

    set(gca,'YLim',[min(dobs) max(dobs)],'XLim',[0 deltt*NTSMP(istn)]);
    xlabel('time (s)');
    set(gca,'XTick',[0:300:1600]);
    if(istn == 1);
      ylabel('Amplitude');
    else;
      set(gca,'YTickLabel',[]);
    end;
    box on;
    jjstn = jjstn + 1;
  end;

%return;