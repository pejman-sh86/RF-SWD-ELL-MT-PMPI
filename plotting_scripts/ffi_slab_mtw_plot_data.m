%function []=ffi_tsunami_plot_data(filebase);
filebase = 'tohoku_data';

isave = 1;
ieps  = 1;
ipng  = 0;
ilin  = 0;

parfile     = strcat(filebase,'_parameter.dat');
edgefile    = strcat(filebase,'_fault_edge.txt');
M0file      = strcat(filebase,'_M0.dat');
datfile     = strcat(filebase,'.hdf5');
datfilecGPS = strcat(filebase,'_cGPS.hdf5');
covfile     = strcat(filebase,'_covmat.mat');

[IMAP,ICOV,I_WP,I_cGPS,I_GPS,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,NGPS]=ffi_read_parfile(parfile);

datfilecGPS  = strcat(filebase,'_cGPS.hdf5');
datfilecGPSs  = strcat(filebase,'_cGPSs.hdf5');

%%
%% Output plot files
%% 
plotfiledatW= strcat(filebase,'_data_W.');
plotfiledatcGPSN = strcat(filebase,'_data_cGPSN.');
plotfiledatcGPSE = strcat(filebase,'_data_cGPSE.');
plotfiledatcGPSZ = strcat(filebase,'_data_cGPSZ.');
plotfiledatcGPSs  = strcat(filebase,'_data_cGPSs.');
plotfiledatGPS  = strcat(filebase,'_data_GPS.');

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc');
rake   = h5readatt(datfile,'/Sensitivity_kernel','rake_comp');
dhyp   = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NTW    = h5readatt(datfile,'/Sensitivity_kernel','N_time_win');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin  = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax  = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
mu     = h5read(datfile,'/Rigidity/mu');
sfgrid = h5read(datfile,'/Sensitivity_kernel/subfaultgrid');
wp_list = h5read(datfile,'/Observed_data/traces_list');

NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
NTW  = cast(NTW,'like',1);

NSF = NRAN*NDEP;
NFPMX = NVMX*NPV*NTW;

if(I_WP == 1);
  rep=dlmread('tohoku_data_replica.dat');
  obs=dlmread('tohoku_data_observed.dat');
  NSTN = rep(1,1);
  NTSMP = cast(rep(2,1:NSTN),'like',1);
  cNTSMP=[0,cumsum(NTSMP)];
  obs(1:2,:)=[];
  rep(1:2,:)=[];
  rep_lin= h5read(datfile,'/Sensitivity_kernel/synthetic_displacements');
  for istn=1:NSTN;
      wp(istn).obs     = obs(cNTSMP(istn)+1:cNTSMP(istn+1));
      wp(istn).rep     = rep(cNTSMP(istn)+1:cNTSMP(istn+1));
      wp(istn).rep_lin = rep_lin(cNTSMP(istn)+1:cNTSMP(istn+1));
      wp(istn).name    = wp_list(istn);
  end;
end;
if(I_cGPS == 1);
  repc=dlmread('tohoku_data_replicacGPS.dat');
  obsc=dlmread('tohoku_data_observedcGPS.dat');
  NcGPS = repc(1,1);
  NTSMPcGPS = repc(2,1:NcGPS);
  cNTSMPcGPS=[0,cumsum(NTSMPcGPS)];
  obsc(1:2,:)=[];
  repc(1:2,:)=[];
  repc_lin= h5read(datfilecGPS,'/Sensitivity_kernel/synthetic_displacements');
  cGPS_list = h5read(datfilecGPS,'/Observed_data/traces_list');

  istn2 = 1;
  for istn=1:NcGPS/3;
    istn2
    cG(istn).obse     = obsc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repe     = repc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repe_lin = repc_lin(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).namee    = cGPS_list(istn2);
    istn2 = istn2 + 1;
    
    istn2
    cG(istn).obsn     = obsc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repn     = repc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repn_lin = repc_lin(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).namen    = cGPS_list(istn2);
    istn2 = istn2 + 1;
    
    istn2
    cG(istn).obsz     = obsc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repz     = repc(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).repz_lin = repc_lin(cNTSMPcGPS(istn2)+1:cNTSMPcGPS(istn2+1));
    cG(istn).namez    = cGPS_list(istn2);
    istn2 = istn2 + 1;
  end;
  repcs=dlmread('tohoku_data_replicacGPSs.dat');
  obscs=dlmread('tohoku_data_observedcGPSs.dat');
  obscs(1,:)=[];
  repcs(1,:)=[];
  repcs_lin= h5read(datfilecGPSs,'/Sensitivity_kernel/synthetic_displacements');
  rescs = obscs-repcs;
  sig_cs(1) = std(rescs(1:3:end));
  sig_cs(2) = std(rescs(2:3:end));
  sig_cs(3) = std(rescs(3:3:end));
end;
if(I_GPS == 1);
  repg=dlmread('tohoku_data_replicaGPS.dat');
  obsg=repg(3,:);
  repg(3,:)=[];
  repg(1,:)=[];
  resg = obsg-repg;
  NGPS = size(repg,2);
  NGPS2 = NGPS/3;
  sig_g(1) = std(resg(1:NGPS2));
  sig_g(2) = std(resg(NGPS2+1:2*NGPS2));
  sig_g(3) = std(resg(2*NGPS2+1:3*NGPS2));
end;

wid = 14;
hei = 10;

if(I_WP == 1);
  nx = 3;
  ny = 3;
  ML = .03; MR = .03; MB = .02; MT = .02;
  SP = .02; PB = .04; PAD = 0; FNT = 14; LNWD = 3;
  figw=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wid hei]);
  idx = [1:ceil(NSTN/(nx*ny)):NSTN];
  for iplt=1:nx*ny;
    subaxis(ny,nx,iplt,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
    clear time;
    time = [0:length(wp(idx(iplt)).obs)-1]*deltt;
    plot(time,wp(idx(iplt)).obs,'-b','LineWidth',LNWD);
    plot(time,wp(idx(iplt)).rep,'--r','LineWidth',LNWD);
    if(ilin == 1);
      plot(time,wp(idx(iplt)).rep_lin,':k','LineWidth',LNWD);
    end;
    set(gca,'XLim',[time(1),time(end)]);
    set(gca,'YTickLabel',[]);
    text(0.02,.9,wp(idx(iplt)).name,'Units','normalized','FontSize',FNT)
    if(iplt > nx*(ny-1));
        xlabel('time (s)');
    end;
  end;
  if(isave == 1)
    if(ipng == 1);
      print(figw,'-painters','-r250',strcat(plotfiledatW,plotext2),'-dpng');
      saveas(figw,strcat(plotfiledatW,plotext1),'fig');
    end
    if(ieps == 1);
      print(figw,'-painters','-r250',strcat(plotfiledatW,plotext3),'-depsc');
    end;
  end;
end;


if(I_cGPS == 1);
  %%
  %% cGPS waveform fits
  %%
  nx = 3;
  ny = 3;
  ML = .03; MR = .03; MB = .02; MT = .02;
  SP = .02; PB = .04; PAD = 0; FNT = 14; LNWD = 3;

  figcGPSe=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wid hei]);
  idx = [1:ceil(NcGPS/3/(nx*ny)):NcGPS/3];
  for iplt=1:nx*ny;
    subaxis(ny,nx,iplt,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
    clear time;
    time = [0:length(cG(idx(iplt)).obse)-1]*deltt;
    plot(time,cG(idx(iplt)).obse,'-b','LineWidth',LNWD);
    plot(time,cG(idx(iplt)).repe,'--r','LineWidth',LNWD);
    if(ilin == 1);
      plot(time,cG(idx(iplt)).repe_lin,':k','LineWidth',LNWD);
    end;
    set(gca,'XLim',[time(1),time(end)]);
    set(gca,'YTickLabel',[]);
    text(0.02,.9,cG(idx(iplt)).namee,'Units','normalized','FontSize',FNT);
    if(iplt > nx*(ny-1));
        xlabel('time (s)');
    end;
  end;
  if(isave == 1)
    if(ipng == 1);
      print(figcGPSe,'-painters','-r250',strcat(plotfiledatcGPSE,plotext2),'-dpng');
      saveas(figcGPSe,strcat(plotfiledatcGPSE,plotext1),'fig');
    end;
    if(ieps == 1);
      print(figcGPSe,'-painters','-r250',strcat(plotfiledatcGPSE,plotext3),'-depsc');
    end;
  end;

  figcGPSn=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wid hei]);
  for iplt=1:nx*ny;
    subaxis(ny,nx,iplt,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
    clear time;
    time = [0:length(cG(idx(iplt)).obsn)-1]*deltt;
    plot(time,cG(idx(iplt)).obsn,'-b','LineWidth',LNWD);
    plot(time,cG(idx(iplt)).repn,'--r','LineWidth',LNWD);
    if(ilin == 1);
      plot(time,cG(idx(iplt)).repn_lin,':k','LineWidth',LNWD);
    end;
    set(gca,'XLim',[time(1),time(end)]);
    set(gca,'YTickLabel',[]);
    text(0.02,.9,cG(idx(iplt)).namen,'Units','normalized','FontSize',FNT);
    if(iplt > nx*(ny-1));
        xlabel('time (s)');
    end;
  end;
  if(isave == 1)
    if(ipng == 1);
      print(figcGPSn,'-painters','-r250',strcat(plotfiledatcGPSN,plotext2),'-dpng');
      saveas(figcGPSn,strcat(plotfiledatcGPSN,plotext1),'fig');
    end;
    if(ieps == 1);
      print(figcGPSn,'-painters','-r250',strcat(plotfiledatcGPSN,plotext3),'-depsc');
    end;
  end;

  figcGPSz=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wid hei]);
  for iplt=1:nx*ny;
    subaxis(ny,nx,iplt,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
    set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
    clear time;
    time = [0:length(cG(idx(iplt)).obsz)-1]*deltt;
    plot(time,cG(idx(iplt)).obsz,'-b','LineWidth',LNWD);
    plot(time,cG(idx(iplt)).repz,'--r','LineWidth',LNWD);
    if(ilin == 1);
      plot(time,cG(idx(iplt)).repz_lin,':k','LineWidth',LNWD);
    end;
    set(gca,'XLim',[time(1),time(end)]);
    set(gca,'YTickLabel',[]);
    text(0.02,.9,cG(idx(iplt)).namez,'Units','normalized','FontSize',FNT);
    if(iplt > nx*(ny-1));
        xlabel('time (s)');
    end;
  end;
  if(isave == 1)
    if(ipng == 1);
      print(figcGPSz,'-painters','-r250',strcat(plotfiledatcGPSZ,plotext2),'-dpng');
      saveas(figcGPSz,strcat(plotfiledatcGPSZ,plotext1),'fig');
    end;
    if(ieps == 1);
      print(figcGPSz,'-painters','-r250',strcat(plotfiledatcGPSZ,plotext3),'-depsc');
    end;
  end;
  
  %%
  %% cGPS static offsets
  %%
  nx = 1;
  ny = 3;
  ML = .03; MR = .03; MB = .02; MT = .02;
  SP = .02; PB = .04; PAD = 0; FNT = 14; LNWD = 3;
  figcGPSs=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12]);
  subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obscs(1:3:end),'o-b','LineWidth',LNWD);
  plot(repcs(1:3:end),'o--r','LineWidth',LNWD);
  if(ilin == 1);
    plot(repcs_lin(1:3:end),'o:k','LineWidth',LNWD);
  end;

  subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obscs(2:3:end),'o-b','LineWidth',LNWD);
  plot(repcs(2:3:end),'o--r','LineWidth',LNWD);
  if(ilin == 1);
    plot(repcs_lin(2:3:end),'o:k','LineWidth',LNWD);
  end;
  
  subaxis(ny,nx,3,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obscs(3:3:end),'o-b','LineWidth',LNWD);
  plot(repcs(3:3:end),'o--r','LineWidth',LNWD);
  if(ilin == 1);
    plot(repcs_lin(3:3:end),'o:k','LineWidth',LNWD);
  end;
  if(isave == 1)
    if(ipng == 1);
      print(figcGPSs,'-painters','-r250',strcat(plotfiledatcGPSs,plotext2),'-dpng');
      saveas(figcGPSs,strcat(plotfiledatcGPSs,plotext1),'fig');
    end;
    if(ieps == 1);
      print(figcGPSs,'-painters','-r250',strcat(plotfiledatcGPSs,plotext3),'-depsc');
    end;
  end;

end;
if(I_GPS == 1);
  %%
  %% cGPS static offsets
  %%
  nx = 1;
  ny = 3;
  ML = .03; MR = .03; MB = .02; MT = .02;
  SP = .02; PB = .04; PAD = 0; FNT = 12; LNWD = 3;
  figGPS=figure;hold on;box on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12]);
  subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obsg(1:NGPS2),'o-b','LineWidth',LNWD);
  plot(repg(1:NGPS2),'o--r','LineWidth',LNWD);
  set(gca,'YLim',[0,30]);
  
  subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obsg(NGPS2+1:2*NGPS2),'o-b','LineWidth',LNWD);
  plot(repg(NGPS2+1:2*NGPS2),'o--r','LineWidth',LNWD);
  set(gca,'YLim',[-15,15]);
  
  subaxis(ny,nx,3,'Spacing',SP,'Padding',PAD,'PaddingBottom',PB,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
  set(gca,'FontSize',FNT,'layer','top','LineWidth',1);hold on;box on;
  plot(obsg(2*NGPS2+1:3*NGPS2),'o-b','LineWidth',LNWD);
  plot(repg(2*NGPS2+1:3*NGPS2),'o--r','LineWidth',LNWD);
  set(gca,'YLim',[-15,15]);
  if(isave == 1)
    if(ipng == 1);
      print(figGPS,'-painters','-r250',strcat(plotfiledatGPS,plotext2),'-dpng');
      saveas(figGPS,strcat(plotfiledatGPS,plotext1),'fig');
    end;
    if(ieps == 1);
      print(figGPS,'-painters','-r250',strcat(plotfiledatGPS,plotext3),'-depsc');
    end;
  end;
end;




%return;