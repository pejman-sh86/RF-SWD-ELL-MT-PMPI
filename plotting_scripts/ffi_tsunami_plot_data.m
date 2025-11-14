function []=ffi_tsunami_plot_data(filebase);

isyn = 1;
datfile  = strcat(filebase,'.hdf5');
parfile  = strcat(filebase,'_parameter.dat');
repfile  = strcat(filebase,'_replica.dat');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = cast(h5readatt(datfile,'/Observed_data','Sample_rate'),'like',1.);
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc_Cart');
dhyp   = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin  = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax  = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
NTW    = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
NTSMP  = cast(h5read(datfile,'/Observed_data/Ntraces'),'like',1);

dat    = h5read(datfile,'/Observed_data/displacements');
NDAT   = length(dat);
rep    = dlmread(repfile);
replin = h5read(datfile,'/Sensitivity_kernel/synthetic_displacements');

if(isyn == 1);
  dat = dat';
  replin = replin';
end;
%%
%% Full Data plots:
%%
figdat = figure();hold on;box on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 12])
nx = 3;
ny = 5;
xim = 0.01/nx;
yim = 0.01/ny;
xymarg = [0.07 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

jj = 1;
for i = 1:NSTN;
  if(i == (nx*ny)+1);
    figdat2 = figure();hold on;box on;
    jj = 1;
  end;
  subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
  set(gca,'FontSize',14,'layer','top','LineWidth',1)

  iend = sum(NTSMP(1:i));
  istart = iend-NTSMP(i)+1;
  plot(deltt*[0:NTSMP(i)-1],dat(1,istart:iend),'LineWidth',1.5);
  plot(deltt*[0:NTSMP(i)-1],rep(3,istart:iend),'-r','LineWidth',1.5);
  %plot(deltt*[0:NTSMP(i)-1],replin(1,istart:iend),'--k','LineWidth',1.5);
  plot([0 deltt*NTSMP(i)],[0 0],'--k','LineWidth',1);

  set(gca,'YLim',[min([min(dat(1,istart:iend)),min(rep(3,istart:iend)),...
      min(replin(1,istart:iend))]) max([max(dat(1,istart:iend)),...
      max(rep(3,istart:iend)),max(replin(1,istart:iend))])],'XLim',[0 deltt*NTSMP(i)]);
  xlabel('time (s)');
  set(gca,'XTick',[0:300:1600]);
  if(i == 1);
    ylabel('Amplitude');
  else;
    set(gca,'YTickLabel',[]);
  end;
  box on;
  jj = jj + 1;
end;
legend('obs','pred','lin');
return;
