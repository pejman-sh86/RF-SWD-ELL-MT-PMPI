
filename = 'sim_illapel_sample.mat';
load(filename);
filebase = strrep(filename,'_sample.mat','');
imagefile = strcat(filebase,'_result.mat'); 

isyn = 1;
ilin = 1;
islip = 1;   % Invert directly for slip (surface wave inversion)

plot_tw1 = strcat(filebase,'_sl1_marg.');
plot_tw2 = strcat(filebase,'_sl2_marg.');
plot_med = strcat(filebase,'_sl_med.');
plot_mean = strcat(filebase,'_sl_mean.');
plot_LP = strcat(filebase,'_logL_logP.');
plotext1    = 'fig';
plotext2    = 'png';

parfile  = strcat(filebase,'_parameter.dat');
datfile  = strcat(filebase,'.hdf5');

NSTN    = h5readatt(datfile,'/Observed_data','N_sta');
deltt   = h5readatt(datfile,'/Observed_data','Sample_rate');
Ntraces = h5read(datfile,'/Observed_data/Ntraces');
dobs = h5read(datfile,'/Observed_data/displacements');
dtru = h5read(datfile,'/Observed_data/displacements_noise_free');
sdpar = h5read(datfile,'/LinearSolution/std_dev_stn');
model = h5read(datfile,'/LinearSolution/model');
replin = h5read(datfile,'/LinearSolution/pred_linear');

NRAN    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx     = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx     = h5readatt(datfile,'/Sensitivity_kernel','max_y');
mu = h5read(datfile,'/Rigidity/mu');
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
dr = Rmx/NRAN;
dz = Zmx/NDEP;
noise = dobs-dtru;
Ntraces = cast(Ntraces,'like',1);


sfarea = dr*dz;
if(islip == 0);
  scale = mu.*1.e6*sfarea; 
  scale = repmat(scale,1,NRAN)./1.e21;
else;
  scale = ones(NDEP,NRAN);
end;
logL = A(:,1);
logP = A(:,1);
lambda = A(:,428);

%%
k=A(1,4);
NSMP = size(A,1);
ist =7;
sl1tmp = A(:,ist:4:ist+(k-1)*4);
sl2tmp = A(:,ist+1:4:ist+1+(k-1)*4);

sl1 = zeros(NDEP,NRAN,NSMP);
sl2 = zeros(NDEP,NRAN,NSMP);

for ismp=1:NSMP;
  ik = 1;
  for ir=1:NRAN;
  for iz=1:NDEP;
    sl1(iz,ir,ismp) = exp(sl1tmp(ismp,ik));
    sl2(iz,ir,ismp) = exp(sl2tmp(ismp,ik));
    sl1_log(iz,ir,ismp) = sl1tmp(ismp,ik);
    sl2_log(iz,ir,ismp) = sl2tmp(ismp,ik);
    ik = ik + 1;
  end;
  end;
  sl1(:,:,ismp) = sl1(:,:,ismp)./scale;
  sl2(:,:,ismp) = sl2(:,:,ismp)./scale;
end;
ik = 1;
for ir=1:NRAN;
  for iz=1:NDEP;
    sltru1(iz,ir) = model(ik);
    sltru2(iz,ir) = model(ik+(NRAN*NDEP));
    ik = ik + 1;
  end;
end;
sltru1 = sltru1./scale;
sltru2 = sltru2./scale;
save(imagefile,'sl1','sl2','logL','logP','lambda');
%%
%% From MCMC inversion
%%
rep=dlmread('sim_illapel_replica.dat');
rep(1:2,:)=[];
res = dobs-rep';

reslin = dobs-replin;

%d(1).dobs = dobs(1:Ntraces(1));
%d(1).noise = noise(1:Ntraces(1));
for ist=1:NSTN;
  d(ist).dobs = dobs(sum(Ntraces(1:ist-1))+1:sum(Ntraces(1:ist)));
  d(ist).noise = noise(sum(Ntraces(1:ist-1))+1:sum(Ntraces(1:ist)));
  d(ist).noise2 = d(ist).noise/(sdpar(ist));
  d(ist).res = res(sum(Ntraces(1:ist-1))+1:sum(Ntraces(1:ist)));
  d(ist).res2 = d(ist).res/(sdpar(ist));
  d(ist).reslin = reslin(sum(Ntraces(1:ist-1))+1:sum(Ntraces(1:ist)));
  d(ist).reslin2 = d(ist).reslin/(sdpar(ist));
  
  norm(ist) = sum(d(ist).noise2.^2);
  norm2(ist) = sum(d(ist).res2.^2);
  normlin(ist) = sum(d(ist).reslin2.^2);
end;

%% Median and mean plots
sltot = sqrt(sl1.^2+sl2.^2);

figmed = figure();
imagesc(median(sltot,3));colorbar;
print(figmed,'-painters','-r250',strcat(plot_med,plotext2),'-dpng');

figmean = figure();
imagesc(mean(sltot,3));colorbar;
print(figmean,'-painters','-r250',strcat(plot_mean,plotext2),'-dpng');

figLP = figure();
h(1)=subaxis(1,3,1,'Spacing',0.1,'Padding',0,...
            'ML', 0.1,'MR',0.02,'MB',.1,'MT',.02);
hold on; box on;
[n1,out]=hist(squeeze(A(:,1)),30);
n1=n1/trapz(out,n1);
n1 = [0, n1, 0];out = [out(1) out out(end)];
[xx,yy]=stairs(out,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(out,n1,'k');
axis tight
xlabel('log(Likelihood)=-1/2\chi^2')

h(2)=subaxis(1,3,2,'Spacing',0.1,'Padding',0,...
            'ML', 0.1,'MR',0.02,'MB',.1,'MT',.02);
hold on; box on;
[n1,out]=hist(squeeze(A(:,2)),30);
n1=n1/trapz(out,n1);
n1 = [0, n1, 0];out = [out(1) out out(end)];
[xx,yy]=stairs(out,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(out,n1,'k');
axis tight
xlabel('log(Prior)=1/\lambda^2 |LM|^2')

h(3)=subaxis(1,3,3,'Spacing',0.1,'Padding',0,...
            'ML', 0.1,'MR',0.02,'MB',.1,'MT',.02);
hold on; box on;
[n1,out]=hist(A(100:end,428),30);
n1=n1/trapz(out,n1);
n1 = [0, n1, 0];out = [out(1) out out(end)];
[xx,yy]=stairs(out,n1,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(out,n1,'k');
axis tight
xlabel('lambda')

print(figLP,'-painters','-r250',strcat(plot_LP,plotext2),'-dpng');

%% Closed-form Uncertainties

if(ilin == 1);
  mu=load('lin_post_mean.txt');
  sigma=load('lin_post_cov.txt');
  xx_ln = [0.:.001:10.];
  ipar = 1;
  for ir = 1:NRAN;
    for iz = 1:NDEP;
      ipar2 = ipar+NRAN*NDEP;
      mu1 = mu(ipar);
      mu2 = mu(ipar2);
      var1 = sigma(ipar)^2;
      var2 = sigma(ipar2)^2;
      lng1(:,iz,ir) = 1./(xx_ln.*sqrt(2.*pi).*sqrt(var1)).*...
                      exp(-1./(2.*var1).*(log(xx_ln)-mu1).^2);
      lng2(:,iz,ir) = 1./(xx_ln.*sqrt(2.*pi).*sqrt(var2)).*...
                      exp(-1./(2.*var2).*(log(xx_ln)-mu2).^2);
      ipar = ipar + 1;
    end;
  end;
else
  xx_ln = [0.:.001:10.];
  lng1 = zeros(length(xx_ln),NDEP,NRAN);
  lng2 = zeros(length(xx_ln),NDEP,NRAN);
end;

xx = [-20.:.01:20.];
[fsl1]=ffi_fix_slip_uncertainty('sim_illapel_sample.mat',isyn,ilin,sl1,xx_ln,scale,lng1,sltru1);
[fsl2]=ffi_fix_slip_uncertainty('sim_illapel_sample.mat',isyn,ilin,sl2,xx_ln,scale,lng2,sltru2);
[fsl1_log]=ffi_fix_slip_uncertainty_log('sim_illapel_sample.mat',isyn,ilin,sl1_log,xx);
[fsl2_log]=ffi_fix_slip_uncertainty_log('sim_illapel_sample.mat',isyn,ilin,sl2_log,xx);

print(fsl1,'-painters','-r250',strcat(plot_tw1,plotext2),'-dpng');
%saveas(fsl2,strcat(plot_tw1,plotext1),'fig');
print(fsl2,'-painters','-r250',strcat(plot_tw2,plotext2),'-dpng');
%saveas(fsl2,strcat(plot_tw1,plotext1),'fig');
