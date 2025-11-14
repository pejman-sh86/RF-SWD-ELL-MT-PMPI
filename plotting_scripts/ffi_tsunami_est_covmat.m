function [] = ffi_tsunami_est_covmat(filebase);
%%
%% Estimate non-stationary residual covariance matrix for 
%% W phase data to be used in FFI.
%%

time1 = [ 100, 300, 300, 300, 300, 300, 300, 1380, 3000, 3900, 600, 900, 300, 300, 250, 250];
time2 = [3800,3000,2400,2400,3600,2800,3200, 3200, 5600, 6400,3600,3600,1200,2000,2000,2000];  
ylimmax = [ 8, 6, 7, 6, 2.5, 6.2, 4.5, 2, 1.0, 1.0, 1.0, 1.0, 5, 5, 5, 5];
ylimmin = [-2,-8,-7,-4,-1.5,-4.0,-2.0,-1,-0.5,-0.5,-0.5,-0.5,-2,-2,-4,-4];
dtime = ceil((time2-time1)/400)*100;
damp = round((ylimmax-ylimmin)/4);

%% Fraction of data set over which to est non stat std dev:
nfrac = 4.;
%% Damping (raises cosine to that power for Axx damping):
damp_power = 1.2;
%% Switch to turn on non-stationary estimates.
inonstat = 1;
%% Switch to turn on unbiased cov est (usually not a good idea)
iunbiased = 0;
%% Switch to save Cd into file
isave = 1;
%% Simulation switch (simulated exp cov mat with sine std dev non stat)
ISIM  = 0;

datfile     = strcat(filebase,'.hdf5');
%datfile2    = strcat(filebase,'_res.hdf5');
repfile     = strcat(filebase,'_replica.dat');
covmatfile  = strcat(filebase,'_icovmat.dat');
sdfile      = strcat(filebase,'_nonstat_sd.dat');
covmatfile2 = strcat(filebase,'_covmat.mat');

dat        = h5read(datfile,'/Observed_data/displacements');
NSTN       = h5readatt(datfile,'/Observed_data','N_sta');
NTSMP      = h5read(datfile,'/Observed_data/Ntraces');
%pred       = h5read(datfile2,'/Observed_data/pred');
pred       = dlmread(repfile);
deltt   = h5readatt(datfile,'/Observed_data','Sample_rate');
NSTN = cast(NSTN,'like',1);

%% Create HDF5 file for tests:
%[a b] = size(dat);
%h5create(datfile2,'/Observed_data/displacements',[a b]);
%h5write(datfile2,'/Observed_data/displacements',dat);
%h5writeatt(datfile2,'/Observed_data','N_sta',NSTN);
%[a b] = size(NTSMP);
%h5create(datfile2,'/Observed_data/Ntraces',[a b]);
%h5write(datfile2,'/Observed_data/Ntraces',NTSMP);
%[a b] = size(pred);
%h5create(datfile2,'/Observed_data/pred',[a b]);
%h5write(datfile2,'/Observed_data/pred',pred);

%dat        = dat';

%fact = [1.9397, 2.0117, 2.3221, 1.8589, 2.0540, 3.0958,...
%        2.3811, 2.6745, 2.9314, 2.3203, 2.5637, 2.4849,...
%        2.4849, 2.5409, 3.1835, 2.9003, 2.7093, 2.4418,...
%        3.3276, 2.1567, 3.6815, 2.7588, 2.6744, 3.3842,...
%        2.4605, 4.0363, 2.5301, 3.3607, 2.9202, 2.9219];
%fact = 1./(2*fact);
%fact = ones(NSTN,1);
%fact = 0.25;

%save tmp.mat

%% Residual errors:
res = dat-pred(3,:);
%std(res)

nx = 4;
ny = 4;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

if(ISIM == 1);
sd_dat = 0.000107944073408472;
ressim = sd_dat*randn(size(dat));
fig_sim1 = figure();hold on;box on;
fig_sim2 = figure();hold on;box on;
for istn = 1:NSTN;

  iend(istn) = sum(NTSMP(1:istn));
  istart(istn) = iend(istn)-NTSMP(istn)+1;
  clear Csim Lsim;
  for idat = 1:NTSMP(istn);
  for jdat = 1:NTSMP(istn);
    Csim(idat,jdat) = exp(-(1./(NTSMP(istn)/10.))*abs(jdat-idat));
  end;end;
  
  %% Make sine shape for std dev non stat behaviour:
  for idat = istart(istn):iend(istn);
    sig(idat) = sd_dat*(sin(pi/NTSMP(istn)*idat));
  end;
  %res(istart(istn):iend(istn)) = sig(istart(istn):iend(istn)).*res(istart(istn):iend(istn));

  for idat = 1:NTSMP(istn)
    idat2 = istart(istn)+idat-1;
    for jdat = 1:NTSMP(istn)
      jdat2 = istart(istn)+jdat-1;
      Csim(idat,jdat) = Csim(idat,jdat);%*sig(idat2)*sig(jdat2);
    end;
  end;  

  Lsim = chol(Csim,'lower');
  res(istart(istn):iend(istn)) = Lsim*ressim(istart(istn):iend(istn)).';

  figure(fig_sim1)
  subplot('Position',[loc(1,istn) loc(2,istn) spw sph]);hold on;box on;
  imagesc(Csim)
  set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);
  set(gca,'YLim',[1 NTSMP(istn)],'YDir','reverse');

  figure(fig_sim2)
  subplot('Position',[loc(1,istn) loc(2,istn) spw sph]);hold on;box on;
  plot(sig(istart(istn):iend(istn)));
  set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);

end;end;

fig_sd = figure();hold on;box on;
jstn = 1;
%%
%% Estimate non-stationary standard deviation:
%%
for istn = 1:NSTN;
  NAVE = round(double(NTSMP(istn))/nfrac);
  if(NAVE > 40);NAVE = 40;end;
  iend(istn) = sum(NTSMP(1:istn));
  istart(istn) = iend(istn)-NTSMP(istn)+1;
  %% Residuals mean removed
  res_mr(istart(istn):iend(istn)) = res(istart(istn):iend(istn));%-mean(res(istart(istn):iend(istn)));
  mat(istn).res_mr = res_mr(istart(istn):iend(istn));
  sd_stat(istn) = std(res_mr(istart(istn):iend(istn)));
  if(inonstat == 1);
    for idat=istart(istn):iend(istn);
      %NAVE1 = idat-NAVE;
      %NAVE2 = idat+NAVE;
      %if(idat <= istart(istn)+NAVE);
      %  NAVE1 = istart(istn);
      %  NAVE2 = istart(istn)+2*NAVE;
      %end;
      %if(idat >= iend(istn)-NAVE-1);
      %  NAVE1 = iend(istn)-2*NAVE;
      %  NAVE2 = iend(istn);
      %end;
      if(idat <= istart(istn)+NAVE);
        NAVE1 = istart(istn);
      else
        NAVE1 = idat-NAVE;
      end;
      if(idat >= iend(istn)-NAVE-1);
        NAVE2 = iend(istn);
      else
        NAVE2 = idat+NAVE;
      end;

      sd_nonstat(idat) = sqrt(mean(res_mr(NAVE1:NAVE2).^2.));
      %disp([idat,NAVE1,NAVE2,sd_nonstat(idat)]);
    end;
  else;
    sd_nonstat(istart(istn):iend(istn)) = sd_stat(istn);
  end;
  mat(istn).sd_nonstat = sd_nonstat(istart(istn):iend(istn));
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(sd_nonstat(istart(istn):iend(istn)),'-k');
  if(ISIM == 1);plot(sig(istart(istn):iend(istn)),'-b');end;
  plot([0 NTSMP(istn)],[sd_stat(istn) sd_stat(istn)],'--k');
  set(gca,'XLim',[1 NTSMP(istn)],'FontSize',14,'LineWidth',1);
  %if(ISIM == 1);set(gca,'YLim',[sigsim-0.5 sigsim+.5]);end
  jstn = jstn + 1;
end;

%%
%% Standardize residuals with non stat std dev:
%%
jstn = 1;
fig_res = figure();
for istn = 1:NSTN;
  %% Residuals standardized with non stat std dev:
  res_ns(istart(istn):iend(istn)) = res(istart(istn):iend(istn))./sd_nonstat(istart(istn):iend(istn));
  sd_resns(istn) = std(res_ns(istart(istn):iend(istn)));
  %% Residuals standardized with stat std dev:
  res_s(istart(istn):iend(istn)) = res(istart(istn):iend(istn))/sd_stat(istn);
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(res_ns(istart(istn):iend(istn)),'-b');
  plot(res_s(istart(istn):iend(istn)),'--r');
  set(gca,'XLim',[1 NTSMP(istn)],'FontSize',14,'LineWidth',1);
  set(gca,'YLim',[-3 3]);
  jstn = jstn + 1;
end;
%disp('Standard deviation:');
%disp(sd_resns);
%disp('Variance:');
%disp(sd_resns.^2.);

%%
%% Compute autocovariance:
%%
fig_Axx = figure();hold on;box on;
jstn = 1;
for istn = 1:NSTN;
  if(iunbiased == 0);
    [mat(istn).Axx,mat(istn).lags] = xcov(res_ns(istart(istn):iend(istn)));
    mat(istn).Axxn= mat(istn).Axx/double(NTSMP(istn));
  else;
    [mat(istn).Axx,mat(istn).lags] = xcov(res_ns(istart(istn):iend(istn)),'unbiased');
    mat(istn).Axxn= mat(istn).Axx;
  end;
  mat(istn).Axxu= xcov(res_ns(istart(istn):iend(istn)),'unbiased');

  %% Damping to help positive definiteness (all eigenvalues positive, chol. decomp. exixst)
  clear damp;
  damp = cos(pi*double([-NTSMP(istn)+1:NTSMP(istn)-1])./(2*double(NTSMP(istn))-1)).^damp_power;
  %gc0 = figure(1);plot(x);
  mat(istn).Axxn = mat(istn).Axxn .* damp(1:end);
  mat(istn).Axxu = mat(istn).Axxu .* damp(1:end);

  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(mat(istn).lags,mat(istn).Axxn);
  plot(mat(istn).lags,mat(istn).Axxu,'--r');
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'FontSize',14,'LineWidth',1);
  set(gca,'YLim',[-.6 1.2]);
  jstn = jstn + 1;
end;
legend('biased','unbiased')

%%
%% Set up covariance matrix in a diagonal form:
%%
nx = 4;
ny = 4;
fig_Cd = figure();hold on;box on;
for istn = 1:NSTN;
  if(rem(istn,5)==0);disp(istn);end;
  mat(istn).Cd =zeros(NTSMP(istn),NTSMP(istn));
  for irow = 1:NTSMP(istn)
    mat(istn).Cd(irow,:) = mat(istn).Axxn(NTSMP(istn)-(irow-1):end-(irow-1));
  end;
  for idat = 1:NTSMP(istn)
    for jdat = 1:NTSMP(istn)
      mat(istn).Cd(idat,jdat) = mat(istn).Cd(idat,jdat)*...
                                mat(istn).sd_nonstat(idat)*mat(istn).sd_nonstat(jdat);
      %if(idat == jdat);
      %  [idat,jdat,mat(istn).sd_nonstat(idat),mat(istn).sd_nonstat(idat)*mat(istn).sd_nonstat(jdat)]
      %end;
    end;
  end;  
  %mat(istn).Cd = mat(istn).Cd*fact(istn);
  mat(istn).Cd = mat(istn).Cd;%*fact;

  F(istn).Cd = mat(istn).Cd;
  %% Choleschy decomposition, L is lower triangular matrix:
  mat(istn).L = chol(mat(istn).Cd,'lower');
  mat(istn).R = chol(mat(istn).Cd);

  subaxis(ny,nx,istn,'Spacing',0.0,'Padding',0,'Paddingbottom',0.03,...
        'Paddingleft',0.03,'ML', 0.04,'MR',0.01,'MB',.04,'MT',.01);
  hold on;box on;axis equal;

  imagesc(deltt*[1 NTSMP(istn)],deltt*[1 NTSMP(istn)],mat(istn).Cd);
  %axis equal;
  set(gca,'XLim',deltt*[1 NTSMP(istn)],'FontSize',14);
  set(gca,'YLim',deltt*[1 NTSMP(istn)],'YDir','reverse','TickDir','out');
  %set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  %set(gca,'XTick',[time1(istn):dtime(istn):time2(istn)]);
  %set(gca,'YTick',[time1(istn):dtime(istn):time2(istn)]);
  if(istn == 1 |istn == 5 |istn == 9 |istn == 13);
      ylabel('Time (s)');
  end;
  if(istn > 12);
  xlabel('Time (s)');    
  end;
  
end;
NTSMP

%%
%% Check residuals by pre-whitening:
%%
nx = 8;
ny = 4;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
jstn = 1;
fig_Axx_check = figure();hold on; box on;
for istn = 1:NSTN;

  %disp('max L*L^-1')
  %max(max(mat(istn).L*inv(mat(istn).L)))
  %max(max(mat(istn).Cd-mat(istn).L*mat(istn).L'))

  mat(istn).res_pw = inv(mat(istn).L)*res_mr(istart(istn):iend(istn))';
  [mat(istn).Axx_pw,mat(istn).lags_pw] = xcov(mat(istn).res_pw);
  mat(istn).Axx_pw = mat(istn).Axx_pw/double(NTSMP(istn));

  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(mat(istn).lags,mat(istn).Axxn/(max(mat(istn).Axxn)),'-k.')
  plot([mat(istn).lags(1) mat(istn).lags(end)],[0 0],'--k');
  set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);

  subplot('Position',[loc(1,jstn+1) loc(2,jstn+1) spw sph]);hold on;box on;
  plot(mat(istn).lags_pw,mat(istn).Axx_pw/(max(mat(istn).Axx_pw)),'-k.')
  plot([mat(istn).lags_pw(1) mat(istn).lags_pw(end)],[0 0],'--k');
  set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);

  jstn = jstn + 2;
end;


jstn = 1;
fig_hist_check = figure();hold on; box on;
x  = [-6.375:.75:6.375]';
xx = [-6.25:.01:6.25]';
nd = 1./sqrt(2*pi)*exp(-(xx.^2)/2.);
for istn = 1:NSTN;

  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  clear n1 xout;
  [n1,xout] = hist(res_ns(istart(istn):iend(istn)),x);
  n1 = n1/trapz(xout,n1);
  [nx,ny]=stairs(xout,n1,'k');
  patch(nx,ny,[0.8,0.8,0.8]);
  stairs(xout,n1,'k');

  plot(xx,nd,'-k','LineWidth',1);
  set(gca,'XTick',[-4 0 4],'XTickLabel',[],'FontSize',14);
  set(gca,'YTick',[0 0.5],'YTickLabel',[],'FontSize',14);
  axis([-6 6 0 0.6]);

  subplot('Position',[loc(1,jstn+1) loc(2,jstn+1) spw sph]);hold on;box on;
  set(gca,'XTick',[-4 0 4],'FontSize',14);
  set(gca,'XTickLabel',[],'FontSize',14);
  set(gca,'YTick',[0 0.5],'FontSize',14);
  set(gca,'YTickLabel',[],'FontSize',14);
  hold on;box on;
  clear n1 xout;
  %mat(istn).res_pw2 = (mat(istn).res_pw-mean(mat(istn).res_pw))/std(mat(istn).res_pw);
  mat(istn).res_pw2 = (mat(istn).res_pw-mean(mat(istn).res_pw));
  [n1,xout] = hist(mat(istn).res_pw2,x);
  n1 = n1/trapz(xout,n1);
  [nx,ny]=stairs(xout,n1,'k');
  patch(nx,ny,[0.8,0.8,0.8]);
  stairs(xout,n1,'k');
  plot(xx,nd,'-k','LineWidth',1)
  axis([-6 6 0 0.6]);

  sd_check_ns(istn) = std(res_ns(istart(istn):iend(istn)));
  sd_check_pw(istn) = std(mat(istn).res_pw);


  jstn = jstn + 2;
end;
fact = sd_check_ns./sd_check_pw
sd_check_ns
sd_check_pw

save covmat.mat F mat;
%
%  Saves both, inverse and orgiginal Cov mat.
%
if(isave == 1)
  save(covmatfile2,'F');
  icovmat = inv(F(1).Cd);
  save(covmatfile,'icovmat','-ascii');
  sdtmp = mat(1).sd_nonstat;
  save(sdfile,'sdtmp','-ascii');
  for istn = 2:NSTN
    clear icovmat;
    icovmat = inv(F(istn).Cd);
    save(covmatfile,'icovmat','-ascii','-append');
    clear sdtmp;
    sdtmp = mat(istn).sd_nonstat;
    save(sdfile,'sdtmp','-ascii','-append');
  end;
  for istn = 1:NSTN
    covmat = F(istn).Cd;
    save(covmatfile,'covmat','-ascii','-append');
  end;

end;
return;
