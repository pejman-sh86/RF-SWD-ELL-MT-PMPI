%function [] = rf_est_covmat(filebase);
filebase = 'HON';
include_rv = 0;
include_swd = 0;
include_ell = 0;
include_mt = 1;
%%
%% Estimate non-stationary residual covariance matrix for 
%% W phase data to be used in CMT inversion.
%%

%% Fraction of data set over which to est non stat std dev:
nfrac = 7.;
MAX_NAVE = 40;
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
%%
repfile     = strcat(filebase,'_replica.dat');
dobsfile    = strcat(filebase,'_dobs.dat');
covmatfile  = strcat(filebase,'_icovmat.dat');
sdfile      = strcat(filebase,'_nonstat_sd.dat');
covmatfile2 = strcat(filebase,'_covmat.mat');

parfile  = strcat(filebase,'_parameter.dat');
%M0file   = strcat(filebase,'_M0.dat');
%datfile  = strcat(filebase,'.hdf5');
% [IMAP,ICOV,NVMX,NPV,ILOC,IAR,IEXCHANGE,...
%  NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
%  minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
%  mindelay,maxdelay]=cmt_read_parfile(parfile);
% [NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

[IMAP,ICOV,iar,i_varpar,irv,itr,iswd, iell, imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL, NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVs,dVpVs, sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile(parfile);

NRF = ntr;
NRF
%%
if include_rv == 0; NRF = 0; end % be sure when include_rv == -1, NRF ~= 0
if include_swd == 0; NMODE = 0; end % be sure when include_swd == 1, NMODE ~= 0
if include_ell == 0; NMODE_ELL = 0; end % be sure when include_ell == 1, NMODE_ELL ~= 0
%%
%NSTN is the number of indipendent dataset such that their logliklihoods
%add. for example: NSTN = number of RF stations*NRF(for each station) +
%number of SWD stations*NMODE(for each stations) + number of MT stations
NSTN = 0;
if include_rv==-1; NSTN = NSTN + NRF; end
if include_swd==1; NSTN = NSTN + NMODE; end
if include_ell==1; NSTN = NSTN + NMODE_ELL; end
if include_mt==1;  NSTN = NSTN + 1; end

NTSMP = zeros(NSTN);
if include_rv==-1; NTSMP(1:NRF) = NTIME; end
if include_swd==1; NTSMP(NRF+1:NRF+NMODE) = NDAT_SWD; end
if include_ell==1; NTSMP(NRF+NMODE+1:NRF+NMODE+NMODE_ELL) = NDAT_ELL; end
if include_mt==1;  NTSMP(NRF+NMODE+NMODE_ELL+1) = NDAT_MT; end
%%
pred       = importdata(repfile);   %row vector

%dat        = dobs';
dat        = importdata(dobsfile);  %row vedctorres

%%
if include_mt == 1
    istart0mt = NRF*NTIME + NMODE*NDAT_SWD +NMODE_ELL*NDAT_ELL;

    real_predZ = pred(istart0mt+1: istart0mt+NDAT_MT);
    imag_predZ = pred(istart0mt+NDAT_MT+1: istart0mt+2*NDAT_MT);
    real_datZ = dat(istart0mt+1: istart0mt+NDAT_MT);
    imag_datZ = dat(istart0mt+NDAT_MT+1: istart0mt+2*NDAT_MT);

    pred(istart0mt+1: istart0mt+NDAT_MT) = complex(real_predZ, imag_predZ);
    pred(istart0mt+NDAT_MT+1: istart0mt+2*NDAT_MT) = [];
    dat(istart0mt+1: istart0mt+NDAT_MT) = complex(real_datZ, imag_datZ);
    dat(istart0mt+NDAT_MT+1: istart0mt+2*NDAT_MT) = [];

end

%% Residual errors:
%res = dat-pred(3,:);
res = dat - pred;
%std(res)

nx = 6;
ny = 6;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

if(ISIM == 1);
sd_dat = 0.000107944073408472;
ressim = sd_dat*randn(size(dat));
fig_cov1 = figure();hold on;box on;
fig_cov2 = figure();hold on;box on;
fig_cov3 = figure();hold on;box on;
fig_sd1 = figure();hold on;box on;
fig_sd2 = figure();hold on;box on;
fig_sd3 = figure();hold on;box on;

jjstn = 0;
jjjstn = 0;
for istn = 1:NSTN;

  jjstn = jjstn + 1;
  jjjstn = jjjstn + 1;
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

  figure(fig_cov1)
  if(istn == (nx*ny)+1);
    figure(fig_cov2);hold on;box on;
    jjstn = 1;
  elseif(istn == (nx*ny)*2+1);
    figure(fig_cov3);hold on;box on;
    jjstn = 1;
  end;
  subplot('Position',[loc(1,jjstn) loc(2,jjstn) spw sph]);hold on;box on;
  imagesc(Csim)
  set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);
  set(gca,'YLim',[1 NTSMP(istn)],'YDir','reverse');

  figure(fig_sd1)
  if(istn == (nx*ny)+1);
    figure(fig_sd2);hold on;box on;
    jjjstn = 1;
  elseif(istn == (nx*ny)*2+1);
    figure(fig_sd3);hold on;box on;
    jjjstn = 1;
  end;
  subplot('Position',[loc(1,jjjstn) loc(2,jjjstn) spw sph]);hold on;box on;
  plot(sig(istart(istn):iend(istn)));
  set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);

end;end;

fig_sdnon1 = figure();hold on;box on;
jstn = 1;
%%
%% Estimate non-stationary standard deviation:
%%
for istn = 1:NSTN;
  if(istn == (nx*ny)+1);
    fig_sdnon2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)*2+1);
    fig_sdnon3 = figure();hold on;box on;
    jstn = 1;
  end;
  NAVE = round(double(NTSMP(istn))/nfrac);
  if(NAVE > MAX_NAVE);NAVE = MAX_NAVE;end;
  iend(istn) = sum(NTSMP(1:istn));
  istart(istn) = iend(istn)-NTSMP(istn)+1;
  %% Residuals mean removed
  res_mr(istart(istn):iend(istn)) = res(istart(istn):iend(istn));%-mean(res(istart(istn):iend(istn)));
  mat(istn).res_mr = res_mr(istart(istn):iend(istn));
  sd_stat(istn) = std(res_mr(istart(istn):iend(istn)),1);
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

      %sd_nonstat(idat) = sqrt(mean(res_mr(NAVE1:NAVE2).^2.));
      sd_nonstat(idat) = sqrt( mean( abs( res_mr(NAVE1:NAVE2) ).^2. ) );
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
fig_res1 = figure();
for istn = 1:NSTN;
  if(istn == (nx*ny)+1);
    fig_res2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)*2+1);
    fig_res3 = figure();hold on;box on;
    jstn = 1;
  end;
  %% Residuals standardized with non stat std dev:
  res_ns(istart(istn):iend(istn)) = res(istart(istn):iend(istn))./sd_nonstat(istart(istn):iend(istn));
  sd_resns(istn) = std(res_ns(istart(istn):iend(istn)));
  %% Residuals standardized with stat std dev:
  res_s(istart(istn):iend(istn)) = res(istart(istn):iend(istn))/sd_stat(istn);
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
%   plot(res_ns(istart(istn):iend(istn)),'-b');
%   plot(res_s(istart(istn):iend(istn)),'--r');
  plot( real(res_ns(istart(istn):iend(istn))),'-b' );
  plot( real(res_s(istart(istn):iend(istn))),'--r' );
  set(gca,'XLim',[1 NTSMP(istn)],'FontSize',14,'LineWidth',1);
  set(gca,'YLim',[-3 3]);
  if include_mt == 1
      subplot('Position',[loc(1,jstn)+2.*spw loc(2,jstn) spw sph]);hold on;box on;
      plot( imag(res_ns(istart(istn):iend(istn))),'-b' );
      plot( imag(res_s(istart(istn):iend(istn))),'--r' );
      set(gca,'XLim',[1 NTSMP(istn)],'FontSize',14,'LineWidth',1);
      set(gca,'YLim',[-3 3]);
  end
  jstn = jstn + 1;
end;
%disp('Standard deviation:');
%disp(sd_resns);
%disp('Variance:');
%disp(sd_resns.^2.);

%%
%% Compute autocovariance:
%%
fig_Axx1 = figure();hold on;box on;
jstn = 1;
for istn = 1:NSTN;
  if(istn == (nx*ny)+1);
    fig_Axx2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)*2+1);
    fig_Axx3 = figure();hold on;box on;
    jstn = 1;
  end;
  if(iunbiased == 0);
    %[mat(istn).Axx,mat(istn).lags] = xcov(res_ns(istart(istn):iend(istn)));
    [mat(istn).Axx,mat(istn).lags] = xcov( conj( res_ns(istart(istn):iend(istn)) ) );
    mat(istn).Axxn= mat(istn).Axx/double(NTSMP(istn));
  else;
    %[mat(istn).Axx,mat(istn).lags] = xcov(res_ns(istart(istn):iend(istn)),'unbiased');
    [mat(istn).Axx,mat(istn).lags] = xcov( conj( res_ns(istart(istn):iend(istn)) ),'unbiased' );
    mat(istn).Axxn= mat(istn).Axx;
  end;
 %mat(istn).Axxu= xcov(res_ns(istart(istn):iend(istn)),'unbiased');
  mat(istn).Axxu= xcov( conj( res_ns(istart(istn):iend(istn)) ),'unbiased' );

  %% Damping to help positive definiteness (all eigenvalues positive, chol. decomp. exixst)
  clear damp;
  damp = cos(pi*[-NTSMP(istn)+1:NTSMP(istn)-1]./(2*NTSMP(istn)-1)).^damp_power;
  %gc0 = figure(1);plot(x);
  mat(istn).Axxn = mat(istn).Axxn .* damp(1:end);
  mat(istn).Axxu = mat(istn).Axxu .* damp(1:end);

  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(mat(istn).lags,real(mat(istn).Axxn));
  plot(mat(istn).lags,real(mat(istn).Axxu),'--r');
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'FontSize',14,'LineWidth',1);
  set(gca,'YLim',[-.6 1.2]);
  if include_mt == 1
      subplot('Position',[loc(1,jstn)+2*spw loc(2,jstn) spw sph]);hold on;box on;
      plot(mat(istn).lags,imag(mat(istn).Axxn));
      plot(mat(istn).lags,imag(mat(istn).Axxu),'--r');
      set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'FontSize',14,'LineWidth',1);
      set(gca,'YLim',[-.6 1.2]); 
  end
  jstn = jstn + 1;
end;
legend('biased','unbiased')

%%
%% Set up covariance matrix in a diagonal form:
%%
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
fig_Cd1 = figure();hold on;box on;
jstn = 1;
for istn = 1:NSTN;
  if(istn == (nx*ny)+1);
    fig_Cd2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)*2+1);
    fig_Cd3 = figure();hold on;box on;
    jstn = 1;
  end;
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

  %clear L C r;
  %L = mat(istn).L;
  %R = mat(istn).R;
  %C = mat(istn).Cd;
  %r = mat(istn).res_mr';

  %x1a= (inv(L)*r);
  %x1b= (inv(L)*r).';
  %x2a= (r*r.');
  %X3 = inv(L)*(L*L.')*inv(L.');
  %X4 = inv(L)*C*inv(L.');


  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  imagesc(real(mat(istn).Cd));
  set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);
  set(gca,'YLim',[1 NTSMP(istn)],'YDir','reverse');
  set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  if include_mt == 1
      subplot('Position',[loc(1,jstn)+2.*spw loc(2,jstn) spw sph]);hold on;box on;
      imagesc(imag(mat(istn).Cd));
      set(gca,'XLim',[1 NTSMP(istn)],'LineWidth',1);
      set(gca,'YLim',[1 NTSMP(istn)],'YDir','reverse');
      set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14)   
  end
  jstn = jstn + 1;
end;

%%
%% Check residuals by pre-whitening:
%%
nx = 16;
ny = 8;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
jstn = 1;
fig_Axx_check1 = figure();hold on; box on;
for istn = 1:NSTN;
  if(istn == (nx*ny)/2+1);
    fig_Axx_check2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)+1);
    fig_Axx_check3 = figure();hold on;box on;
    jstn = 1;
  end;

  %disp('max L*L^-1')
  %max(max(mat(istn).L*inv(mat(istn).L)))
  %max(max(mat(istn).Cd-mat(istn).L*mat(istn).L'))


  mat(istn).res_pw = inv(mat(istn).L)*res_mr(istart(istn):iend(istn))';
  [mat(istn).Axx_pw,mat(istn).lags_pw] = xcov(mat(istn).res_pw);
  mat(istn).Axx_pw = mat(istn).Axx_pw/NTSMP(istn);

  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(mat(istn).lags,real(mat(istn).Axxn)/(max(real(mat(istn).Axxn))),'-k.')
  plot([mat(istn).lags(1) mat(istn).lags(end)],[0 0],'--k');
  set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);

  subplot('Position',[loc(1,jstn+1) loc(2,jstn+1) spw sph]);hold on;box on;
  plot(mat(istn).lags_pw,real(mat(istn).Axx_pw)/(max(real(mat(istn).Axx_pw))),'-k.')
  plot([mat(istn).lags_pw(1) mat(istn).lags_pw(end)],[0 0],'--k');
  set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
  set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);
  
  if include_mt == 1
      subplot('Position',[loc(1,jstn)+4*spw loc(2,jstn) spw sph]);hold on;box on;
      plot(mat(istn).lags,imag(mat(istn).Axxn)/(max(imag(mat(istn).Axxn))),'-b.')
      plot([mat(istn).lags(1) mat(istn).lags(end)],[0 0],'--b');
      set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
      set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);

      subplot('Position',[loc(1,jstn+1)+4*spw loc(2,jstn+1) spw sph]);hold on;box on;
      plot(mat(istn).lags_pw,imag(mat(istn).Axx_pw)/(max(imag(mat(istn).Axx_pw))),'-b.')
      plot([mat(istn).lags_pw(1) mat(istn).lags_pw(end)],[0 0],'--b');
      set(gca,'YTickLabel',[],'XTickLabel',[],'FontSize',14);
      set(gca,'XLim',[-NTSMP(istn) NTSMP(istn)],'YLim',[-.6 1.1]);
  end

  jstn = jstn + 2;
end;


jstn = 1;
fig_hist_check1 = figure();hold on; box on;
x  = [-6.375:.75:6.375]';
xx = [-6.25:.01:6.25]';
nd = 1./sqrt(2*pi)*exp(-(xx.^2)/2.);
for istn = 1:NSTN;

  if(istn == (nx*ny)/2+1);
    fig_hist_check2 = figure();hold on;box on;
    jstn = 1;
  elseif(istn == (nx*ny)+1);
    fig_hist_check3 = figure();hold on;box on;
    jstn = 1;
  end;
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  clear n1 xout;
  [n1,xout] = hist(res_ns(istart(istn):iend(istn)),x);
  n1 = n1/trapz(xout,n1);
  [nnx,nny]=stairs(xout,n1,'k');
  patch(nnx,nny,[0.8,0.8,0.8]);
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
  [nnx,nny]=stairs(xout,n1,'k');
  patch(nnx,nny,[0.8,0.8,0.8]);
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
istn0mt = NRF + NMODE + NMODE_ELL;
if(isave == 1)
  save(covmatfile2,'F');
  icovmat = inv(F(1).Cd);
  if 1 > istn0mt
      real_icovmat = real(icovmat);
      imag_icovmat = imag(icovmat);
      icovmat = [real_icovmat, imag_icovmat];
  end
  save(covmatfile,'icovmat','-ascii');
  sdtmp = mat(1).sd_nonstat;
  save(sdfile,'sdtmp','-ascii');
  for istn = 2:NSTN
    clear icovmat;
    icovmat = inv(F(istn).Cd);
    if istn > istn0mt
        real_icovmat = real(icovmat);
        imag_icovmat = imag(icovmat);
        icovmat = [real_icovmat, imag_icovmat];
    end
    save(covmatfile,'icovmat','-ascii','-append'); 
    clear sdtmp;
    sdtmp = mat(istn).sd_nonstat;
    save(sdfile,'sdtmp','-ascii','-append');
  end;
   for istn = 1:NSTN
   covmat = F(istn).Cd;
   if istn > istn0mt
       real_covmat = real(covmat);
       imag_covmat = imag(covmat);
       covmat = [real_covmat; imag_covmat];
  end
      save(covmatfile,'covmat','-ascii','-append');
   end;

end;
%return;
