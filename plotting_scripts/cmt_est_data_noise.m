function [] = cmt_est_data_noise(filebase);
%%
%% Estimate non-stationary residual covariance matrix for 
%% W phase data to be used in CMT inversion.
%%

%% Switch for removing outliers:
irout = 1;
%% Switch to save Cd into file
isave = 1;

noisefile   = strcat(filebase,'_filt_rec.mat');
noiseoutfile   = strcat(filebase,'_noise.mat');
%noisefile   = strcat(filebase,'_filt_rec.dat');
polyfile   = strcat(filebase,'_poly_fit.dat');
covmatfile  = strcat(filebase,'_icovmat.dat');
sdfile      = strcat(filebase,'_stat_sd.dat');
covmatfile2 = strcat(filebase,'_covmat.mat');

%% For Maule remove tome window with event:
%%  201002262031A RYUKYU ISLANDS, JAPAN
%%  Date: 2010/ 2/26   Centroid Time: 20:31:29.7 GMT
%%  Lat=  25.86  Lon= 128.61
%%  Depth= 18.0   Half duration= 7.4
%%  Centroid time minus hypocenter time:  2.7
%%  Moment Tensor: Expo=26  -0.754 0.413 0.342 0.795 0.800 3.360 
%%  Mw = 7.0    mb = 6.7    Ms = 7.0   Scalar Moment = 3.49e+26
%%  Fault plane:  strike=91    dip=80   slip=9
%%  Fault plane:  strike=360    dip=81   slip=170
t1_cut = 36000;
t2_cut = 44500;

parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');

deltt = 4;
[IMAP,ICOV,NVMX,NPV,ILOC,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay]=cmt_read_parfile(parfile);
[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

dat = dobs';

%% 10 hours of data noise at 4-s sampling:
load(noisefile);
%noise(:,t1_cut:t2_cut) = [];
noise = noise(:,1:deltt:end);
time = [1:size(noise,2)]*deltt;
lags = [fliplr(-time(1:end-1)),0,time(1:end-1)];

%poly=load(polyfile);
%poly(:,t1_cut:t2_cut) = [];
%poly = poly(:,1:deltt:end);

nx = 8;
ny = 6;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.04 0.02 0.08];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

jstn = 1;
fig_n(1) = figure();
fig_Ann(1) = figure();
fig_Ann2(1) = figure();
ifig = 1;
for istn = 1:NSTN;
  if(rem(istn,10)==0);
    fprintf(1,'\n Station No. %4i\n\n',istn);
  end;
  %% Compute autocovariance:
  Ann(istn,:) = xcov(noise(istn,:),'biased');
  sd(istn)    = std(noise(istn,:));

  %% Find outliers:
  n(istn).noise = noise(istn,:);
  n(istn).time = [1:size(n(istn).noise,2)]*deltt;
  if(irout == 1);
  for iter = 1:6;
    n(istn).sd  = std(n(istn).noise);
    idx = find(abs(n(istn).noise)>5.*n(istn).sd);
    idx1 = idx-10;
    idx2 = idx+10;
    clear idx;
    idx = 1;
    for iel = 1:length(idx1);
      idx = [idx,[idx1(iel):idx2(iel)]];
    end;
    idx = unique(idx);
    idx(find(idx < 2)) = [];
    idx(find(idx > length(n(istn).noise))) = [];
    n(istn).noise(idx)=[];
    n(istn).time = [1:size(n(istn).noise,2)]*deltt;
    clear idx idx1 idx2;
  end
  end
  P=polyfit(n(istn).time,n(istn).noise,2);
  n(istn).noise = n(istn).noise - polyval(P,n(istn).time);
  clear P;

  n(istn).lags = [fliplr(-n(istn).time(1:end-1)),0,n(istn).time(1:end-1)];

  n(istn).Ann = xcov(n(istn).noise,'biased');
  n(istn).sd  = std(n(istn).noise);
  sd2(istn)   =n(istn).sd;

  %% Cut out center according to W-phase length at that station:
  ilag0(istn) = floor(length(n(istn).lags)/2)+1;
  ilag1(istn) = floor(length(n(istn).lags)/2)-NTSMP(istn);
  ilag2(istn) = floor(length(n(istn).lags)/2)+NTSMP(istn);
  %disp([length([ilag0(istn):ilag2(istn)]),NTSMP(istn)])

  %% Damping with double expoential to help positive definiteness 
  %% (all eigenvalues positive, chol. decomp. exixst, ensures invertability)
  %% Increase damping until positive definite matrix achieved.
  damp_power = 0.001;
  for idamp = 1:1000;
    damp_power = damp_power * 1.05;
    clear damp;
    %damp = cos(pi*[-NTSMP(istn)+1:NTSMP(istn)-1]./(2*NTSMP(istn)-1)).^damp_power;
    damp = exp(-abs([-NTSMP(istn)+1:NTSMP(istn)-1])*damp_power);
    n(istn).AnnW = n(istn).Ann(ilag1(istn)+1:ilag2(istn)-1) .* damp;
    n(istn).Cd = zeros(NTSMP(istn),NTSMP(istn));
    for irow = 1:NTSMP(istn)
      n(istn).Cd(irow,:) = n(istn).AnnW(NTSMP(istn)-(irow-1):end-(irow-1));
    end;
    %% Test if matrix is postive definite with Cholesky decomposition:
    [~,ipos] = chol(n(istn).Cd);
    if(ipos == 0);
      fprintf(1,'Converged after %4i steps with damping_power of %10.4f\n',idamp,damp_power);
      break;
    end;
  end;
  if(ipos > 0);
    [ipos,idamp,damp_power]
    fprintf(1,'Not converged after %i steps!\n',idamp);
    %error('ERROR: not converged');
    %return;
  end;

  %stop
  if(istn == (nx*ny)+1);
    fig_n(2) = figure();hold on;box on;
    fig_Ann(2) = figure();hold on;box on;
    fig_Ann2(2) = figure();hold on;box on;
    jstn = 1;
    ifig = ifig + 1;
  elseif(istn == (nx*ny)*2+1);
    fig_n(3) = figure();hold on;box on;
    fig_Ann(3) = figure();hold on;box on;
    fig_Ann2(3) = figure();hold on;box on;
    jstn = 1;
  end;

  figure(fig_Ann(ifig));
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(n(istn).lags,n(istn).Ann,'.k');
  set(gca,'XLim',[n(istn).lags(ilag1(istn)),n(istn).lags(ilag2(istn))]);

  figure(fig_n(ifig));
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  plot(n(istn).time,n(istn).noise,'-k');


  nsets = floor(length(n(istn).noise)/NTSMP(istn));
  n(istn).noise2 = zeros(nsets,NTSMP(istn));
  n(istn).time2 = [1:NTSMP(istn)]*deltt;
  n(istn).lags2 = [fliplr(-n(istn).time2(1:end-1)),0,n(istn).time2(1:end-1)];
  figure(fig_Ann2(ifig));
  subplot('Position',[loc(1,jstn) loc(2,jstn) spw sph]);hold on;box on;
  for iset = 1:nsets;
    n(istn).noise2(iset,:) = n(istn).noise((iset-1)*NTSMP(istn)+1:iset*NTSMP(istn));
    n(istn).Ann2 = xcov(n(istn).noise2(iset,:),'biased');
    plot(n(istn).lags2,n(istn).Ann2,'-k');
  end
  n(istn).Cd_samp = cov(n(istn).noise2);


  jstn = jstn + 1;
end;
figure();hold on;
plot(sd,'b');plot(sd2,'r');

%
%  Saves both, inverse and orgiginal Cov mat.
%
if(isave == 1)
  save(noiseoutfile,'n');
  icovmat = inv(n(1).Cd);
  save(covmatfile,'icovmat','-ascii');
  sdtmp = n(1).sd;
  save(sdfile,'sdtmp','-ascii');
  for istn = 2:NSTN
    clear icovmat;
    icovmat = inv(n(istn).Cd);
    save(covmatfile,'icovmat','-ascii','-append');
    clear sdtmp;
    sdtmp = n(istn).sd;
    save(sdfile,'sdtmp','-ascii','-append');
  end;
  for istn = 1:NSTN
    covmat = n(istn).Cd;
    save(covmatfile,'covmat','-ascii','-append');
  end;
end;
return;
