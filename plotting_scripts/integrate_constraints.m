function []=integrate_constraints();

%% Site 1 Malta:
%minlim = [0.150, log(1.e6), log(0.001), log(4.53e-5),log(1.e6) ];
%maxlim = [0.910, log(5.e8), log(0.147), log(0.4),    log(5.e8) ];
%% Site 2 Malta:
minlim = [0.300, log(1.e6), log(0.001), log(3.35e-4),log(1.e6) ];
maxlim = [0.910, log(5.e8), log(0.147), log(0.4),    log(5.e8) ];

%% Site 3 Tyrrhenian:
%minlim = [0.100, log(1.e6),  log(0.001), log(8.94e-4), log(1.e6)]
%maxlim = [0.910, log(1.e11), log(0.147), log(0.4),     log(1.e11)]

maxpert = maxlim-minlim;
mstart= [0.85   17.00   0.040  -4.00   14.0];
minlimmisc = [3.55e10, 2.3526e9,2710.,1028.];
maxlimmisc = [3.65e10, 2.3586e9,2740.,1031.];
maxpertmisc = maxlimmisc-minlimmisc;
miscstart = [3.6400000e+10   2.3530000e+09   2.7200000e+03   1.0290000e+03];

minlimcr = [1450. 1.2    0. 0.  0.];
maxlimcr = [6000. 2.6 2500. 2. 20.];

NSMP = 5e5;
NAVEF = 6;
%bands = [400 504  635  800 1008 1270 1600 2016 2540 3200];
%fr = [400 504  635  800 1008 1270 1600 2016 2540 3200];
%fr = [400 504  635  800 1008 1270 1600];
%% Site 2:
fr = [300., 400., 504., 635., 800.,1008.,1270.,1600.,2000.,2500.,3150.];

sqrt2 = sqrt(2.);

%for iband=1:length(bands);
% flo = bands(iband) / (2.^(1./30.)); % Use 1/3 octave
% fhi = bands(iband) * (2.^(1./30.)); %
% fstep = (fhi-flo)/(NAVEF-1);
% fr((iband-1)*NAVEF+1:iband*NAVEF) = flo + ([1:NAVEF]-1) .* fstep;
%end
NFREQ = length(fr);

mkeep = zeros(NSMP,5);
ekeep = zeros(NSMP,5);
ikeep = 1;
for ismp=1:NSMP;
  
  if(rem(ismp,100000)==0);disp([ismp,ikeep,ikeep/ismp]);end;
  ranu = rand(1,9);
  m = minlim + ranu(1:5) .* maxpert;
  misc = minlimmisc + ranu(6:9) .* maxpertmisc;

  ireject = 0;
  %% Enforce Ks<Kp:
  if(1.5*m(2)<m(5));ireject = 1;end;
  if(ireject == 0);
    [cp,alfp,cs,alfs,rho]=VGSlambda_no_tau_s(misc(1),misc(2),misc(3),misc(4),...
    m(1),exp(m(2)),exp(m(5)),exp(m(3)),exp(m(4)),fr);
    %[cp,alfp,cs,alfs,rho]=VGSlambda_no_tau_s(misc(1),misc(2),misc(3),misc(4),...
    %m(1),m(2),m(5),exp(m(3)),exp(m(4)),fr);

    rho = rho/1000.;
    if(ireject == 0);
      if(rho<minlimcr(2) | rho>maxlimcr(2));ireject = 1;end;
    end;
    %% 
    %% Check fro rejection...
    %% 
    if(ireject == 0);
      cl=(1.54-0.907*rho+0.3695*rho.^1.88)*1.5004*1000.;
      ch=(1.62-0.907*rho+0.3695*rho.^2.05)*1.5014*1000.;
      if(cl<minlimcr(1)); cl=minlimcr(1);end;
      if(ch>maxlimcr(1)); ch=maxlimcr(1);end;
      if(rho>2.000); ch=maxlimcr(1);end;
      for ifr = 1:NFREQ;
        %% Poisson:
        if(cp(ifr)<cs(ifr)*sqrt2);ireject = 1;end;
        %% Alpha_P, Alpha_S:
        %if(ireject == 0);
        %  if(alfp(ifr)>maxlimcr(4));ireject = 1;end;
        %end;
        %if(ireject == 0);
        %  if(alfs(ifr)>maxlimcr(5));ireject = 1;end;
        %end;
        %% Hamilton:
        if(ireject == 0);
          if(cp(ifr)>ch);ireject = 1;end;
        end;
        if(ireject == 0);
          if(cp(ifr)<cl);ireject = 1;end;
        end;
      end;
    end;
  end;
  if(ireject == 0);
    mkeep(ikeep,:) = m;
    misckeep(ikeep,:) = misc;
    ekeep(ikeep,:) = [mean(cp),rho,mean(cs),mean(alfp),mean(alfs)];
    ikeep = ikeep + 1;
  end;
end;
mkeep(ikeep:end,:)=[];
ekeep(ikeep:end,:)=[];
misckeep(ikeep:end,:)=[];
ikeep
save tmp.mat mkeep ekeep misckeep ikeep
figure();hold on;box on;
for ipar = 1:5;
  subplot(1,5,ipar);hold on;
  hist(mkeep(:,ipar),400);
  set(gca,'XLim',[minlim(ipar) maxlim(ipar)]);
end;
figure();hold on;box on;
for ipar = 1:4;
  subplot(1,4,ipar);hold on;
  hist(misckeep(:,ipar),400);
  set(gca,'XLim',[minlimmisc(ipar) maxlimmisc(ipar)]);
end;
figure();hold on;box on;
for ipar = 1:5;
  subplot(1,5,ipar);hold on;
  hist(ekeep(:,ipar),400);
  %set(gca,'XLim',[minlimcr(ipar) maxlimcr(ipar)]);
end;
return;
