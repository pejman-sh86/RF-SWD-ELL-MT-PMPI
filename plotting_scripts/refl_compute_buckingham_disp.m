%function [fr, z, nfc, nfa, nfr, meac, mear, meaa] = compute_buckingham_disp(samplefile);
%samplefile = 'p0006_pave010_0988_2513_sample.mat';
%samplefile = 'x_s01_1_400_2000_vgs_sample.mat';
%samplefile = 'blde4_b_bl_ping_id2_nb2_15deg_400_1600_sample.mat';
samplefile = 'blge5_3_bl_ping_id1_7oct_verbose_sample.mat';
%samplefile = 'blge5_3_bl_ping_id1_7oct_verbose_05Ang_sample.mat';
%samplefile = 'blge5_3_bl_ping_id1_IvgsEven_5Fres1d5_sample.mat';

set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg = 1;
idispmarg = 1;
ipltdisp  = 1;
ivoro = 0;
ilogalf = 1;
inorm = 0;
i_vref = 1; %% Assumes reference model is employed in inversion
isyn  = 0;
ipltcs= 0;  %% This is to only plot shear profiles
ishear= 0;  %% This is for 5 panel plot with shear and comp
isave = 1;
idatfit = 0;
ibucking = 3; % 1 = MVF
              % 2 = GS
              % 3 = VGS
ivgsnosh= 0;  % For VGS, strain hard as misc parameter?

filebase  = strrep(samplefile,'sample.mat','');
datafile  = strrep(samplefile,'_sample.mat','.txt');
vreffile  = strrep(samplefile,'_sample.mat','_vel_ref.txt');
if(ipltcs == 0);
  plotfile2 = strcat(filebase,'transdim_disp_prof_cp.');
  plotfile3 = strcat(filebase,'transdim_disp_curves_cp.');
else;
  plotfile2 = strcat(filebase,'transdim_disp_prof_cs.');
  plotfile3 = strcat(filebase,'transdim_disp_curves_cs.');
end;
plotext1  = 'fig';
plotext2  = 'png';
plotext3  = 'eps';
IEPS = 0;

NLMX = 20;
NAVEF = 6;
IGA = 0;
if(ibucking == 3 & ivgsnosh == 0);
  NPL  = 5;
  if(ishear == 1);
    NPL  = 6;
  end;
else;
  NPL = 4;
end;
NPLR = 4;
%if(ishear == 1 | ipltcs==1);
%  NPLR  = 6;
%end;
frbw = 1.0/15.0;
FBW = 50.0;
%% AUV
%bands = [988., 1113., 1288., 1913., 2263., 2513.];

%% MUD
%bands = [400. 460. 528. 607. 697. 801. 920. 1057. 1214.];

%% Site 1  8 m packet
%bands = [400 504  635  800 1008 1270];
%bands = [400., 504.,  635.,  800., 1008., 1270., 1600., 2016.];

%% Site 02  4 m packet
%bands = [300.,400.,504.,635.,800.,1008.,1270,1600];
%bands = [300.,400.,504.,635.,800.,1008.,1270.,1600.,2000.,2500.,3150.];
%bands = [315.,400.,500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.,220000.];

%% Site 02 13 m Malta
%bands = [400.,600.,800.,1000.,1200.,1400.,1600.];

%% Clutter Site:
%bands = [1601.4, 2017.6, 2542, 3202.8, 4035.2, 5084.1, 6405.5];
%bands = [1677.4, 2017.6, 2542, 3202.8, 4035.2, 5084.1];
%bands = [2017.6, 2542, 3202.8, 4035.2, 5084.1, 5800.];
%bands = [2212.96, 2542.02, 3202.75, 4035.21, 5084.05, 5840.04];

%% Tyrrhenian Sea Site 3:
%bands = [400 504  635  800 1008 1270 1600];
%bands = [400 504  635  800 1008 1270 1600 2016 2540 3200];
%%  6m packet
%bands = [504  635  800 1008 1270 1600 2016];

%% Tyrrhenian Sea Site 5; 14.8 m:
%bands = [464.  566.  690.  929.  1132.  1250.];
bands = [464. 690. 928. 1250. 1682.];

%% For plotting labels:
%bands2 = [ 315.,400., 504., 800.,1270.,2000.,3150.,4035.];
%bands2 = [ 400., 504., 800.,1270.,1600.];
%bands2 = [400. 460. 528. 607. 697. 801. 920. 1057. 1214.];
%bands2 = [988., 1113., 1288., 1913., 2263., 2513.];
%bands2 = [400., 600.,  800., 1000., 1200., 1400., 1600.];
bands2 = [464.  566.  690.  929.  1132.  1250.];

if(isyn==0);
  icore = 1;
  %% Tyrrhenian Site 05
  zplt = [2., 4., 6.5, 9., 13., 15.];
else;
  icore = 0;
  zplt = [0.6, 2.0, 3.2, 4.8, 6.0];
end;
nzplt = length(zplt);
NBAND = length(bands);
NFREQ = NBAND*NAVEF;

%% Shear VGS simulations:
if(ibucking == 3);
  ktru = 4;
  mtru = [0.15, 0.800, 18.000, -6.6036, -4.173,...
          0.60, 0.700, 19.000, -6.6036, -1.037,...
          1.40, 0.600, 20.000, -6.6036, -6.203,...
          4.445,0.500, 21.000, -6.6036, -5.210,...
                0.400, 22.000, -3.9380, -5.210];
  mtrumisc = [35612723556.66, 2355895918.32, 2724.76, 1029.88];
  %mtru = [0.15, 0.800, 18.000, -6.6036, -4.173, 14.000,...
  %        0.60, 0.700, 19.000, -6.6036, -1.037, 14.000,...
  %        1.40, 0.600, 20.000, -6.6036, -6.203, 16.000,...
  %        4.445,0.500, 21.000, -6.6036, -5.210, 18.600,...
  %              0.400, 22.000, -3.9380, -5.210, 20.000];
  %mtrumisc = [35612723556.66, 2355895918.32, 2724.76, 1029.88];
  if(ivoro > 0);
    %% True parameters for simulation:
    vorotru(1,:) = [0.000, 0.800, 18.000, -6.6036, -4.173, 14.000];
    vorotru(2,:) = [0.290, 0.700, 19.000, -6.6036, -1.037, 14.000];
    vorotru(3,:) = [0.910, 0.600, 20.000, -6.6036, -6.203, 14.000];
    vorotru(4,:) = [1.890, 0.500, 21.000, -6.6036, -5.210, 18.600];
    vorotru(5,:) = [3.500, 0.400, 22.000, -3.9380, -5.210, 20.000];
    voroidxtru(1,:) = [1, 1, 1, 1, 1, 1];
    voroidxtru(2,:) = [1, 1, 1, 0, 1, 0];
    voroidxtru(3,:) = [1, 1, 1, 0, 0, 0];
    voroidxtru(4,:) = [1, 1, 1, 0, 1, 1];
    voroidxtru(5,:) = [1, 1, 1, 0, 0, 1];
     if(ivoro == 1);
       [mtru,nuniquetru] = refl_voro_to_lay(vorotru,voroidxtru,5,6,5);
     elseif(ivoro == 2);
       [mtru,nuniquetru] = refl_laynode_to_lay(vorotru,voroidxtru,5,6,5);
     end;
     % 0.1500      0.8000     18.0000     -6.6036     -4.1730     14.0000
     % 0.4500      0.7000     19.0000     -6.6036     -1.0370     14.0000
     % 0.8000      0.6000     20.0000     -6.6036     -6.2030     16.0000
     % 3.0450      0.5000     21.0000     -6.6036     -5.2100     18.6000
     %             0.4000     22.0000     -3.9380     -5.2100     20.0000
     mtrumisc = [36234377700., 2353069330., 2737.9595, 1029.39201];
     %nunique
  end;
elseif(ibucking == 0);
  ktru = 4;
  mtru = [0.15, 1474.1, 1.2841, 0.001283,   35.22, 1.569,...
          0.60, 1530.9, 1.7553, 0.048731,   45.88, 2.386,...
          1.40, 1774.5, 2.0146, 0.126840,  201.89, 2.6720,...
          6.40, 2120.2, 2.0451, 0.163191,  407.03, 1.9810,...
                3349.6, 2.3451, 0.146134, 1414.42, 0.30403];
end;

%% Site07:
%% Site07:
pmin = [1400 1.0 -4.]';
pmax = [2100 2.5  1.5]';
pminplt = [1450 1.2 -4.]';
pmaxplt = [2100 2.5 1.5]';
ylimcmin = [1450 1450 1550 1550 1500 1500]';
ylimcmax = [1550 1550 1650 1650 1650 1650]';
ylimamin = pmin(3)*ones(6,1);
ylimamax = pmax(3)*ones(6,1);
if(1 == 1);
    %% Syn var par
    %pmin = [1450 1.20 (0.0)]';
    %pmax = [2100 2.40 (1.0)]';
    %% Site 02:
    %pmin = [1450 1.15 log10(1e-6)]';
    %pmax = [1650 1.85 log10(1.0)]';
    pmin2 = [ 0  1.2 -2.5]';
    pmax2 = [75  2.1  1.0]';
    if(ishear == 1 | ipltcs==1);
      %% Syn var par
      %pmin = [1400 1.10 (0.0)    0.  0.]';
      %pmax = [2500 2.60 (2.0) 1500. 20.]';
      %% Site 02:
      %pmin = [1450 1.20 (0.0)    0.  0.]';
      %pmax = [1600 2.00 (1.0)  500. 20.]';
      pmin = [1450 1.20 log10(1e-4)    0. log10(1e-2)]';
      pmax = [1820 2.20 log10(2.0)   300. log10(1000.)]';
      %pmin = [1450 1.20 (1e-4)    0. (1e-2)]';
      %pmax = [1620 2.00 (0.5)   300. (20.)]';
      %% Tyrrhenian
      %pmin = [1450 1.20 log10(1e-4)    0.  log10(1e-1)]';
      %pmax = [2200 2.50 log10(2.0)   1400. log10(250.)]';
    end;
    ylimcmin = [1450 1450 1450 1600]';
    ylimcmax = [1550 1550 1600 1800]';
    ylimamin = [-3 -3 -3 -3];
    ylimamax = [0.5 0.5 0.5 .5];
    %ylimamin = pmin(3)*ones(nzplt,1);
    %ylimamax = pmax(3)*ones(nzplt,1);
end
hmax = 15.;
hmax2 = 15.;
if(idatfit == 1)
   tmp = dlmread(datafile);
   z_t    = tmp(1,1);
   cw     = tmp(2,1);
   rw     = tmp(3,1);
   hmax   = tmp(4,1)+.2;
   %hmax   = 18.;
   dobs   = tmp(5:length(bands)+4,:)';
   angobs = tmp(length(bands)+5,:);
   NANG   = length(angobs);
%   rex = tmp(length(bands)+6:end,:)';
end

load(samplefile);
%A = A(1:2:end,:);
NPROF = size(A,1);
if(NPROF>100000);
  BURNIN = 20000;
else;
  BURNIN = ceil(NPROF/8);
end;
A = A(BURNIN:end,:);
NPROF = size(A,1);
if(NPROF > 150000);
  idx = randperm(NPROF);
  idx = idx(1:150000);
  idx = sort(idx);
  A=A(idx,:);
  NPROF = size(A,1);
end;
k = A(:,4);
logL = A(:,1);
NFP = (k*NPL)+(NPL-1);
m = A(:,5:end-5);
NPROF  = size(m,1);
NPROF2 = round(NPROF/1);
if(NPROF2>30000);NPROF2 = 30000;end;
%disp('Sample size: '),disp(size(m));

for iband=1:NBAND;
  if(NAVEF > 1)
    if(IGA == 1)
      flo = bands(iband) - bands(iband)*frbw;
      fhi = bands(iband) + bands(iband)*frbw;
    elseif(IGA == 2)
      flo = bands(iband) / (2.^(1./15.)); % Use 1/3 octave
      fhi = bands(iband) * (2.^(1./15.)); %
    elseif(IGA == 3)
      flo = bands(iband) - bands(iband)*frbw;
      fhi = bands(iband) + bands(iband)*frbw;
    else
      flo = bands(iband) - FBW;
      fhi = bands(iband) + FBW;
    end
    fstep = (fhi-flo)/(NAVEF-1);
    fr((iband-1)*NAVEF+1:iband*NAVEF) = flo + ([1:NAVEF]-1) .* fstep;
  else
    fr = bands(iband);
  end
end;
if(isyn == 1);
  if(ibucking == 3);
    if(ivgsnosh == 0);
      if(ishear == 0);
        for ilay = 1:ktru+1;
          ipar = (ilay-1)*NPL+2;
          if(ilay == ktru+1);ipar = (ilay-1)*NPL+1;end;
          [cptru(ilay,:),alfptru(ilay,:),cstru(ilay,:),alfstru(ilay,:),...
           rhotru(ilay,1)] = ...
          VGSlambda_Mod(mtrumisc(1),mtrumisc(2),mtrumisc(3),mtrumisc(4),...
          mtru(ipar),10.^(mtru(ipar+1)),10.^(mtru(ipar+2)),10.^(mtru(ipar+3)),fr);
        end;
      end;
    end;
  end;
end;
%save('sim_VGS.mat','cptru','alfptru','cstru','alfstru','rhotru','mtrumisc','mtru','fr');
%disp('Computing VGS...');
cp   = zeros(NLMX+1,NFREQ,NPROF);
cs   = zeros(NLMX+1,NFREQ,NPROF);
rho  = zeros(NLMX+1,NPROF);
alfp = zeros(NLMX+1,NFREQ,NPROF);
alfs = zeros(NLMX+1,NFREQ,NPROF);
por  = zeros(NPROF,NLMX+1);
if(i_vref == 1);
    vel_ref = dlmread(vreffile);
    vel_ref(1,:) = [];
end;

for imod = 1:NPROF;
%  imod = 1
  %imod
  %k(imod)
  %tic;
  if(i_vref == 0);
    for ilay = 1:k(imod)+1;
      ipar = (ilay-1)*NPL+2;
      if(ilay == k(imod)+1);ipar = (ilay-1)*NPL+1;end;
      por(imod,ilay)= m(imod,ipar);
    end;
  else
    for ilay = 1:k(imod)+1;
      ipar = (ilay-1)*NPL+2;
      if(ilay == k(imod)+1);ipar = (ilay-1)*NPL+1;end;
      if(ilay <= k(imod));
        z = m(imod,ipar-1);
      else;
        z = m(imod,ipar-NPL);
      end;
      [pref]=refl_getrefpor(z,vel_ref);
      por(imod,ilay)= pref + m(imod,ipar);
    end;
  end;
  %toc
  for ilay = 1:k(imod)+1;
    ipar = (ilay-1)*NPL+2;
    if(ilay == k(imod)+1);ipar = (ilay-1)*NPL+1;end;
    [cp(ilay,:,imod),alfp(ilay,:,imod),...
     cs(ilay,:,imod),alfs(ilay,:,imod),rho(ilay,imod)] = ...
    VGSlambda_Mod(m(imod,end-4),m(imod,end-3),m(imod,end-2),...
    m(imod,end-1),por(imod,ilay),10.^(m(imod,ipar+1)),...
    10.^(m(imod,ipar+2)),10.^(m(imod,ipar+3)),fr);
    if(imag(cp(ilay,:,imod)) > 0.);
        stop
    end;
  end;
  %toc
end;
rho = rho/1000.;

figure();hold on;
for i=1:NFREQ;
  clear tmp;
  tmp = squeeze(cp(:,i,:));
  cp2 = tmp(:);
  rho2= rho(:);
  plot(rho2,cp2,'.k');
end;
%disp('...done computing VGS.');

%save model.mat model fr


NZI = 500;  % NZ for Interfaces
NZ  = 300;  % NZ for pcolor panels
nsmooth = ceil(NZ/80.);
NC = 500;
NC2 = 500;
NR = 500;
NA = 500
%NZI = 100;  % NZ for Interfaces
%NZ  = 100;  % NZ for pcolor panels
%nsmooth = ceil(NZ/80.);
%NC = 100;
%NC2 = 100;
%NR = 100;
%NA = 100;
if(imarg == 1)
  disp('Plot marginal profiles...');
  cp2   = zeros(NLMX+1,NPROF);
  alfp2 = zeros(NLMX+1,NPROF);

  ipltidx2 = 0;
  for ipltidx = 1:NAVEF:NFREQ+1;
    ipltidx2 = ipltidx2 + 1;
    tic;
    disp(ipltidx)
    if(ipltidx == NFREQ+1);ipltidx = ipltidx-1;end;
    cp2 = cp(:,ipltidx,:);
    alfp2 = alfp(:,ipltidx,:);

    plotfile1 = strcat(filebase,'f_',num2str(round(fr(ipltidx)),'%05i'),...
                       '_transdim_marg_cra_prof_cp.');
    if(isyn == 1);
      for ilay = 1:ktru+1;
        ipar1 = (ilay-1)*NPL+2;
        ipar2 = (ilay-1)*NPLR+2;
        if(ilay == ktru+1);
          ipar1 = (ilay-1)*NPL+1;
          ipar2 = (ilay-1)*NPLR+1;
        end;
        if(ilay < ktru+1);mtru2(ipar2-1)  = mtru(ipar1-1);end;
        %if(ipltcs==0)
          mtru2(ipar2)    = cptru(ilay,ipltidx);
          if(ilogalf == 0);
            mtru2(ipar2+2)  = alfptru(ilay,ipltidx);
          else;
            mtru2(ipar2+2)  = log10(alfptru(ilay,ipltidx));
          end;
          %if(ishear==1 | ipltcs==1)
          %  mtru2(ipar2+3)  = cstru(ilay,ipltidx);
          %  if(ilogalf == 0);
          %    mtru2(ipar2+4)  = alfstru(ilay,ipltidx);
          %  else;
          %    mtru2(ipar2+4)  = log10(alfstru(ilay,ipltidx));
          %  end
          %end
        %else
        %  mtru2(ipar2)    = cstru(ilay,ipltidx);
        %  mtru2(ipar2+2)  = alfstru(ilay,ipltidx);
        %end;
        mtru2(ipar2+1)  = rhotru(ilay)/1000.;
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %% COMPUTE PROFILE MARGINALS
    %%
    dc = (pmax(1)-pmin(1))/NC;
    dr = (pmax(2)-pmin(2))/NR;
    da = (pmax(3)-pmin(3))/NA;
    clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
    rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
    alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
    %if(ishear == 1 | ipltcs==1);
    %  cslim = pmin(4)+cumsum((pmax(4)-pmin(4))/NA*ones(1,NC2));
    %  aslim = pmin(5)+cumsum((pmax(5)-pmin(5))/NA*ones(1,NA));
    %end;
    %t1 = toc
    dz  = hmax/(NZ-1);
    z   = cumsum(dz*ones(1,NZ))-dz;
    dzi = hmax/(NZI-1);
    zi  = cumsum(dzi*ones(1,NZI))-dzi;

    h = zeros(1,NZI);  %% h is interface probability
    c = zeros(NPROF,NZ);
    r = zeros(NPROF,NZ);
    a = zeros(NPROF,NZ);
    %if(ishear == 1 | ipltcs==1);
    %  cs = zeros(NPROF,NZ);
    %  as = zeros(NPROF,NZ);
    %end;
    for iprof = 1:NPROF
        %iprof
       ki = k(iprof);
       %if(rem(iprof,10000)==0)
       %   fprintf(1,'%8i',iprof)
       %end
       idxh = zeros(NLMX+1,1);
       prof = zeros(NLMX+1,4);
       %clear prof;
       %% Find index for current model
       %if(ki > 0)
          idxh(1:ki) = (([1:ki]-1)*NPL)+1;
          idxh(ki+1) = idxh(ki);
          %stop
       %else
       %   idxh = [];
       %end
       %% Compute the profile for current model
       %if(ki > 0)
          %prof(1:ki,1) = cumsum(m(iprof,idxh(1:end-1)),2);
          prof(1:ki,1) = m(iprof,idxh(1:ki));
          prof(ki+1,1) = prof(ki,1)+m(iprof,idxh(ki+1));
       %else
       %   prof(1,1) = hmax;
       %end
       prof(1:ki+1,2) = cp2(1:ki+1,iprof);
       if(ilogalf == 0);
         prof(1:ki+1,4) = alfp2(1:ki+1,iprof);
       else;
         prof(1:ki+1,4) = log10(alfp2(1:ki+1,iprof));
       end;
       %if(ishear == 1 | ipltcs==1);
       %  prof(:,5) = cs(1:ki+1,ipltidx,iprof);
       %  if(ilogalf == 0);
       %    prof(:,6) = (alfs(1:ki+1,ipltidx,iprof));
       %  else;
       %    prof(:,6) = log10(alfs(1:ki+1,ipltidx,iprof));
       %  end;
       %end;
       prof(1:ki+1,3) = rho(1:ki+1,iprof);

       for ilay=2:ki+1  %% k is # layers of current model
          idxzi = round(prof(ilay-1,1)/dzi);
          if(idxzi == 0);idxzi = 1;end;
          if(idxzi > NZI);idxzi = NZI;end;
          h(idxzi)     = h(idxzi) + 1;
       end;
       c(iprof,:) = prof(1,2);
       r(iprof,:) = prof(1,3);
       a(iprof,:) = prof(1,4);
       %if(ishear == 1 | ipltcs==1);
       %  cs(iprof,:) = prof(1,5);
       %  as(iprof,:) = prof(1,6);
       %end;
       for ilay=2:ki+1  %% k is # layers of current model
         idxz = round(prof(ilay-1,1)/dz);
         if(idxz == 0);idxz = 1;end;
         c(iprof,idxz:end) = prof(ilay,2);
         r(iprof,idxz:end) = prof(ilay,3);
         a(iprof,idxz:end) = prof(ilay,4);
         %if(ishear == 1 | ipltcs==1);
         %  cs(iprof,idxz:end) = prof(ilay,5);
         %  as(iprof,idxz:end) = prof(ilay,6);
         %end;
       end;
     end;
    %fprintf(1,'\n')

    %
    % Compute histograms for each depth
    %
    %t2 = toc
    %disp('Computing histograms...');
    for iz=1:NZ
       [Nc(iz,:),binsc] = hist(c(:,iz),clim);
       [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
       [Na(iz,:),binsa] = hist(a(:,iz),alim);
       %if(ishear == 1 | ipltcs==1);
       %  [Ncs(iz,:),binscs] = hist(cs(:,iz),cslim);
       %  [Nas(iz,:),binsas] = hist(as(:,iz),aslim);
       %end;
       %if(ipltidx == 36);
         [nfc(ipltidx2,iz,:)] = hpd(c(:,iz),100,95);
         [nfr(ipltidx2,iz,:)] = hpd(r(:,iz),100,95);
         [nfa(ipltidx2,iz,:)] = hpd(a(:,iz),100,95);
         meac(ipltidx2,iz) = median(c(:,iz));
         mear(ipltidx2,iz) = median(r(:,iz));
         meaa(ipltidx2,iz) = median(a(:,iz));
       %end;
     end;
    %
    % Normalize Histograms
    %
    %t3 = toc
    %disp('Normalizing histograms...');
    if(inorm == 0)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/NPROF;
         Nr(iz,:) = Nr(iz,:)/NPROF;
         Na(iz,:) = Na(iz,:)/NPROF;
         %if(ishear == 1 | ipltcs==1);
         %  Ncs(iz,:) = Ncs(iz,:)/NPROF;
         %  Nas(iz,:) = Nas(iz,:)/NPROF;
         %end;
       end;
    elseif(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
         Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
         %if(ishear == 1 | ipltcs==1);
         %  Ncs(iz,:) = Ncs(iz,:)/trapz(binscs,Ncs(iz,:));
         %  Nas(iz,:) = Nas(iz,:)/trapz(binsas,Nas(iz,:));
         %end;
       end;
    elseif(inorm == 2)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         Na(iz,:) = Na(iz,:)/max(Na(iz,:));
         %if(ishear == 1 | ipltcs==1);
         %  Ncs(iz,:) = Ncs(iz,:)/max(Ncs(iz,:));
         %  Nas(iz,:) = Nas(iz,:)/max(Nas(iz,:));
         %end;
       end;
    end;
    %if(ipltidx == 36);
    %  save hpd.mat fr ipltidx z nfc nfa nfr meac mear meaa;
    %end;

    %save tmp.mat

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %% PLOT PROFILE MARGINALS
    %%
    opts = struct('bounds','tight','LockAxes',1, ...
                  'Width',8,'Height',4.8,'Color','cmyk',...
                  'Renderer','painters','Format','png',...
                  'FontMode','fixed','FontSize',14,'FontEncoding','adobe');
    loc = [0.08,  0.18, 0.4533, 0.7266;...
          0.14, 0.14, 0.14, 0.14];
    spw1 = 0.09;
    spw2 = 0.2633;
    sph = 0.82;
    if(ishear==1 | ipltcs==1);
      loc = [0.075, 0.163, 0.331, 0.499, 0.667, 0.835;...
             0.14,   0.14,  0.14,  0.14,  0.14, 0.14];
      spw1 = 0.08;
      spw2 = 0.16;
      sph  = 0.82;
    end;
    %t4 = toc
%    disp('plotting...');
    fig1 = figure();hold on; box on;
    figw = 18;
    figh = 10;
    %set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
    fig1.PaperUnits = 'inches';
    fig1.PaperPosition = [0 0 figw figh];
    %set(fig1, 'renderer', 'painters')
    h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
    hold on; box off;
    h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
    hold on; box off;
    h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
    hold on; box off;
    h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
    hold on; box off;
    if(ishear==1 | ipltcs==1);
      h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]);
      hold on; box off;
      h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph]);
      hold on; box off;
    end;

    subplot(h4);hold on;box on;
    h = h/sum(h);
    h = [0,h,0];
    zi = [zi(1),zi,zi(end)];
    [xx,yy]=stairs(h,zi,'-k');
    patch(xx,yy,[0.8,0.8,0.8]);
    [xx,yy]=stairs(h,zi,'-k');
    set(gca,'Fontsize',14,'XLim',[0 0.08],'YLim',[0 hmax2]);
    set(gca,'YDir','reverse','TickDir','out');
    xlabel('Interface probability');
    ylabel('Depth (m)');
    box on;

    subplot(h1)
    set(gca,'FontSize',14);
    pcolor(clim,z,Nc);shading flat;
    set(h1,'layer','top')
    set(gca,'Fontsize',14,'XLim',[pminplt(1) pmaxplt(1)],'YLim',[0 hmax2]);
    set(gca,'YDir','reverse','TickDir','out');
    xlabel('V_P (m/s)');
    if(isyn == 1)
        plprof(mtru2,hmax,NPLR,'-w',1,0);
        plprof(mtru2,hmax,NPLR,'--k',1,0);
    end;
    set(gca,'YTickLabel',[]);
    set(gca,'XTick',[0:100:5000]);
    box on;
    subplot(h2)
    set(gca,'FontSize',14);
    pcolor(rlim,z,Nr);shading flat;
    %surf(rlim,z,Nr);shading flat;
    set(h2,'layer','top')
    set(gca,'Fontsize',14,'XLim',[pminplt(2) pmaxplt(2)],'YLim',[0 hmax2]);
    set(gca,'YDir','reverse','TickDir','out');
    xlabel('Density (g/ccm)');
    if(isyn == 1)
       plprof(mtru2,hmax,NPLR,'-w',2,0);
       plprof(mtru2,hmax,NPLR,'--k',2,0);
    end
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[1.0:0.2:3.0]);
    set(gca,'XTick',[1.0:0.2:3.0]);
    box on;
    subplot(h3)
    set(gca,'FontSize',14);
    pcolor(alim,z,Na);shading flat;
    %surf(alim,z,Na);shading flat;
    set(h3,'layer','top')
    set(gca,'Fontsize',14,'XLim',[pminplt(3) pmaxplt(3)],'YLim',[0 hmax2]);
    set(gca,'YDir','reverse','YTickLabel',[],'TickDir','out');
    set(gca,'XTick',[-10:2:0]);
    xlabel('log10(\alpha_P (dB/m/kHz))');
    if(isyn == 1)
        plprof(mtru2,hmax,NPLR,'-w',3,0);
        plprof(mtru2,hmax,NPLR,'--k',3,0);
    end;
    cmap = colormap(jet);
    colormap(cmap);
    box on;
    if(ishear == 1 | ipltcs==1);
      subplot(h5)
      set(gca,'FontSize',14);
      pcolor(cslim,z,Ncs);shading flat;
      set(h5,'layer','top')
      set(gca,'Fontsize',14,'XLim',[0 100],'YLim',[0 hmax2]);
    %  set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax2]);
      set(gca,'YDir','reverse','TickDir','out');
      set(gca,'YTickLabel',[]);
      cmap = colormap(jet);
      colormap(cmap);
      box on;
      if(isyn == 1)
        plprof(mtru2,hmax,NPLR,'-w',4,0);
        plprof(mtru2,hmax,NPLR,'--k',4,0);
      end;
      xlabel('V_S (m/s)');
      %set(gca,'XTick',[0:100:2500]);
      set(gca,'XTick',[0:50:2500]);

      subplot(h6)
      set(gca,'FontSize',14);
      pcolor(aslim,z,Nas);shading flat;
      set(h6,'layer','top')
      set(gca,'Fontsize',14,'XLim',[pminplt(5) pmaxplt(5)],'YLim',[0 hmax2]);
      set(gca,'YDir','reverse','TickDir','out');
      set(gca,'YTickLabel',[]);
      cmap = colormap(jet);
      colormap(cmap);
      box on;
      xlabel('log10(\alpha_S (dB/m/kHz))');
      %set(gca,'XTick',[0:5:30]);
      set(gca,'XTick',[-10:1:3]);
      if(isyn == 1)
        plprof(mtru2,hmax,NPLR,'-w',5,0);
        plprof(mtru2,hmax,NPLR,'--k',5,0);
      end;
    end;
    if(icore == 1)
      load core.mat
      subplot(h1);
      %stairs([c1(:,2);c1(end,2)],[c1(:,1);hmax2],'k','Linewidth',1.5);
      %if(site == 2)
       plot(c1(:,2),c1(:,1),'-k','Linewidth',1);
       plot(c1(:,2),c1(:,1),'--w','Linewidth',1);
       plot(c2(:,2),c2(:,1),'-k','Linewidth',1);
       plot(c2(:,2),c2(:,1),'--w','Linewidth',1);
       plot(c3(:,2),c3(:,1),'-k','Linewidth',1);
       plot(c3(:,2),c3(:,1),'--w','Linewidth',1);
       %plot(c4(:,2),c4(:,1),'-k','Linewidth',1);
       %plot(c4(:,2),c4(:,1),'-.w','Linewidth',1);
       %end
       subplot(h2);
       %     stairs([r1(:,2);r1(end,2)],[r1(:,1);hmax2],'k','Linewidth',1.5);
       %     if(site == 2)
       plot(r1(:,2),r1(:,1),'-k','Linewidth',1);
       plot(r1(:,2),r1(:,1),'--w','Linewidth',1);
       plot(r2(:,2),r2(:,1),'-k','Linewidth',1);
       plot(r2(:,2),r2(:,1),'--w','Linewidth',1);
       plot(r3(:,2),r3(:,1),'-k','Linewidth',1);
       plot(r3(:,2),r3(:,1),'--w','Linewidth',1);
       %end
       %subplot(h3);
       %stairs([a1(:,2);a1(end,2)],[a1(:,1);hmax2],'k','Linewidth',1.5);
    end;
    %t5 = toc
    if(isave == 1)
        %saveas(fig1,strcat(plotfile1,plotext1),'fig');
        %saveas(fig1,strcat(plotfile1,plotext2),'png');
        print(fig1,strcat(plotfile1,plotext2),'-dpng','-r0');
        if(IEPS == 1)
          print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
        end;
    end;
    %t6 = toc
    close(fig1)
end; % end ipltidx

end; %imarg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% COMPUTE PROFILE MARGINALS for dispersion values
%%
if(idispmarg == 1)
  %disp('Plot dispersion profiles...');
  clear h c r a;
  if(isyn == 1);
    for ilay = 1:ktru+1;
      ipar1 = (ilay-1)*NPL+2;
      ipar2 = (ilay-1)*NPLR+2;
      if(ilay == ktru+1);
        ipar1 = (ilay-1)*NPL+1;
        ipar2 = (ilay-1)*NPLR+1;
      end;
      if(ilay < ktru+1);mtru3(ipar2-1)  = mtru(ipar1-1);end;
      mtru3(ipar2)   = cptru(ilay,end)-cptru(ilay,1);
      if(ilogalf==0)
        mtru3(ipar2+2) = (abs(alfptru(ilay,end)-alfptru(ilay,1)));
      else;
        mtru3(ipar2+2) = log10(abs(alfptru(ilay,end)-alfptru(ilay,1)));
      end;
      mtru3(ipar2+1) = rhotru(ilay)/1000.;
    end;
  end;
  nsmooth = ceil(NZ/80.);
  clim = pmin2(1)+cumsum((pmax2(1)-pmin2(1))/NC*ones(1,NC));
  rlim = pmin2(2)+cumsum((pmax2(2)-pmin2(2))/NR*ones(1,NR));
  alim = pmin2(3)+cumsum((pmax2(3)-pmin2(3))/NA*ones(1,NA));

  dz = hmax/(NZ-1);
  z = cumsum(dz*ones(1,NZ))-dz;

  h = zeros(1,NZ);
  c = zeros(NPROF,NZ);
  r = zeros(NPROF,NZ);
  a = zeros(NPROF,NZ);
  cp1 = squeeze(cp(:,1,:));
  cpe = squeeze(cp(:,end,:));
  alfp1 = squeeze(alfp(:,1,:));
  alfpe = squeeze(alfp(:,end,:));

  for iprof = 1:NPROF
    ki = k(iprof);
    %if(rem(iprof,10000)==0)
    %   fprintf(1,'%8i',iprof)
    %end
    idxh = zeros(NLMX+1,1);
    prof = zeros(NLMX+1,4);
    %% Find index for current model
    %if(ki > 0)
       idxh(1:ki) = (([1:ki]-1)*NPL)+1;
       idxh(ki+1) = idxh(ki);
    %else
    %   idxh = [];
    %end
    %% Compute the profile for current model
    if(ki > 0)
       %prof(1:ki,1) = cumsum(m(iprof,idxh(1:end-1)),2);
       prof(1:ki,1) = m(iprof,idxh(1:ki));
       prof(ki+1,1) = prof(ki,1)+m(iprof,idxh(ki+1));
    else
       prof(1,1) = hmax;
    end
    prof(1:ki+1,2) = cpe(1:ki+1,iprof)-cp1(1:ki+1,iprof);
    prof(1:ki+1,3) = rho(1:ki+1,iprof);
    prof(1:ki+1,4) = alfpe(1:ki+1,iprof)-alfp1(1:ki+1,iprof);
    c(iprof,:) = prof(1,2);
    r(iprof,:) = prof(1,3);
    if(ilogalf==0)
      a(iprof,:) = abs(prof(1,4));
    else;
      a(iprof,:) = log10(abs(prof(1,4)));
    end;
    for ilay=2:ki+1  %% k is # layers of current model
       idxz = round(prof(ilay-1,1)/dz);
       if(idxz == 0);idxz=1;end;
       %h(idxz)     = h(idxz) + 1;
       c(iprof,idxz:end) = prof(ilay,2);
       r(iprof,idxz:end) = prof(ilay,3);
       if(ilogalf==0)
         a(iprof,idxz:end) = abs(prof(ilay,4));
       else;
         a(iprof,idxz:end) = log10(abs(prof(ilay,4)));
       end;
    end;
  end;
  %fprintf(1,'\n')

   %
   % Compute histograms for each depth
   %
   for iz=1:NZ
%      Nc(iz,:) = histc_tstar(c(:,iz),logL,clim,Tstar);
%      Nr(iz,:) = histc_tstar(r(:,iz),logL,rlim,Tstar);
%      Na(iz,:) = histc_tstar(a(:,iz),logL,alim,Tstar);
      Nc(iz,:) = histc(c(:,iz),clim);
      Nr(iz,:) = histc(r(:,iz),rlim);
      Na(iz,:) = histc(a(:,iz),alim);
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
         Na(iz,:) = Na(iz,:)/max(Na(iz,:));
      end;
   elseif(inorm == 2)
      Nc = Nc/max(max(Nc));
      Nr = Nr/max(max(Nr));
      Na = Na/max(max(Na));
   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT PROFILE MARGINALS
%%
%nx = 4;
%ny = 1;
%xim = 0.01;
%yim = 0.06;
%xymarg = [0.1 0.04 0.04 0.14];
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',14,'FontEncoding','adobe');
%[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

loc = [0.08,  0.18, 0.4533, 0.7266;...
      0.14, 0.14, 0.14, 0.14];
spw1 = 0.09;
spw2 = 0.2633;
sph = 0.82;

fig2 = figure();hold on; box on;
figw = 18;
figh = 10;
%set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 figw figh];
set(fig2, 'renderer', 'painters')
h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
hold on; box off;
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;

subplot(h4);hold on;box on;
if(idispmarg == 1)
   h = h/sum(h);
   [xx,yy]=stairs(h,z,'-k');
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,z,'-k');
   set(gca,'Fontsize',14,'XLim',[0 0.04],'YLim',[0 hmax2]);
   set(gca,'YDir','reverse');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;

subplot(h1)
set(gca,'FontSize',14);
if(idispmarg == 1)
  pcolor(clim,z,Nc);shading flat;
  for izz=1:length(zplt);
    plot([pmin2(1) pmax2(1)],[zplt(izz) zplt(izz)],'-k');
    plot([pmin2(1) pmax2(1)],[zplt(izz) zplt(izz)],'--w');
  end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin2(1) pmax2(1)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Vp Dispersion (m/s)');

if(isyn == 1)
if(idispmarg == 1)
   plprof(mtru3,hmax,NPLR,'-w',1,0);
   plprof(mtru3,hmax,NPLR,'--k',1,0);
end;
end;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[0:25:100]);
set(gca,'XTick',[0:25:100]);
box on;
subplot(h2)
set(gca,'FontSize',14);
if(idispmarg == 1)
   pcolor(rlim,z,Nr);shading flat;
end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin2(2) pmax2(2)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('Density (g/ccm)');
if(isyn == 1)
if(idispmarg == 1)
   plprof(mtru3,hmax,NPLR,'-w',2,0);
   plprof(mtru3,hmax,NPLR,'--k',2,0);
end
end
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1.4 1.6 1.8]);
set(gca,'XTick',[1.4 1.6 1.8]);
box on;
subplot(h3)
set(gca,'FontSize',14);
if(idispmarg == 1)
  pcolor(alim,z,Na);shading flat;
  for izz=1:length(zplt);
    plot([pmin2(3) pmax2(3)],[zplt(izz) zplt(izz)],'-k');
    plot([pmin2(3) pmax2(3)],[zplt(izz) zplt(izz)],'--w');
  end;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin2(3) pmax2(3)],'YLim',[0 hmax2]);
set(gca,'YDir','reverse');
xlabel('log10(Atten. change (dB/m/kHz))');
if(isyn == 1)
if(idispmarg == 1)
   plprof(mtru3,hmax,NPLR,'-w',3,0);
   plprof(mtru3,hmax,NPLR,'--k',3,0);
end
end;
set(gca,'YTickLabel',[]);
set(gca,'XTick',[log10(0.0025),log10(0.018),log10(0.135)]);
set(gca,'XTickLabel',[0.0025,0.018,0.135]);
%cmap = colormap(flipud(gray));
cmap = colormap(jet);
%cmap(1,:) = [1 1 1];
colormap(cmap);
%colorbar('peer',h1,'location','WestOutside');
box on;
if(icore == 1)
  load core.mat
  subplot(h2);
  plot(r1(:,2),r1(:,1),'-k','Linewidth',1);
  plot(r1(:,2),r1(:,1),'--w','Linewidth',1);
%  plot(r2(:,2),r2(:,1),'-k','Linewidth',1);
%  plot(r2(:,2),r2(:,1),'--w','Linewidth',1);
%  plot(r3(:,2),r3(:,1),'-k','Linewidth',1);
%  plot(r3(:,2),r3(:,1),'--w','Linewidth',1);
end;

if(isave == 1)
   saveas(fig2,strcat(plotfile2,plotext2),'png');
   if(IEPS == 1)
     print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
   end;
end;

end; %idispmarg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT DISPERSION CURVES
%%
if(ipltdisp == 1)
disp('Plot dispersion curves...');
clear c a h;
h = zeros(1,NZ);
c = zeros(NPROF2,NZ,NFREQ);
a = zeros(NPROF2,NZ,NFREQ);
cp2   = zeros(NLMX+1,NPROF);
alfp2 = zeros(NLMX+1,NPROF);

for ipltidx = 1:NFREQ+1;
  if(ipltidx == NFREQ+1);ipltidx = ipltidx-1;end;
  cp2 = squeeze(cp(:,ipltidx,:));
  alfp2 = squeeze(alfp(:,ipltidx,:));
  dc = (pmax(1)-pmin(1))/NC;
  da = (pmax(3)-pmin(3))/NA;
  clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
  alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
  dz = hmax/(NZ-1);
  z = cumsum(dz*ones(1,NZ))-dz;

  if(isyn == 1);
    for ilay = 1:ktru+1;
      ipar1 = (ilay-1)*NPL+2;
      ipar2 = (ilay-1)*NPLR+2;
      if(ilay == ktru+1);
        ipar1 = (ilay-1)*NPL+1;
        ipar2 = (ilay-1)*NPLR+1;
      end;
      if(ilay < ktru+1);
        mtru4(ipar2-1)  = mtru(ipar1-1);
        ztru(ilay) = mtru(ipar1-1);
      else
        ztru(ilay) = hmax;
      end;
      mtru4(ipar2)    = cptru(ilay,ipltidx);
      mtru4(ipar2+1)  = rhotru(ilay)/1000.;
      if(ilogalf==0)
        mtru4(ipar2+2)  = (alfptru(ilay,ipltidx));
      else;
        mtru4(ipar2+2)  = log10(alfptru(ilay,ipltidx));
      end;
    end;
    for ilay = 1:length(zplt)
      idx = find(ztru-zplt(ilay)>0);
      if(isempty(idx));ilaytru(ilay)=length(ztru);else;ilaytru(ilay)=idx(1);end;
      clear idx;
    end;
  end;

   for iprof2 = 1:NPROF2;
%     iprof = idxran(iprof2);
     iprof = iprof2;
     ki = k(iprof);
     %if(rem(iprof2,10000)==0)
     %   fprintf(1,'%8i',iprof2)
     %end
     idxh = zeros(NLMX,1);
     prof = zeros(NLMX,4);
     %% Find index for current model
     %if(ki > 0)
        idxh(1:ki) = (([1:ki]-1)*NPL)+1;
        idxh(ki+1) = idxh(ki);
     %else
     %   idxh = [];
     %end
     %% Compute the profile for current model
     if(ki > 0)
        prof(1:ki,1) = m(iprof,idxh(1:ki));
        prof(ki+1,1) = prof(ki,1)+m(iprof,idxh(ki+1));
     else
        prof(1,1) = hmax;
     end
     prof(1:ki+1,2) = cp2(1:ki+1,iprof);
     prof(1:ki+1,4) = alfp2(1:ki+1,iprof);

     c(iprof2,:,ipltidx) = prof(1,2);
     if(ilogalf==0)
       a(iprof2,:,ipltidx) = (prof(1,4));
     else;
       a(iprof2,:,ipltidx) = log10(prof(1,4));
     end;
     for ilay=2:ki+1  %% k is # layers of current model
        idxz = round(prof(ilay-1,1)/dz);
        if(idxz == 0);idxz = 1;end;
        h(idxz)     = h(idxz) + 1;
        c(iprof2,idxz:end,ipltidx) = prof(ilay,2);
        if(ilogalf==0)
          a(iprof2,idxz:end,ipltidx) = (prof(ilay,4));
        else;
          a(iprof2,idxz:end,ipltidx) = log10(prof(ilay,4));
        end;
     end;
   end;
   %fprintf(1,'\n')
end;

clear Nc Na clim alim;
fig3 = figure();hold on; box on;
figw = 16;
figh = 10;
%set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 figw figh];
set(fig3, 'renderer', 'painters')
jj = 1;
izidx = round(zplt/dz);
if(isyn == 0);
  clear cptru alpftru cstru alfstru rhotru;
end;
for izz = 1:length(zplt);
  iz = izidx(izz);
  %
  % Compute histograms for each depth
  %
  NV = 400;
  cmin = ylimcmin(izz);
  cmax = ylimcmax(izz);
  clim = cmin+cumsum((cmax-cmin)/NV*ones(1,NV));
  amin = ylimamin(izz);
  amax = ylimamax(izz);
  alim = amin+cumsum((amax-amin)/NV*ones(1,NV));
  for ifr=1:NFREQ
     Nc(:,ifr) = histc(c(:,iz,ifr),clim);
     Na(:,ifr) = histc(a(:,iz,ifr),alim);
  end;
  %
  % Normalize Histograms
  %
  if(inorm == 1)
     for ifr=1:NFREQ
        Nc(:,ifr) = Nc(:,ifr)/max(Nc(:,ifr));
        Na(:,ifr) = Na(:,ifr)/max(Na(:,ifr));
     end;
  elseif(inorm == 2)
     Nc = Nc/max(max(Nc));
     Na = Na/max(max(Na));
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% PLOT PROFILE MARGINALS
  %%
  nx = 2;
  ny = length(zplt);
  xim = 0.05;
  yim = 0.05/ny;
  xymarg = [0.07 0.01 0.01 0.1];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
  difffr = diff(fr)/2;
  frplt = [fr(1)-2*difffr(1),fr-[difffr,difffr(end)],...
           fr(end)+2*difffr(end)];

  subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
  set(gca,'FontSize',14);
  if(izz==3);ylabel('Velocity (m/s)');end;
  set(gca,'XTick',bands2);
  set(gca,'XTickLabel',[]);
  if(izz == length(zplt));
    set(gca,'XTickLabel',bands2);
    xlabel('Frequency (Hz)');
  end;
  set(gca,'XLim',[fr(1) fr(end)]);
%  set(gca,'YLim',[min(min(c(:,iz,:))) max(max(c(:,iz,:)))]);
  set(gca,'YLim',[ylimcmin(izz) ylimcmax(izz)]);
  set(gca,'FontSize',14);
  set(gca,'LineWidth',1);
  imagesc(frplt,clim,Nc);shading flat;
  if(isyn == 1);
    semilogx(fr,cptru(ilaytru(izz),:),'-k');
    semilogx(fr,cptru(ilaytru(izz),:),'--w');
  end;
  set(gca,'XScale','log');
%  text(400,max(max(c(:,iz,:)))-(max(max(c(:,iz,:)))-min(min(c(:,iz,:))))/10,[num2str(round(100.*z(iz))/100.) ' m'],'FontSize',12,'Color',[1,1,1])
  text(2000,ylimcmax(izz)-(ylimcmax(izz)-ylimcmin(izz))/10,[num2str(round(100.*z(iz))/100.) ' m'],'FontSize',16,'Color',[1,1,1])

  subplot('Position',[loc(1,jj+1) loc(2,jj+1) spw sph]);hold on;box on;
  set(gca,'FontSize',14);
  if(izz == 3);ylabel('log10(Attenuation (dB/m/kHz))');end;
  set(gca,'XTick',bands2);
  set(gca,'XTickLabel',[]);
  if(izz == length(zplt));
    set(gca,'XTickLabel',bands2);
    xlabel('Frequency (Hz)');
  end;
  set(gca,'XLim',[fr(1) fr(end)]);
%  set(gca,'YLim',[min(min(a(:,iz,:))) max(max(a(:,iz,:)))]);
  set(gca,'YLim',[ylimamin(izz) ylimamax(izz)]);
  set(gca,'YTick',[-5,-3,-1,1,3,5]);
  set(gca,'FontSize',14);
  set(gca,'LineWidth',1);
  imagesc(frplt,alim,Na);shading flat;
  if(isyn == 1);
    semilogx(fr,log10(alfptru(ilaytru(izz),:)),'-k');
    semilogx(fr,log10(alfptru(ilaytru(izz),:)),'--w');
  end;
  set(gca,'XScale','log');
  clear Nc Na clim alim;

  jj = jj + 2;
end; % end izz

if(isave == 1)
%  saveas(fig3,strcat(plotfile3,plotext1),'fig');
  saveas(fig3,strcat(plotfile3,plotext2),'png');
  if(IEPS == 1)
    print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
  end;
end;
end; %ipltdisp
disp('done');

%return;
