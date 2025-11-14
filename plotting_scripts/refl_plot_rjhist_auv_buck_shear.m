%function [p_mean,p_mead] = plot_rjhist_auv_buck_shear(filename);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

%filename = 'blge5_3_bl_ping_id1_IvgsEven_5Fres1d5_sample.mat';
filename = 'blge5_3_bl_ping_id1_7oct_verbose_sample.mat';
p_mean=0;
p_mead=0;
imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 0; %% Normalize profile marginals line by line
isyn    = 0;
ivoro   = 0; %% 1 -> Voronoi Cells; 2 -> Layer nodes
imap    = 0;
imead   = 0;
imean   = 0;
imax    = 0;
isave   = 1;
idatfit = 0;
irex    = 0;
idathist= 0;
iboudreau= 0;
i_vref  = 1;   %% Assumes reference model is employed in inversion
ishear  = 0;
ivpvs   = 0;
iar     = 1;
isd     = 1;
imagscale = 1;
ibucking= 3;
%ibucking= 0;
ivgsnosh= 0;
IEPS = 0;
NVMX = 20;
NMISC = 1;
if(ibucking == 0);
  NPL = 4;
  NPV = 4;
  if(ishear == 1);NPL = 6;NPV=6;end;
elseif(ibucking == 3 & ivgsnosh == 0);
  NPL = 5;
  NPV = 5;
  if(ishear == 1);NPL = 6;NPV=6;end;
  if(ivgsnosh == 1)
    NMISC = 5;
  else;
     NMISC = 4;
  end;
else;
  NPL = 4;
end;

%smpext      = '_particles.mat';
smpext      = '_sample.mat';
filebase    = strrep(filename,smpext,'');
datafile    = strcat(filebase,'.txt');
vreffile    = strcat(filebase,'_vel_ref.txt');
%vreffile    = 'p0004_pave010_0988_2513_vel_ref.txt';
vorofile    = strcat(filebase,'_voro_sample.mat');
repfile     = strcat(filebase,'_rep.dat');
mapfile     = strcat(filebase,'_map.dat');
if(ivoro >= 1);mapfile = strcat(filebase,'_map_voro.dat');end;
profilefile = strcat(filebase,'_profile.mat');
resfile     = strcat(filebase,'_res_ens.dat');
resarfile   = strcat(filebase,'_resar_ens.dat');
repensfile  = strcat(filebase,'_rep_ens.dat');

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';
plotfile0   = strcat(filebase,'complexity.');
plotfile1   = strcat(filebase,'chains_logl.');
plotfile2   = strcat(filebase,'ar_par.');
plotfile3   = strcat(filebase,'sigma.');
plotfile33  = strcat(filebase,'misc.');
plotfile4   = strcat(filebase,'number_interfaces.');
plotfile12  = strcat(filebase,'number_interfaces_varpar.');
plotfile5   = strcat(filebase,'logL.');
plotfile6   = strcat(filebase,'transdim_marg_prof.');
plotfile66  = strcat(filebase,'transdim_map_prof.');
plotfile10  = strcat(filebase,'transdim_node_dens.');
plotfile7   = strcat(filebase,'datafit.');
plotfile77  = strcat(filebase,'datafit2.');
plotfile8   = strcat(filebase,'axx.');
plotfile9   = strcat(filebase,'reshist.');
corefile    = 'core.mat';

NAVEF = 7;
IGA = 0;
oct = 15.;
frbw = 1.0/15.0;
%FBW  = 37.5;
FBW  = 50.;
icore   = 0;
ibl     = 0;
nffile = 'nf_band3.mat';
hmin = 0.001;

%% AUV
%bands = [988., 1113., 1288., 1913., 2263., 2513.];

%% MUD site
%bands = [400. 460. 528. 607. 697. 801. 920. 1057. 1214.];


%% Site 1  8 m packet
%bands = [400 504  635  800 1008 1270];
%bands = [400 504  635  800 1008 1270 1600 2016];

%% Site 02 8.4 m Elba
%bands = [500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3203.];
%bands = [300.,400.,500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.];
%bands = [318. 349. 400. 439. 504. 553.  636. 697. 801. 878. 1009. 1106. 1271. 1394. 1601. 1756. 2018. 2213. 2542. 2788. 3203. 3513. 4035.];

%% Site 04 13 m Malta
%bands = [400.,600.,800.,1000.,1200.,1400.,1600.];

%% Site 04 13 m Malta
%bands = [550.,700.,900.,1100.,1400.,1800.];

%% Site 13  8 m packet
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000 6300];

%% Tyrrhenian Site 03
%bands = [400 504  635  800 1008 1270 1600];
%bands = [400 504  635  800 1008 1270 1600 2016 2540 3200];
%% 6m packet
%bands = [504  635  800 1008 1270 1600 2016];

%% Tyrrhenian Site 05
bands = [464.37  566.08  690.06  928.75  1132.15  1250];
%bands = [464.37 690.06 928.75 1250.0 1682.38];

%% Site 16  8 m packet
%bands = [504  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  600  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032 5000 6300];
%bands = [504  635  800 1008 1270 1600 2016 2540 3200 4032];
%bands = [504  635  800 1008 1270 1600 2016 2540];

%% Site 07
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250.,1600.,2000.,2500.,3150.];
%bands = [ 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000.];
%bands = [ 100., 125., 160., 200., 250., 315., 400., 500., 630., ...
%          800.,1000.,1250., 1600.,2000.,2500.,3150.,4000.,5000.,6300.,8000.,10000.];
%% Site 16  8 m packet
%bands = [500.,630.,800.,1000.,1250.,1600.,2000.,2500.,3150.,4000.,5000.];
NBAND = length(bands);
NSD = NBAND;
%% AUV:
%pmin = [1420 1.2 0.0]';
%pmax = [1800 2.4 1.0]';
%% Site07:
%pmin = [1550 1.4 0.0]';
%pmax = [1800 2.3 1.0]';
%% Bucking simulation:
if(ibucking == 0)
  %pmin = [1450 1.15 -3.]';
  %pmax = [1610 2.00 -.3]';
  pmin = [1450 1.1 -3.0]';
  pmax = [2000 2.5  1.0]';
  if(isyn == 1)
    pmin = [1450 1.20 0  ]';
    pmax = [1750 2.20 1.0]';
    %pmin = [1400 1.2 0    log(sqrt(2))  0]';
    %pmax = [1700 2.0 0.1  log(100)  5]';
  end
elseif(ibucking == 1)
  pmin = [0.25 1.0 log(50.00)]';
  pmax = [0.90 2.1 log(2.e6)]';
elseif(ibucking == 2)
  pmin = [0.25 log(2.e7) 0.05]';
  pmax = [0.90 log(8.e8) 0.15]';
elseif(ibucking == 3)
  if(ivgsnosh == 0)
    %% Site 03 Tyrrhenian Sea:
    %pmin = [0.10 log(1.e6)  log(0.001) log(4.53e-5) log(1.e6)]';
    %pmax = [0.91 log(1.e11) log(0.05) log(0.4)   log(1.e11)]';
    %% Site 02 Malta:
    %pmin = [0.30 log(1.e6)  log(0.001) log(3.35e-4) log(1.e6)]';
    %pmax = [0.91 log(5.e8) log(0.147) log(0.4)   log(5.e8)]';
    %% Var-Par sim:
    %pmin = [0.10 log(1.e6)  log(0.001) log(8.94e-5) log(5.e4)]';
    %pmax = [0.91 log(5.e9) log(0.147) log(0.4)   log(5.e9)]';
    %% VGS AUV:
    pmin = [0.00 6.9  -1.4 -4.5]';
    pmax = [0.95 8.75 -0.4 -1.3]';
  else;
    pmin = [0.20 log(2.e7) -10]';
    pmax = [0.91 log(8.e8) -3]';
  end;
elseif(ibucking == 4)
  pmin = [0.20 1.0 log(1.e-17)]';
  pmax = [0.91 2.0 log(1.e-9)]';
end;
if(ibucking > 0)
  miscmin = [3.55e10, 2.3510e9,2630.,1028.]';
  miscmax = [3.65e10, 2.3520e9,2700.,1031.]';
%  miscmin = [3.25e10, 2.0000e9,2500.,1010.]';
%  miscmax = [3.95e10, 2.5000e9,2800.,1040.]';
%  miscmin = [2.095e10, 2.3400e9,2680.,1025.]';
%  miscmax = [2.105e10, 2.3500e9,2720.,1030.]';
end;
%   NDAVE   = ceil(NPROF/10);
%NDAVE   = 500;
thinstep = 1;

%% Charles' VGS simulations:
%mtrumisc = [3.6136468e10,2.3780926e9,2722.1565,1030.8898];
%mtru = [0.4,  0.76,  log(3.0e7),  0.07,  log(1.2E-2),...
%        1.0,  0.43,  log(8.5e7),  0.1,   log(2.2E-4),...
%        2.2,  0.58,  log(3.5e7),  0.075, log(3E-3),...
%        3.2,  0.385, log(3.8e8),  0.12,  log(1.2E-4),...
%              0.47,  log(5.9e7),  0.08,  log(3E-4)];

%%
%% VGS shear simulation
%%
if(ibucking == 3);
  mtru = [0.29, 0.800, 16.000, -6.6036, -4.173,...
          0.91, 0.700, 16.000, -6.6036, -4.173,...
          1.89, 0.600, 18.000, -6.6036, -4.173,...
          3.50, 0.500, 21.000, -6.6036, -5.210,...
                0.400, 21.000, -6.6036, -5.210 ];
  mtrumisc = [  36234377700., 2353069330., 2737.9595, 1029.39201];

  if(ivoro >= 1);
  if(ishear== 1);
    %% True parameters for simulation:
    vorotru(1,:) = [0.000, 0.800, 18.000, -6.6036, -4.173, 14.000];
    vorotru(2,:) = [0.500, 0.700, 19.000, -6.6036, -1.037, 14.000];
    vorotru(3,:) = [1.000, 0.600, 20.000, -6.6036, -6.203, 14.000];
    vorotru(4,:) = [1.500, 0.500, 21.000, -6.6036, -5.210, 18.600];
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
  else;
    %% True parameters for simulation:
    vorotru(1,:) = [0.000, 0.800, 16.000, -6.6036, -4.173];
    vorotru(2,:) = [0.290, 0.700, 16.000, -6.6036, -4.173];
    vorotru(3,:) = [0.910, 0.600, 18.000, -6.6036, -4.173];
    vorotru(4,:) = [1.890, 0.500, 21.000, -6.6036, -5.210];
    vorotru(5,:) = [3.500, 0.400, 21.000, -6.6036, -5.210];
    voroidxtru(1,:) = [1, 1, 1, 1, 1];
    voroidxtru(2,:) = [1, 1, 0, 0, 0];
    voroidxtru(3,:) = [1, 1, 1, 0, 0];
    voroidxtru(4,:) = [1, 1, 1, 0, 1];
    voroidxtru(5,:) = [1, 1, 0, 0, 0];
    if(ivoro == 1);
      [mtru,nuniquetru] = refl_voro_to_lay(vorotru,voroidxtru,5,5,5);
    elseif(ivoro == 2);
      [mtru,nuniquetru] = refl_laynode_to_lay(vorotru,voroidxtru,5,5,5);
    end;
    % 0.1500      0.8000     18.0000     -6.6036     -4.1730     14.0000
    % 0.4500      0.7000     19.0000     -6.6036     -1.0370     14.0000
    % 0.8000      0.6000     20.0000     -6.6036     -6.2030     16.0000
    % 3.0450      0.5000     21.0000     -6.6036     -5.2100     18.6000
    %             0.4000     22.0000     -3.9380     -5.2100     20.0000
    mtrumisc = [36234377700., 2353069330., 2737.9595, 1029.39201];


  end;
  end;
elseif(ibucking == 0);
  if(ivoro > 0);
  if(ishear == 1);
    vorotru(1,:) = [0.000, 1550., 1.5000, 0.0100,  300.0, 0.400];
    vorotru(2,:) = [0.500, 1600., 1.5000, 0.0100,  300.0, 0.400];
    vorotru(3,:) = [1.000, 1900., 1.9000, 0.0200,  800.0, 0.400];
    vorotru(4,:) = [1.500, 2500., 2.6000, 0.0200, 1600.0, 0.400];
    vorotru(5,:) = [3.500, 2800., 2.6000, 0.0200, 1600.0, 0.400];
    voroidxtru(1,:) = [1, 1, 1, 1, 1, 1];
    voroidxtru(2,:) = [1, 1, 0, 0, 0, 0];
    voroidxtru(3,:) = [1, 1, 1, 1, 1, 0];
    voroidxtru(4,:) = [1, 1, 1, 0, 1, 0];
    voroidxtru(5,:) = [1, 1, 0, 0, 0, 0];
    if(ivpvs == 1);
      vorotru(:,5) = vorotru(:,2)./exp(vorotru(:,5));
    end;
    if(ivoro == 1);
      [mtru,nuniquetru] = refl_voro_to_lay(vorotru,voroidxtru,5,6,5);
    elseif(ivoro == 2);
      [mtru,nuniquetru] = refl_laynode_to_lay(vorotru,voroidxtru,5,6,5);
    end;
  else;
    vorotru(1,:) = [0.000, 1500., 1.4000, 0.22];
    vorotru(2,:) = [0.200, 1530., 1.6000, 0.22];
    vorotru(3,:) = [0.500, 1580., 1.6000, 0.22];
    vorotru(4,:) = [1.500, 1550., 1.6000, 0.22];
    vorotru(5,:) = [3.000, 1600., 1.8000, 0.22];
    voroidxtru(1,:) = [1, 1, 1, 1];
    voroidxtru(2,:) = [1, 1, 1, 0];
    voroidxtru(3,:) = [1, 1, 0, 0];
    voroidxtru(4,:) = [1, 1, 0, 0];
    voroidxtru(5,:) = [1, 1, 1, 0];
    if(ivoro == 1);
      [mtru,nuniquetru] = refl_voro_to_lay(vorotru,voroidxtru,5,4,5);
    elseif(ivoro == 2);
      [mtru,nuniquetru] = refl_laynode_to_lay(vorotru,voroidxtru,5,4,5);
    end;
  end;
  else
    %mtru = [0.200, 1500.0, 1.40, 0.2200,...
    %        0.500, 1530.0, 1.60, 0.2200,...
    %        1.500, 1580.0, 1.60, 0.2200,...
    %        3.000, 1550.0, 1.60, 0.2200,...
    %               1600.0, 1.80, 0.2200];
    mtru = [0.94, 1480., 1.30, 0.01, ...
            1.44, 1606., 1.75, 0.20, ...
            3.10, 1553., 1.51, 0.20, ... 
                  1681., 1.95, 0.10]
  end;
end;
disp('Converting voro to lay');
if(ivoro == 1);
  refl_convert_sample_voro_to_lay(vorofile,NPV,NVMX,NMISC,NBAND)
elseif(ivoro == 2);
  refl_convert_sample_laynode_to_lay(vorofile,NPV,NVMX,NMISC,NBAND)
end;
disp('Done converting voro to lay');
if(ivoro > 0);
  load(vorofile);
  B=A;
end;
load(filename);
if(ivoro == 0);B=A;end;

logLmin = min(A(:,1))-10;
logLmax = max(A(:,1))+10;
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NBAND);
   armax = 1.*ones(1,NBAND);
else
   order = 1;
end;
datafile
tmp = dlmread(datafile);
z_t    = tmp(1,1);
cw     = tmp(2,1);
rw     = tmp(3,1);
hmax   = tmp(4,1)+.5;
%hmax   = 5.;
dobs   = tmp(5:length(bands)+4,:)';
angobs = tmp(length(bands)+5,:);
NANG   = length(angobs);
if(irex == 1);
 rex    = tmp(length(bands)+6:end-1,:)';
else;
 rex    = ones(length(bands),NANG)';
end;
if(icore == 1)
  disp('Loading core data.');
  load(corefile);
end;
NPROF = size(A,1);

%%
%% Random permutation of sample is only applied to data fit computation
%%
disp('Done loading data.');

k = A(:,4);
%% Compute prior volume
logP = zeros(size(k));
for i=1:length(k);
  vol = 1.;
  for ik=1:k(i);
    vol = vol * (hmax-(ik-1.)*hmin)*0.489*(pmax(end)-pmin(end));
  end;
  vol = vol * 0.489 * (pmax(end)-pmin(end));  %% This is half-space volume
  logP(i) = -log(vol);
end;

NFP = (k*NPL)+NPL-1;
m = A(:,5:end-5);
%% Read in ref velocity
if(i_vref == 1);
  vel_ref = dlmread(vreffile);
  vel_ref(1,:) = [];
  for ismp=1:length(m);
    for ik=1:k(ismp)+1;
      if(ik <= k(ismp));
        ipar = (ik-1)*NPL+1;
        z = m(ismp,ipar);
      else;
        ipar = (ik-1)*NPL;
        z = m(ismp,ipar-NPL+1);
      end;
      if(ibucking == 0);
          [vref,rref]=refl_getref(z,vel_ref);
          m(ismp,ipar+1) = vref + m(ismp,ipar+1);
          m(ismp,ipar+2) = rref + m(ismp,ipar+2);
      else;
          [pref]=refl_getrefpor(z,vel_ref);
          m(ismp,ipar+1) = pref + m(ismp,ipar+1);
      end;
    end;
  end;
end;

%% Voronoi nodes sample file
if(ivoro >= 1);
  kv = B(:,4);
  gv = zeros(sum(kv),1);
  kvpar = zeros(size(kv,1),NPV);
  mv = B(:,5:end-5);
  jj1 = 1;jj2 = 1;jj3 = 1;jj4 = 1;jj5 = 1;
  igv = 1;
  for i=1:size(B,1);
    if(rem(i,10000)==0)
       fprintf(1,'counted nodes in %8i samples\n',i)
    end
    for j=1:kv(i);
      if(NPV > 1);
      if(mv(i,(j-1)*NPV+2) > -99.);
        por(jj1,:)   = mv(i,[(j-1)*NPV+1,(j-1)*NPV+2]);
        kvpar(i,1) = kvpar(i,1) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj1 = jj1 + 1;
      end;end;
      if(NPV > 2);
      if(mv(i,(j-1)*NPV+3) > -99.);
        Kp(jj2,:)  = mv(i,[(j-1)*NPV+1,(j-1)*NPV+3]);
        kvpar(i,2) = kvpar(i,2) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj2 = jj2 + 1;
      end;end;
      if(NPV > 3);
      if(mv(i,(j-1)*NPV+4) > -99.);
        strain(jj3,:) = mv(i,[(j-1)*NPV+1,(j-1)*NPV+4]);
        kvpar(i,3) = kvpar(i,3) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj3 = jj3 + 1;
      end;end;
      if(NPV > 4);
      if(mv(i,(j-1)*NPV+5) > -99.);
        taup(jj4,:)   = mv(i,[(j-1)*NPV+1,(j-1)*NPV+5]);
        kvpar(i,4) = kvpar(i,4) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj4 = jj4 + 1;
      end;end;
      if(NPV > 5);
      if(mv(i,(j-1)*NPV+6) > -99.);
        Ks(jj5,:) = mv(i,[(j-1)*NPV+1,(j-1)*NPV+6]);
        kvpar(i,5) = kvpar(i,5) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj5 = jj5 + 1;
      end;end;
      %% Position parameter always present:
      kvpar(i,NPV) = kvpar(i,NPV) + 1;
      igv = igv + 1;
    end;
  end;
  gv(find(gv==0))=[];
  nparv = sum(kvpar,2);
  disp('Done counting nodes.');

  %%
  %% Histograms of No. nodes, No. unique layers, No. parameters,
  %% No. parameters per node
  %%
  fig0 = figure();
  nx = 1;
  ny = 4;
  xim = 0.00;
  yim = 0.1;
  xymarg = [0.1 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
  hold on; box on;set(gca,'FontSize',14);
  [n1,lim]=hist(kv,[1:NVMX]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
  n1 = n1/sum(n1);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n1,'k');
  clear n1 lim;
  xlabel('No. Layer nodes');
  set(gca,'XLim',[0 NVMX],'TickDir','out');
  if(isyn == 1);
    plot([size(voroidxtru,1) size(voroidxtru,1)],[0 1],'-w');
    plot([size(voroidxtru,1) size(voroidxtru,1)],[0 1],'--k');
  end;

  h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
  hold on; box on;set(gca,'FontSize',14);
  [n1,lim]=hist(A(:,4),[min(A(:,4)):max(A(:,4))]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
  n1 = n1/sum(n1);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n1,'k');
  clear n1 lim;
  xlabel('No. unique layer interfaces');
  set(gca,'XLim',[min(A(:,4))-1 max(A(:,4))+1],'TickDir','out');
  if(isyn == 1);
    plot([nuniquetru nuniquetru],[0 1],'-w');
    plot([nuniquetru nuniquetru],[0 1],'--k');
  end;

  h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
  hold on; box on;set(gca,'FontSize',14);
  [n1,lim]=hist(A(:,3),[min(A(:,3)):max(A(:,3))]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
  n1 = n1/sum(n1);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n1,'k');
  clear n1 lim;
  xlabel('No. parameters');
  set(gca,'XLim',[min(A(:,3))-1 max(A(:,3))+1],'TickDir','out');
  if(isyn == 1);
    plot([sum(sum(voroidxtru)) sum(sum(voroidxtru))],[0 1],'-w');
    plot([sum(sum(voroidxtru)) sum(sum(voroidxtru))],[0 1],'--k');
  end;

  h4 = subplot('Position',[loc(1,4) loc(2,4) spw sph]);
  hold on; box on;set(gca,'FontSize',14);
  [n1,lim]=hist(gv,[1:NPV]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
  n1 = n1/sum(n1);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n1,'k');
  clear n1 lim;
  xlabel('No. parameters per node');
  set(gca,'XLim',[0 NPV],'TickDir','out');

end;
%%
%% Save MAP for max(P(k)) to file for starting new inversion
%%
kmin = min(k)
kmax = max(k)
if(kmax > 0);
  for i = 1:kmax-kmin+1;
    ik = kmin+i-1;
    idx = find(A(:,4) == ik);
    nmod_k(i) = length(idx);
    [a,b] = max(A(idx,1));
    mapk(i).par=A(idx(b),5:4+A(idx(b),4)*4+3);
    clear idx;
  end
  [a,b] = max(nmod_k);
  kpeak = b+kmin - 1
  idx = find(A(:,4) == kpeak);
  [a,b] = max(A(idx,1));map=A(idx(b),4:end-6);
  clear idx;
else
  [a,b] = max(A(:,1));map=A(b,4:end-6);
end;
save(mapfile,'map','-ascii');

if(iar == 1)
   (size(B,2)-6-order*NBAND-NMISC+1:size(B,2)-NMISC-6)
   alpha = B(:,end-6-(order*NBAND)-NMISC+1:end-NMISC-6);
end
if(ibucking > 0)
   (size(B,2)-6-NMISC-order*NBAND+1:size(B,2)-6-NMISC)
   bulkrho = B(:,end-6-(NMISC)+1:end-6);
end

logL = A(:,1);
if(imap == 1);
   [logLmap,jmap] = max(logL);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

if(kmax > 0);
  h = zeros(size(A,1),max(A(:,4)));
  dep = zeros(size(A,1),max(A(:,4)));
  for i = 1:size(A,1);
    if(A(:,4)>0);
      idxh = [0:k(i)-1]*NPV+1;
      dep(i,1:k(i)) = m(i,idxh);
      h(i,1:k(i)) = diff([0,dep(i,1:k(i))]);
      clear idxh;
    end;
  end;
end;
save tmp.mat h dep;
if(isd == 1)
   sd(:,1:NSD) = A(:,end-6-(order*NBAND)-NSD-NMISC+1:end-6-(order*NBAND)-NMISC);
   disp('Done getting sigma.');
end

%stop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
% fig1=figure;
% NTH = max(A(:,end));
% NPT = max(A(:,end-1));
% nx = ceil(sqrt(NTH*NPT))
% ny = nx
% xim = 0.01;
% yim = 0.06;
% xymarg = [0.1 0.04 0.04 0.14];
% [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
% 
% ith = 1;
% ipt = 1;
% jj=0;
% for i=1:NTH;
% for j=1:NPT;
%    jj=jj+1;
%    h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
%    hold on; box off;
%    set(gca,'FontSize',14);
% 
%    clear idx1 idx2;
%    if(find(A(:,end)==i));
%    idx1 = find(A(:,end)==i);
%    idx2 = find(A(idx1,end-1)==j);
%    if(ivoro == 0);
%      [AX,H1,H2] = plotyy([1:length(idx2)],...
%                 A(idx1(idx2),1),[1:length(idx2)],...
%                 A(idx1(idx2),4));
%    else
%      [AX,H1,H2] = plotyy([1:length(idx2)],...
%                 B(idx1(idx2),1),[1:length(idx2)],...
%                 B(idx1(idx2),4));
%    end;
%    set(AX(1),'YTick',[0:40:1000]);
%    if(rem(i-1,nx)==0);
%       set(get(AX(1),'Ylabel'),'String','logL') 
%       set(AX(1),'YTickLabel',[0:40:1000]);
%    else;
%       set(AX(1),'YTickLabel',[]);
%    end;
%    set(AX(2),'YTick',[0:2:100]);
%    if(rem(i,nx)==0);
%       set(get(AX(2),'Ylabel'),'String','No. interfaces') 
%       set(AX(2),'YTickLabel',[0:2:100]);
%    else;
%       set(AX(2),'YTickLabel',[]);
%    end;
%    if(i>NTH-((nx*ny)-NTH+1));
%       xlabel('rjMCMC step');
%    end;
%    set(AX(2),'XTickLabel',[],'YLim',[kmin-1 kmax+1]);
%    set(gca,'YLim',[logLmin logLmax]);
%    end;
% end;
% end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Posterior CHAIN PLOTS
%%
% fig11   = figure;
% NTH = max(A(:,end));
% NPT = max(A(:,end-1));
% nx = ceil(sqrt(NTH*NPT));
% ny = nx;
% xim = 0.01;
% yim = 0.06;
% xymarg = [0.1 0.04 0.04 0.14];
% [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
% 
% ith = 1;
% ipt = 1;
% jj=0;
% for i=1:NTH;
% for j=1:NPT;
%    jj=jj+1;
%    h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
%    hold on; box off;
%    set(gca,'FontSize',14);
% 
%    clear idx1 idx2;
%    if(find(A(:,end)==i));
%    idx1 = find(A(:,end)==i);
%    idx2 = find(A(idx1,end-1)==j);
%    plot([1:length(idx2)],A(idx1(idx2),1)+logP(idx1(idx2)));
%    set(gca,'YTick',[0:40:1000]);
%    if(rem(i-1,nx)==0);
%       set(gca,'YTickLabel',[0:40:1000]);
%    else;
%       set(gca,'YTickLabel',[]);
%    end;
%    if(i>NTH-((nx*ny)-NTH+1));
%       xlabel('rjMCMC step');
%    end;
%    %set(gca,'YLim',[logLmin+max(logP) logLmax+min(logP)]);
%    end;
% end;
% end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
if(isd == 1)
   fig3=figure;
   nx = NSD;
   ny = 2;
   xim = 0.01;
   yim = 0.06;
   xymarg = [0.1 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NSD;
      h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14);
      h2 = subplot('Position',[loc(1,i+NSD) loc(2,i+NSD) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14);
      subplot(h1);hold on;box on;
      plot([1:thinstep:length(sd(:,i))],sd(1:thinstep:end,i),'k');
      if(i == 1);ylabel('Data error standard deviation');end;
      xlabel('rjMCMC step');
      set(gca,'XLim',[0 length(sd(:,i))])
      if(imagscale == 1);set(gca,'YLim',[0.0 0.30],'YTick',[0.0:0.05 0.3]);end;
      if(imagscale == 0);set(gca,'YLim',[0.0 0.10],'YTick',[0.0:0.02 0.1]);end;
      if(i > 1);set(gca,'YTickLabel',[]);end;

      subplot(h2);set(gca,'Layer','top');hold on;
      [n1,lim]=hist(sd(:,i),100);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
      n1 = n1/sum(n1);
      lim = lim - (lim(3)-lim(2))/2;
      [xx,yy]=stairs(n1,lim,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(n1,lim,'k');
      clear n1 lim;
      if(i == 1);ylabel('Data error standard deviation');end;
%      xlabel('Probability');
      if(imagscale == 1);set(gca,'YLim',[0.0 0.30],'YTick',[0.0:0.05 0.3]);end;
      if(imagscale == 0);set(gca,'YLim',[0.0 0.10],'YTick',[0.0:0.02 0.1]);end;
      if(i > 1);set(gca,'YTickLabel',[]);end;
      box on;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Scatter plot porosity vs grain-shearing modulus
%%
fig333=figure;
b_len = sum(k+1);
bb = zeros(b_len,1);
gg = zeros(b_len,1);
l = 0;
for i=1:length(k)
  for j=1:k(i)+1;
    l = l + 1;
    ipar = (j-1)*NPL+2;
    if(j == k(i)+1); ipar = ipar - 1;end;
    bb(l) = m(i,ipar);
    gg(l) = m(i,ipar+1);
  end;
end


%x=load('fort.77');
%hold on; box on;
%set(gca,'FontSize',14);
%plot(x(:,2),log10(x(:,1)),'.');
%plot(bb(:),gg(:),'.k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Bulk Moduli and Density (pore fluid/grain material) PLOT
%%
if(ibucking > 0)
   fig33=figure;
   nx = ceil(NMISC/2);
   ny = 2;
   xim = 0.01;
   yim = 0.10;
   xymarg = [0.1 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NMISC;
      h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
      hold on; box off;
      set(gca,'FontSize',14,'Layer','top');
      miscmead(i) = median(bulkrho(:,i));
      [n,lim]=hist(bulkrho(:,i),100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
      n = n/sum(n);
      lim = lim - (lim(3)-lim(2))/2;
      [xx,yy]=stairs(lim,n,'k');
      patch(xx,yy,[0.8,0.8,0.8]);
      stairs(lim,n,'k');
      if(isyn == 1);
        plot([mtrumisc(i) mtrumisc(i)],[0 .04],'-w');
        plot([mtrumisc(i) mtrumisc(i)],[0 .04],'--k');
      end;
      clear n lim;
      set(gca,'XLim',[miscmin(i) miscmax(i)])
      if(i == 1);xlabel('Kg (mineral grains) [Pa]');end;
      if(i == 2);xlabel('Kf (interstitial fluid) [Pa]');end;
      if(i == 3);xlabel('Density grains (kg/m^3)');end;
      if(i == 4);xlabel('Density fluid (kg/m^3)');end;
      if(i == 5);xlabel('Strain hardening');end;
%      set(gca,'YLim',[0.0 0.06],'YTick',[0.0 0.02 0.04 0.06])
      set(gca,'YTickLabel',[]);
      box on;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ALPHA PLOT
%%
if(iar == 1)
   figw = 12;
   figh = 8;
   fig2=figure('visible','on');
   set(fig2,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = NBAND;
   ny = 3;
   xim = 0.01;
   yim = 0.01;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   j = 1;
   for i=1:NBAND
     h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);
     h2 = subplot('Position',[loc(1,i+NBAND) loc(2,i+NBAND) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);

     subplot(h1);hold on;box on;
     [n,lim]=hist(alpha(:,(i-1)*order+j),60);n = [0, n, 0];lim = [lim(1) lim lim(end)];
     n = n/sum(n);
     lim = lim - (lim(3)-lim(2))/2;
     [xx,yy]=stairs(n,lim,'k');
     patch(xx,yy,[0.8,0.8,0.8]);
     stairs(n,lim,'k');
     clear n lim;
     if(j==order);xlabel('Probability');end;
     if(i==1);ylabel('AR coefficient');end;
     if(i>1);set(gca,'YTickLabel',[]);end;
     if(j<order);set(gca,'XTickLabel',[]);end;
     set(gca,'XLim',[0 0.22],'XAxisLocation','top');
     set(gca,'YLim',[-.6 1]);

     idxar = ones(size(alpha(:,(i-1)*order+j)));
     idx = find(alpha(:,(i-1)*order+j) < -0.5);
     idxar(idx) = 0;
     subplot(h2);hold on;box on;
     lim = [-1:1:2];
     [n]=hist(idxar,lim);%n = [0, n, 0];lim = [lim(1) lim lim(end)];
     n = n/sum(n);
     lim = lim - (lim(3)-lim(2))/2;
     [xx,yy]=stairs(lim,n,'k');
     patch(xx,yy,[0.8,0.8,0.8]);
     stairs(lim,n,'k');
     clear n lim;
     if(i==1);ylabel('Prob. AR(1) off/on');end;
     xlabel('');
     if(i>1);set(gca,'YTickLabel',[]);end;
     xticklabel = ({'off','on'});
     set(gca,'XTick',[0,1],'XTickLabel',xticklabel);
     set(gca,'XLim',[-1 2])
     set(gca,'YLim',[0 1.1])
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT
%%
fig4=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
if(ivoro == 0);
  plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
  set(gca,'XLim',[0 length(k)],'TickDir','out');
  set(gca,'YLim',[min(k)-1 max(k)+1],'TickDir','out');
else;
  plot([1:thinstep:length(kv)],kv(1:thinstep:end),'k')
  set(gca,'XLim',[0 length(kv)],'TickDir','out');
  set(gca,'YLim',[min(kv)-1 max(kv)+1],'TickDir','out');
end;
ylabel('No. nodes');
xlabel('rjMCMC step');

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14);
if(ivoro == 0);
  [n,lim]=hist(k,[0:kmax+1]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
else;
  [n,lim]=hist(kv,[0:kmax+1]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
end;
n = n/sum(n);
lim = lim - (lim(3)-lim(2))/2;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. nodes');
ylabel('Probability');
set(gca,'XLim',[0  NVMX],'TickDir','out');
set(gca,'YLim',[ 0.  1.0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% No. parameters PLOT
%%
if(ivoro >= 1);
fig13=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
plot([1:thinstep:length(nparv)],nparv(1:thinstep:end),'k')
ylabel('No. parameters');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(nparv)],'TickDir','out');
set(gca,'YLim',[min(nparv)-1 max(nparv+1)],'TickDir','out');

subplot(2,1,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(nparv,[0:max(nparv)]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim - (lim(3)-lim(2))/2;
[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(lim,n,'k');
clear n lim;
xlabel('No. parameters');
ylabel('Probability');
set(gca,'XLim',[min(nparv)-1 max(nparv)+1],'TickDir','out');
set(gca,'YLim',[ 0. 0.5]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% K PLOT per dimension
%%
if(ivoro >= 1);
fig12=figure;
figw = 12;
figh = 6;
set(fig12,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = NPV-1;
ny = 2;
xim = 0.01;
yim = 0.1 ;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for ipar=1:NPV-1;
  h1 = subplot('Position',[loc(1,ipar) loc(2,ipar) spw sph]);
  hold on; box off;
  set(gca,'FontSize',12);
  h2 = subplot('Position',[loc(1,ipar+NPV-1) loc(2,ipar+NPV-1) spw sph]);
  hold on; box off;
  set(gca,'FontSize',12);

  subplot(h1);
  set(gca,'FontSize',12);
  plot([1:thinstep:length(kvpar(:,ipar))],kvpar(1:thinstep:end,ipar),'k')
  if(ipar == 1);ylabel('No. nodes');
  else;set(gca,'YTickLabel',[]);end;
  xlabel('rjMCMC step');
  set(gca,'YLim',[0 NVMX],'TickDir','out');
  set(gca,'XLim',[0 length(k)],'TickDir','out');
  box on;

  subplot(h2);
  set(gca,'FontSize',12);
  [n,lim]=hist(kvpar(:,ipar),[0:NVMX]);n = [0, n, 0];lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;
  if(ibucking > 0);
    if(ipar == 1);xlabel('No. nodes with porosity');end;
    if(ipar == 2);xlabel('No. nodes with \gamma_P');end;
    if(ipar == 3);xlabel('No. nodes with n');end;
    if(ipar == 4);xlabel('No. nodes with \tau_P');end;
    if(ipar == 5);xlabel('No. nodes with \gamma_S');end;
  else;
    if(ipar == 1);xlabel('No. nodes with velocity');end;
    if(ipar == 2);xlabel('No. nodes with density');end;
    if(ipar == 3);xlabel('No. nodes with attenuation');end;
  end;
  if(ipar == 1);ylabel('Probability');
  else;set(gca,'YTickLabel',[]);end;
  if(isyn == 1);
    plot([sum(voroidxtru(:,ipar+1)) sum(voroidxtru(:,ipar+1))],[0 1],'-w');
    plot([sum(voroidxtru(:,ipar+1)) sum(voroidxtru(:,ipar+1))],[0 1],'--k');
  end;
  set(gca,'XLim',[0 NVMX],'TickDir','out');
  set(gca,'YLim',[ 0. 0.6]);
  box on;
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL PLOT
%%
fig5=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14);
plot(logL,'k');
ylabel('log Likelihood');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(logL)])
set(gca,'YLim',[logLmin logLmax])
subplot(1,2,2);hold on;box on;
set(gca,'FontSize',14);
[n,lim]=hist(logL,100);n = [0, n, 0];lim = [lim(1) lim lim(end)];
n = n/sum(n);
lim = lim - (lim(3)-lim(2))/2;
[xx,yy]=stairs(n,lim,'k');
patch(xx,yy,[0.8,0.8,0.8]);
stairs(n,lim,'k');
clear n lim;
xlabel('Probability');
set(gca,'YLim',[logLmin logLmax])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% COMPUTE PROFILE MARGINALS
%%
if(imarg == 1)
   NZ = 300;
   NZI= 500;
   nsmooth = ceil(NZ/80.);
   NC = 200;
   NR = 150;
   NA = 150;
   NA2 = 100;
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
   if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
     tlim = pmin(4)+cumsum((pmax(4)-pmin(4))/NA2*ones(1,NA2));
     if(ishear == 1);
     tlim2 = pmin(5)+cumsum((pmax(5)-pmin(5))/NR*ones(1,NR));end;
   end;

   dz  = hmax/(NZ-1);
   z   = cumsum(dz*ones(1,NZ))-dz;
   dzi = hmax/(NZI-1);
   zi  = cumsum(dzi*ones(1,NZI))-dzi;

   h = zeros(1,NZI);
   nlo = zeros(NZ,1);
   nhi = zeros(NZ,1);
   c = zeros(NPROF,NZ);
   r = zeros(NPROF,NZ);
   a = zeros(NPROF,NZ);
   if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
     t = zeros(NPROF,NZ);
     if(ishear == 1);t2 = zeros(NPROF,NZ);end;
   end;
   disp('Sample size: '),disp(size(m))
   for iprof = 1:NPROF

     if(rem(iprof,10000)==0)
        fprintf(1,'done %8i samples\n',iprof)
     end
     clear idxh idxc idxr idxa prof;
     %% Find index for current model
     if(k(iprof) > 0)
        idxh = (([1:k(iprof)]-1)*NPL)+1;
        idxc = (([1:k(iprof)]-1)*NPL)+2;
        idxr = (([1:k(iprof)]-1)*NPL)+3;
        idxa = (([1:k(iprof)]-1)*NPL)+4;
        if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
          idxt = (([1:k(iprof)]-1)*NPL)+5;
          if(ishear == 1);
            idxt2 = (([1:k(iprof)]-1)*NPL)+6;
          end;
        end;
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+NPL-1];
        idxr = [idxr idxr(end)+NPL-1];
        idxa = [idxa idxa(end)+NPL-1];
        if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
          idxt = [idxt idxt(end)+NPL-1];
          if(ishear == 1);
            idxt2 = [idxt2 idxt2(end)+NPL-1];
          end;
        end;
     else
        idxh = [];
        idxc = [1];
        idxr = [2];
        idxa = [3];
        if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
          idxt = [4];
          if(ishear == 1);idxt2 = [5];end;
        end;
     end

     %% Compute the profile for current model
     if(k(iprof) > 0)
       prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
       prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
       prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
       prof(:,2) = m(iprof,idxc);
       prof(:,3) = m(iprof,idxr);
       prof(:,4) = m(iprof,idxa);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
           prof(:,5) = m(iprof,idxt);
         if(ivpvs == 1 & ibucking == 0);
           prof(:,5) = m(iprof,idxc)./exp(m(iprof,idxt));
         end;
         if(ishear == 1);prof(:,6) = m(iprof,idxt2);end;
       end;

       for ilay=2:k(iprof)+1  %% k is # layers of current model
          idxzi = round(prof(ilay-1,1)/dzi);
          if(idxzi == 0);idxzi = 1;end;
          h(idxzi) = h(idxzi) + 1;
       end;
       c(iprof,:) = prof(1,2);
       r(iprof,:) = prof(1,3);
       a(iprof,:) = prof(1,4);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         t(iprof,:) = prof(1,5);
         if(ishear == 1);t2(iprof,:) = prof(1,6);end;
       end;
       for ilay=2:k(iprof)+1  %% k is # layers of current model
          idxz = round(prof(ilay-1,1)/dz);
          if(idxz == 0)idxz = 1;end;
          c(iprof,idxz:end) = prof(ilay,2);
          r(iprof,idxz:end) = prof(ilay,3);
          a(iprof,idxz:end) = prof(ilay,4);
          if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
            t(iprof,idxz:end) = prof(ilay,5);
            if(ishear == 1);t2(iprof,idxz:end) = prof(ilay,6);end;
          end;
       end;

     else
       prof(:,2) = m(iprof,idxc);
       prof(:,3) = m(iprof,idxr);
       prof(:,4) = m(iprof,idxa);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
           prof(:,5) = m(iprof,idxt);
         if(ivpvs == 1 & ibucking == 0);
           prof(:,5) = m(iprof,idxc)./exp(m(iprof,idxt));
         end;
         if(ishear == 1);prof(:,6) = m(iprof,idxt2);end;
       end;
       c(iprof,:) = prof(1,2);
       r(iprof,:) = prof(1,3);
       a(iprof,:) = prof(1,4);
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         t(iprof,:) = prof(1,5);
         if(ishear == 1);t2(iprof,:) = prof(1,6);end;
       end;
     end
   end;
   fprintf(1,'\n')
   disp('Done with profiles.');

   %
   % Compute histograms for each depth
   %
   disp('Starting histograms...');
   for iz=1:NZ
      [Nc(iz,:),binsc] = hist(c(:,iz),clim);
      [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
      [Na(iz,:),binsa] = hist(a(:,iz),alim);
      if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
        [Nt(iz,:),binst] = hist(t(:,iz),tlim);
        if(ishear == 1);[Nt2(iz,:),binst2] = hist(t2(:,iz),tlim2);end;
      end;
      if(iboudreau == 1);
        nlo(iz) = sqrt(1-0.36*log(max(c(:,iz))^2))-0.06;
        nhi(iz) = sqrt(1-1.6*log(min(c(:,iz))^2))+0.2;
      end;
      [nfc(iz,:)] = hpd(c(:,iz),100,95);
      [nfr(iz,:)] = hpd(r(:,iz),100,95);
      [nfa(iz,:)] = hpd(a(:,iz),100,95);
      if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
        [nft(iz,:)] = hpd(t(:,iz),100,95);
        if(ishear == 1);[nft2(iz,:)] = hpd(t2(:,iz),100,95);end;
      end;
      meac(iz) = median(c(:,iz));
      mear(iz) = median(r(:,iz));
      meaa(iz) = median(a(:,iz));
      if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
        meat(iz) = median(t(:,iz));
        if(ishear == 1);meat2(iz) = median(t2(:,iz));end;
      end;
   end;
   %
   % Normalize Histograms
   %
   if(inorm == 0)
     for iz=1:NZ
       Nc(iz,:) = Nc(iz,:)/NPROF;
       Nr(iz,:) = Nr(iz,:)/NPROF;
       Na(iz,:) = Na(iz,:)/NPROF;
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         Nt(iz,:) = Nt(iz,:)/NPROF;
         if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/NPROF;end;
       end;
     end;
   elseif(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:));
         Na(iz,:) = Na(iz,:)/trapz(binsa,Na(iz,:));
        if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
          Nt(iz,:) = Nt(iz,:)/trapz(binst,Nt(iz,:));
          if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/trapz(binst2,Nt2(iz,:));end;
        end;
      end;
   elseif(inorm == 2)
     for iz=1:NZ
       Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
       Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:));
       Na(iz,:) = Na(iz,:)/max(Na(iz,:));
       if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
         Nt(iz,:) = Nt(iz,:)/max(Nt(iz,:));
         if(ishear == 1);Nt2(iz,:) = Nt2(iz,:)/max(Nt2(iz,:));end;
       end;
     end;
   end;
   disp('Done histograms.');

%if(ibucking == 3 & ivgsnosh == 0);
%save tmp.mat z nfc nfa nfr nft meac mear meaa meat miscmead;
%else;
%save tmp.mat z nfc nfa nfr meac mear meaa miscmead;
%end;

c_mean = mean(c);
r_mean = mean(r);
a_mean = mean(a);
if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
  t_mean = mean(t);
end;
c_mead = median(c);
r_mead = median(r);
a_mead = median(a);
if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));
  t_mead = median(t);
end;

[ntmp,idxcmax] = max(Nc,[],2);
[ntmp,idxrmax] = max(Nr,[],2);
[ntmp,idxamax] = max(Na,[],2);
if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));[ntmp,idxtmax] = max(Nt,[],2);end;
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   r_max(iz) = rlim(idxrmax(iz));
   a_max(iz) = alim(idxamax(iz));
   if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));t_max(iz) = tlim(idxtmax(iz));end;
end;
NAVE = 8;
for i=1:length(c_max)
   if(i <= NAVE)
       NAVE1 = 1;
   else
       NAVE1 = i-NAVE;
   end
   if(i >= length(c_max)-NAVE)
       NAVE2 = length(c_max);
   else
       NAVE2 = i+NAVE;
   end
   c_max_sm(i) = sqrt(mean(c_max(NAVE1:NAVE2).^2));
   r_max_sm(i) = sqrt(mean(r_max(NAVE1:NAVE2).^2));
   a_max_sm(i) = sqrt(mean(a_max(NAVE1:NAVE2).^2));
   if((ibucking == 3 & ivgsnosh == 0) | (ibucking == 0 & ishear == 1));t_max_sm(i) = sqrt(mean(t_max(NAVE1:NAVE2).^2));end;
end

end; % end imarg

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
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
%[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

if((ibucking==3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
  spw1 = 0.06;
  spw2 = 0.2;
  sph = 0.82;
  loc(:,1) = [0.06; 0.14];
  loc(:,2) = [loc(1,1)+spw1+.01; 0.14];
  loc(:,3) = [loc(1,2)+spw2+.01; 0.14];
  loc(:,4) = [loc(1,3)+spw2+.01; 0.14];
  loc(:,5) = [loc(1,4)+spw2+.01; 0.14];  
  if(ishear==1);
    loc = [0.075, 0.163, 0.331, 0.499, 0.667, 0.835;...
           0.14,   0.14,  0.14,  0.14,  0.14, 0.14];
    spw1 = 0.08;
    spw2 = 0.16;
    sph  = 0.82;
  end;

else
  loc = [0.08,  0.18, 0.4533, 0.7266;...
         0.14, 0.14, 0.14, 0.14];
  spw1 = 0.09;
  spw2 = 0.2633;
  sph = 0.82;
end;

fig6 = figure;hold on; box on;
figw = 18;
figh = 10;
fig6.PaperUnits = 'inches';
fig6.PaperPosition = [0 0 figw figh];
set(fig6, 'renderer', 'painters')
h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);
hold on; box off;
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;
if((ibucking == 3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
  h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]);
  hold on; box off;
  if(ishear == 1);
    h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph]);
    hold on; box off;
  end;
end;

subplot(h4);hold on;box on;
if(imarg == 1)
   h = h/sum(h);
   [xx,yy]=stairs(h,zi,'-k');
   xx = [0;xx;0];
   yy = [yy(1);yy;yy(end)];
   patch(xx,yy,[0.8,0.8,0.8]);
   [xx,yy]=stairs(h,zi,'-k');
   set(gca,'Fontsize',14,'YLim',[0 hmax],'XLim',[0 0.1]);
   set(gca,'YDir','reverse','TickDir','out');
   xlabel('Interface probability');
   ylabel('Depth (m)');
   box on;
end;

subplot(h1)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
   if(imead == 1);plot(c_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(c_mean,z,'.k','Linewidth',2);end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
if(isyn == 1)
   h1=plprof(mtru,hmax,NPL,'-w',1,0);
   set(h1,'LineWidth',2);
   h1=plprof(mtru,hmax,NPL,'--k',1,0);
   set(h1,'LineWidth',2);
end;
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',1,0);
end;
if(imax == 1)
  for i=1:length(z);
     [cmx,j] = max(Nc(i,:));
     c_max(i) = clim(j);
  end;
  plot(c_max,z,'w');
end;
if(ibucking > 0);
  set(gca,'XTick',[0.25 0.5 0.75]);
  xlabel('Porosity');
else; 
  set(gca,'XTick',[1450:50 :4000]);
  xlabel('Comp. Vel. (m/s)');
end;
box on;

subplot(h2)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(rlim,z,Nr);shading flat;
   if(imead == 1);plot(r_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(r_mean,z,'.k','Linewidth',2);end;
end;
if(iboudreau == 1);plot(nlo,z,'-w');end;
if(iboudreau == 1);plot(nhi,z,'-w');end;
%surf(rlim,z,Nr);shading flat;
set(h2,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
if(ibucking == 0);
  xlabel('Density (g/ccm)');
  set(gca,'XTick',[1.3:.2:3]);
elseif(ibucking == 1);
  xlabel('Tortuosity');
elseif(ibucking == 3);
  xlabel('log10(\gamma_P)');
  set(gca,'XTick',[1:2:26]);
end;
if(isyn == 1)
   h1=plprof(mtru,hmax,NPL,'-w',2,0);
   set(h1,'LineWidth',2);
   h1=plprof(mtru,hmax,NPL,'--k',2,0);
   set(h1,'LineWidth',2);
end
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',2);
end
if(imax == 1)
  for i=1:length(z);
     [rmx,j] = max(Nr(i,:));
     r_max(i) = rlim(j);
  end;
  plot(r_max,z,'w');
end;
box on;

subplot(h3)
set(gca,'FontSize',14);
if(imarg == 1)
   pcolor(alim,z,Na);shading flat;
   if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
   if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
end;
%surf(alim,z,Na);shading flat;
set(h3,'layer','top')
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
if(ibucking == 0);
  xlabel('log_{10} \alpha_P (dB/m/kHz)');
  %set(gca,'XTick',[0:0.3:1]);
  set(gca,'XTick',[-10:1:10]);
elseif(ibucking == 1);
  xlabel('log Transition freq.');
elseif(ibucking == 2);
  xlabel('strain hardening idx');
elseif(ibucking==3 & ivgsnosh==0);
  xlabel('log10(Mat. idx)');
  set(gca,'XTick',[-10:1:2]);
elseif(ibucking==3 & ivgsnosh==1);
  xlabel('log10(\tau_P)');
end;
if(isyn == 1)
   h1=plprof(mtru,hmax,NPL,'-w',3,0);
   set(h1,'LineWidth',2);
   h1=plprof(mtru,hmax,NPL,'--k',3,0);
   set(h1,'LineWidth',2);
end;
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',3);
end;
if(imax == 1)
  for i=1:length(z);
     [amx,j] = max(Na(i,:));
     a_max(i) = alim(j);
  end;
  plot(a_max,z,'w');
  save('max_model.mat','z','c_max','r_max','a_max');
end;
cmap = colormap(jet);
colormap(cmap);
box on;

if((ibucking==3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
  subplot(h5)
  set(gca,'FontSize',14);
  if(imarg == 1)
    pcolor(tlim,z,Nt);shading flat;
    if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
    if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
  end;
  set(h5,'layer','top')
  set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax]);
  set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
  cmap = colormap(jet);
  colormap(cmap);
  box on;
  if(isyn == 1)
    h1=plprof(mtru,hmax,NPL,'-w',4,0);
    set(h1,'LineWidth',2);
    h1=plprof(mtru,hmax,NPL,'--k',4,0);
    set(h1,'LineWidth',2);
  end;
  if(ibucking == 0);
    xlabel('V_S (m/s)');
    set(gca,'XTickLabel',[0:750:2500]);
    set(gca,'XTick',[0:750:2500]);
  elseif(ibucking == 3);
    xlabel('log10(\tau_P)');
    set(gca,'XTick',[-15:1:0]);
  end;

  if(ishear==1);
    subplot(h6)
    set(gca,'FontSize',14);
    if(imarg == 1)
      pcolor(tlim2,z,Nt2);shading flat;
      if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
      if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
    end;
    set(h6,'layer','top')
    set(gca,'Fontsize',14,'XLim',[pmin(5) pmax(5)],'YLim',[0 hmax]);
    set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
    cmap = colormap(jet);
    colormap(cmap);
    box on;
    if(ibucking == 0);
      xlabel('\alpha_S (dB/m/kHz)');
      set(gca,'XTickLabel',[0:1:5]);
      set(gca,'XTick',[0:1:5]);
    elseif(ibucking == 3);
      xlabel('log(\gamma_S)');
      set(gca,'XTickLabel',[0:5:25]);
      set(gca,'XTick',[0:5:25]);
    end;
    if(isyn == 1)
      h1=plprof(mtru,hmax,NPL,'-w',5,0);
      set(h1,'LineWidth',2);
      h1=plprof(mtru,hmax,NPL,'--k',5,0);
      set(h1,'LineWidth',2);
    end;
  end;
end;

if(icore == 1)
   barstep1 = 5;
   barstep1r = 5;
   barstep2 = 10;
   barstep2r = 1;
   barstep3 = 10;

   subplot(h1);
   plot(c1(:,2),c1(:,1),'-w','Linewidth',1);
%   errorbarxy(c1(1:barstep1:end,2),c1(1:barstep1:end,1), ...
%              10.*ones(size(c1(1:barstep1:end,2))),...
%              zeros(size(c1(1:barstep1:end,2))),'w','w');
   plot(c2(:,2),c2(:,1),'-g','Linewidth',1);
%   errorbarxy(c2(1:barstep2:end,2),c2(1:barstep2:end,1), ...
%             10.*ones(size(c2(1:barstep2:end,2))),...
%              zeros(size(c2(1:barstep2:end,2))),'g','g');
   plot(c3(:,2),c3(:,1),'-r','Linewidth',1);
%   errorbarxy(c3(1:barstep3:end,2),c3(1:barstep3:end,1), ...
%              10.*ones(size(c3(1:barstep3:end,2))),...
%              zeros(size(c3(1:barstep3:end,2))),'r','r');
   %plot(c4(:,2),c4(:,1),'c','Linewidth',1);
   subplot(h2);
   plot(r1(:,2),r1(:,1),'-w','Linewidth',1);
%   errorbarxy(r1(1:barstep1r:end,2),r1(1:barstep1r:end,1), ...
%              2./100.*r1(1:barstep1r:end,2),...
%              zeros(size(r1(1:barstep1r:end,2))),'w','w');
   plot(r2(:,2),r2(:,1),'-g','Linewidth',1);
%   errorbarxy(r2(1:barstep2r:end,2),r2(1:barstep2r:end,1), ...
%              2./100.*r2(1:barstep2r:end,2),...
%              zeros(size(r2(1:barstep2r:end,2))),'w','w');
   plot(r3(:,2),r3(:,1),'-r','Linewidth',1);
%   errorbarxy(r3(1:barstep3:end,2),r3(1:barstep3:end,1), ...
%              2./100.*r3(1:barstep3:end,2),...
%              zeros(size(r3(1:barstep3:end,2))),'r','r');
end;

%%
%%  Node density of voronoi cells:
%%
if(ivoro >= 1);
fig10 = figure;hold on; box on;
set(fig10, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,4) loc(2,4) spw2 sph]);
hold on; box off;
if((ibucking == 3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
  h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]);
  hold on; box off;
  if(ishear == 1);
    h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph]);
    hold on; box off;
  end;
end;

subplot(h1)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax]);
[hdens]=cloudPlot(por(:,2),por(:,1),[pmin(1) pmax(1) 0 hmax],true,[1001 201]);
colormap( jet(1000) );
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out');
set(gca,'CLim',[0 1.5],'FontSize',14)
if(ibucking >0);
  xlabel('Porosity');
  set(gca,'XTick',[0.25 0.5 0.75]);
else;
  xlabel('Comp. vel. (m/s)');
  set(gca,'XTick',[1450:50:3000]);
end;
ylabel('Depth (m)');
box on;

subplot(h2)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax]);
[hdens]=cloudPlot(Kp(:,2),Kp(:,1),[pmin(2) pmax(2) 0 hmax],true,[1001 201]);
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(gca,'CLim',[0 1.5],'FontSize',14)
if(ibucking >0);
  xlabel('log(K_P)');
  set(gca,'XTick',[1:2:26]);
else;
  xlabel('Density (g/cm^3)');
  set(gca,'XTick',[1.3:0.2:3.0]);
end;
box on;

subplot(h3)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
[hdens]=cloudPlot(strain(:,2),strain(:,1),[pmin(3) pmax(3) 0 hmax],true,[1001 201]);
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(gca,'CLim',[0 1.5],'FontSize',14)
if(ibucking >0);
  xlabel('material index');
  set(gca,'XTick',[-10:1:2]);
else;
  xlabel('log_{10} \alpha_P (dB/m/kHz)');
  set(gca,'XTick',[-10:1:2]);
end;
box on;

if((ibucking == 3 & ivgsnosh==0) | (ibucking == 0 & ishear == 1));
subplot(h5)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax]);
[hdens]=cloudPlot(taup(:,2),taup(:,1),[pmin(4) pmax(4) 0 hmax],true,[1001 201]);
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(gca,'CLim',[0 1.5],'FontSize',14)
xlabel('log(\tau_P)');
set(gca,'XTick',[-15:1:0]);
box on;
end;

if(ishear == 1);
subplot(h6)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(5) pmax(5)],'YLim',[0 hmax]);
[hdens]=cloudPlot(Ks(:,2),Ks(:,1),[pmin(5) pmax(5) 0 hmax],true,[1001 201]);
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(gca,'CLim',[0 1.5],'FontSize',14)
xlabel('log(Ks)');
set(gca,'XTick',[0:5:25]);
box on;
end;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  COMPUTE DATA FIT
%%
if(idatfit == 1)

   for i=1:length(bands);
      stat(i).idxnan  = find(rex(:,i) == 0);
      stat(i).idxgood = find(rex(:,i) ~= 0);
      stat(i).rnan = rex(:,i);
      stat(i).rnan(stat(i).idxnan) = NaN;
   end;
   size(dobs)
   size(rex)
   dobs = dobs .* rex;
   for iband=1:length(bands);
     if(NAVEF > 1)
       if(IGA == 1)
         flo = bands(iband) - bands(iband)*frbw;
         fhi = bands(iband) + bands(iband)*frbw;
       elseif(IGA == 2)
         flo = bands(iband) / (2.^(1./(2.*oct))); % Use 1/3 octave
         fhi = bands(iband) * (2.^(1./(2.*oct))); %
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
   end
   NFREQ = length(fr);

   disp('load ensemble files');
   restmp=dlmread(resfile);
   if(iar == 1);resartmp=dlmread(resarfile);end;
   reptmp=dlmread(repensfile);
   NDAVE = length(reptmp)/NBAND

   ipass = zeros(NBAND,1);
   ifail = zeros(NBAND,1);
   for j=1:NDAVE;
     if(rem(j,100)==0)
        fprintf(1,'%8i',j);
     end
     for iband=1:length(bands);

       ref(:,iband,j) = reptmp(NBAND*(j-1)+iband,:);
       res(:,iband,j) = restmp(NBAND*(j-1)+iband,:);
       res(:,iband,j) = res(:,iband,j)-mean(res(:,iband,j));

       %% Do runs test
       [h(j),p(j)]=runstest(res(:,iband,j),0.);
       if(h(j) == 1);ifail(iband) = ifail(iband)+1;end;
       if(h(j) == 0);ipass(iband) = ipass(iband)+1;end;

       stat(iband).res(:,j) = res(:,iband,j);
       [stat(iband).axx(:,j),stat(iband).axxlags(:,j)] = ...
            xcorr(stat(iband).res(:,j),'coeff');
     %%if(iband==12);save tmp.mat refl fr r_ave angobs;end;
     end;
   end;

   ipass'
   ifail'
   ipass'+ifail'
   ipass'./(ipass'+ifail')
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT DATA MISFIT
   %%
   nx = 3;
   ny = ceil(length(bands)/nx);
   if(length(bands)>12);ny = 3;end;
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig7 = figure;hold on; box on;
   if(length(bands)>12);fig77 = figure;hold on; box on;end;
   ii = 1;
   diffang = diff(angobs)/2;
   angplt = [angobs(1)-2*diffang(1),angobs-[diffang,diffang(end)],...
             angobs(end)+2*diffang(end)];
   for i=1:length(bands);
      figure(fig7);
      jj = i;
      if(i>12);
        figure(fig77);
        jj=i-12;
      end;
      subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      refmean = mean(ref,3);
%      set(gca,'XLim',[angobs(1)-2 angobs(end)+2]);
      set(gca,'XLim',[angplt(1) angplt(end)]);
      set(gca,'FontSize',14);
      set(gca,'LineWidth',1);
      if(idathist == 1);
         NV = 800;
         vmin = 0.0;
         vmax = 1.1;
         vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
         for iang=1:NANG
            Nd(:,iang) = histc(ref(iang,i,:),vlim);
         end;
         %
         % Normalize Histograms
         %
         if(inorm == 1)
            for iz=1:NANG
               Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
            end;
         elseif(inorm == 2)
            Nd = Nd/max(max(Nd));
         end;
         Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
         pcolor(angplt,vlim,Nd);shading flat;
         clear Nd,vlim;
         plot(angobs,dobs(:,i),'o','MarkerEdgeColor','w',...
                   'MarkerFaceColor','k',...
                   'MarkerSize',4,'Linewidth',1);
         plot(angobs,dobs(:,i),'-w','Linewidth',1.5);
         plot(angobs,dobs(:,i),'--k','Linewidth',1.5);
      else;
         for j=1:length(angobs);
                y = ref(j,i,:);
            [nf(j,:)] = hpd(y,100,95);
         end;
         for j=1:NDAVE;
            plot(angobs,ref(:,i,j),':r');
         end;
%         plot(angobs,refmean(:,i),'-r');
         plot(angobs,nf(:,1),'-w');
         plot(angobs,nf(:,1),'--k');
         plot(angobs,nf(:,2),'-w');
         plot(angobs,nf(:,2),'--k');
%         plot(angobs,dobs(:,i),'-k','Linewidth',1.5);
         plot(angobs,dobs(:,i),'xk','Linewidth',1.5);
      end;
      if(ibl == 0)
         if(max(max(dobs))> 1.)
            ylim  = [0,1.3];
            ytick = [0,0.3,0.6,0.9,1.2];
            text(60,1.,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         elseif(max(max(dobs))> 0.64)
            ylim  = [0,1.01];
            ytick = [0,0.3,0.6,0.9];
            text(60,0.85,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         else
            ylim  = [0,0.66];
            ytick = [0,0.2,0.4,0.6];
            text(60,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12,'Color',[0,0,0])
         end
         set(gca,'YLim',ylim,'XLim',[angplt(1) angplt(end)]);
         if ((i == (ii-1)*nx+1))
             set(gca,'YTickLabel',ytick,'YTick',ytick);
             ylabel('Refl. Coeff.');
             ii = ii + 1;
         else
              set(gca,'YTickLabel',[],'YTick',ytick);
         end
      else
         text(48,15,[num2str(bands(i)) ' Hz'],'FontSize',12)
         if(max(max(dobs))< 20)
            ylim  = [0,20];
            ytick = [0,5,10,15];
         else
            ylim  = [10,40];
            ytick = [10,20,30];
         end
         set(gca,'YLim',ylim);
         if ((i == 1) | (i == nx+1) | (i == (2*nx)+1))
             set(gca,'YTickLabel',ytick,'YTick',ytick);
             ylabel('Bottom Loss (dB)');
         else
             set(gca,'YTickLabel',[],'YTick',ytick);
         end
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[10:10:90],'XTick',[10:10:90]);
          xlabel('Angle (deg.)');
      else
          set(gca,'XTickLabel',[],'XTick',[10:10:90]);
      end
   end;
   clear nf;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT Axx
   %%
   nx = 4;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   fig8=figure();
   ii = 1;
   for i=1:length(bands); 
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      clear axxmean;
      axxmean = mean(stat(i).axx,2);
      plot([stat(i).axxlags(1,1) stat(i).axxlags(end,1)],[0 0],'--k');
      plot(stat(i).axxlags(:,1),axxmean,'.-k');
      for j=1:length(axxmean);
         y = stat(i).axx(j,:);
         [nf(j,:)] = hpd(y,100,98);
         clear y;
      end;
      plot(stat(i).axxlags(:,1),nf(:,1),'--k');
      plot(stat(i).axxlags(:,1),nf(:,2),'--k');
      text(stat(i).axxlags(end,1)-20,0.8,[num2str(bands(i)) ' Hz'],'FontSize',12)
      set(gca,'XLim',[stat(i).axxlags(1,1) stat(i).axxlags(end,1)],'YLim',[-.45 1.1]);
      if ((i == (ii-1)*nx+1))
          set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
          ylabel('Autocorrelation');
          ii = ii+1;
      else
           set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[-20 0 20],'XTick',[-20 0 20]);
          xlabel('Lag');
      else
          set(gca,'XTickLabel',[],'XTick',[-20 0 20]);
      end
      clear nf;
   end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT res hist 
   %%
   fig9=figure();
   x = -16.375:.75:16.375;
   xx = -16.25:.01:16.25;
   ii = 1;
   for i=1:length(bands); 
      resmean = mean(stat(i).res(stat(i).idxgood),2);
      subplot('Position',[loc(1,i) loc(2,i) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      resmean = resmean - mean(resmean);
      if(isd == 1)
         [n1,xout] = hist(resmean./mean(sd(:,i)),x);
      else
         [n1,xout] = hist(resmean./std(resmean),x);
      end
      narea = sum(n1) * (xout(2)-xout(1));
      n1 = n1/narea;
      xout = xout - (xout(2)-xout(1))/2;
      stairs(xout,n1,'k');
      nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
      plot(xx,nd,'--k')

      text(1,0.55,[num2str(bands(i)) ' Hz'],'FontSize',12)
      axis([-6 6 0 0.6]);
      set(gca,'XTick',[-4 0 4],'FontSize',14);
      set(gca,'XTickLabel',[],'FontSize',14);
      if ((i == (ii-1)*nx+1))
          set(gca,'YTickLabel',[0 0.2 0.4],'YTick',[0 0.2 0.4]);
          ii = ii+1;
      else
           set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end
      if (i > length(bands)-nx)
          set(gca,'XTickLabel',[-4 0 4],'XTick',[-4 0 4]);
          xlabel('Res. (std.dev.)');
      else
          set(gca,'XTickLabel',[],'XTick',[-4 0 4]);
      end
      clear resmean;
   end;
end;

if(isave == 1)
  if(ivoro >= 1);
    saveas(fig0,strcat(plotfile0,plotext2),'png');
    if(IEPS == 1);
    print(fig0,'-painters','-r250',strcat(plotfile0,plotext3),'-depsc');
    end;
  end;
%   saveas(fig1,strcat(plotfile1,plotext1),'fig');
%   saveas(fig1,strcat(plotfile1,plotext2),'png');
%   if(IEPS == 1);
%   print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
%   end;
   if(iar == 1)
%     saveas(fig2,strcat(plotfile2,plotext1),'fig');
     saveas(fig2,strcat(plotfile2,plotext2),'png');
   if(IEPS == 1);
     print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
   end;
   end
   if(isd == 1)
%      saveas(fig3,strcat(plotfile3,plotext1),'fig');
      saveas(fig3,strcat(plotfile3,plotext2),'png');
   if(IEPS == 1);
      print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
   end;
   end
   if(ibucking > 0)
%      saveas(fig3,strcat(plotfile3,plotext1),'fig');
      saveas(fig33,strcat(plotfile33,plotext2),'png');
   if(IEPS == 1);
      print(fig33,'-painters','-r250',strcat(plotfile33,plotext3),'-depsc');
   end;
   end
%   saveas(fig4,strcat(plotfile4,plotext1),'fig');
   saveas(fig4,strcat(plotfile4,plotext2),'png');
   if(IEPS == 1);
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');
   end;
%   saveas(fig5,strcat(plotfile5,plotext1),'fig');
   saveas(fig5,strcat(plotfile5,plotext2),'png');
   if(IEPS == 1);
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');
   end;
%   saveas(fig6,strcat(plotfile6,plotext1),'fig');
%   saveas(fig6,strcat(plotfile6,plotext2),'png');
   print(fig6,strcat(plotfile6,plotext2),'-dpng','-r0')
   if(IEPS == 1);
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');
   end;
   if(ivoro >= 1);
%    saveas(fig10,strcat(plotfile10,plotext1),'fig');
     saveas(fig10,strcat(plotfile10,plotext2),'png');
     if(IEPS == 1);
     print(fig10,'-painters','-r250',strcat(plotfile10,plotext3),'-depsc');
     end;
%    saveas(fig12,strcat(plotfile12,plotext1),'fig');
     saveas(fig12,strcat(plotfile12,plotext2),'png');
     if(IEPS == 1);
       print(fig12,'-painters','-r250',strcat(plotfile12,plotext3),'-depsc');
     end;
   end;
   if(idatfit == 1);
%      saveas(fig7,strcat(plotfile7,plotext1),'fig');
      saveas(fig7,strcat(plotfile7,plotext2),'png');
   if(IEPS == 1);
      print(fig7,'-painters','-r250',strcat(plotfile7,plotext3),'-depsc');
   end;
     if(length(bands)>12);
%       saveas(fig77,strcat(plotfile77,plotext1),'fig');
       saveas(fig77,strcat(plotfile77,plotext2),'png');
   if(IEPS == 1);
       print(fig77,'-painters','-r250',strcat(plotfile77,plotext3),'-depsc');
   end;
     end;
%      saveas(fig8,strcat(plotfile8,plotext1),'fig');
      saveas(fig8,strcat(plotfile8,plotext2),'png');
   if(IEPS == 1);
      print(fig8,'-painters','-r250',strcat(plotfile8,plotext3),'-depsc');
   end;
%      saveas(fig9,strcat(plotfile9,plotext1),'fig');
      saveas(fig9,strcat(plotfile9,plotext2),'png');
   if(IEPS == 1);
      print(fig9,'-painters','-r250',strcat(plotfile9,plotext3),'-depsc');
   end;
   end;
%   if(imarg == 1)
%      save(profilefile,'z', 'c', 'r', 'a','c_mean','r_mean','a_mean');
%   end;
end;

return;
