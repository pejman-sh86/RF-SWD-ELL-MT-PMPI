%----------------------------------------------------------
% Compute Marginals from PPD
%----------------------------------------------------------
function [] = plot_hist();
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

i_save = 1;

isd = 0;
itrans = 1;      % Transition Layer switch
IBFIX = 0;       % Fix the basement in transition layer? 0=yes, 1=no
i_scaley = 0;
Tstar = 1;
Tstar2 = 1;
itrue = 0;
nx = 4;
ny = 6;

% INPUT:
%ppd     = 'x_s21_1_5_40_4lay_sample.mat';
%data    = 'x_s21_1_5_40_4lay.mat';
%ppd     = 'x_s20_1_8_40_0lay_sample.mat';
%data    = 'x_s20_1_8_40_0lay.mat';
%ppd     = 'x_s19_1_8_25_0lay_sample.mat';
%data    = 'x_s19_1_8_25_0lay.mat';
%ppd     = 'x_s16_1_5_25_6lay_sample.mat';
%data    = 'x_s16_1_5_25_6lay.mat';
%ppd     = 'x_s13_2_10_50_5layb_sample.mat';
%data    = 'x_s13_2_10_50_5layb.mat';
%ppd     = 'tmp_sample.mat';
%ppd     = 'x_s07_1_1_100_alf_2lay_sample.mat';
%data    = 'x_s07_1_1_100_alf_2lay.mat';
%ppd     = 'x_s07_1_1_100_4lay_sample.mat';
%data    = 'x_s07_1_1_100_4lay.mat';
ppd     = 'x_s05_1_8_32_0lay_sample.mat';
data    = 'x_s05_1_8_32_0lay.mat';
%ppd     = 'x_s05_2_8_25_4lay_sample.mat';
%data    = 'x_s05_2_8_25_4lay.mat';
%ppd     = 'x_s04_1_3_16_0lay_sample.mat';
%data    = 'x_s04_1_3_16_0lay.mat';
%ppd     = 'x_s02_1_3_25_6lay_sample.mat';
%data    = 'x_s02_1_3_25_6lay.mat';
%ppd     = 'x_s01_1_3_20_5lay_sample.mat';
%data    = 'x_s01_1_3_20_5lay.mat';
%ppd     = 'sim_T_1_4_17_0lay_sample.mat';
%data    = 'sim_T_1_4_17_0lay.mat';
%ppd     = 'sim_C_1_7_40_6lay_sample.mat';
%data    = 'sim_C_1_7_40_6lay.mat';
%ppd     = 'sim_B_1_10_50_1lay_sample.mat';
%data    = 'sim_B_1_10_50_1lay.mat';
%ppd     = 'sim_A_1_3_16_2lay_sample.mat';
%data    = 'sim_A_1_3_16_2lay.mat';

% OUTPUT:
if(i_save == 1)
%   mapfile  = 'x_s21_1_5_40_4lay_fgsmap.txt';
%   meanfile = 'x_s21_1_5_40_4lay_mean.txt';
%   hpdfile  = 'x_s21_1_5_40_4lay_hpds.txt';
%   plotfile = 'x_s21_1_5_40_4lay_marg_a.png';
%   mapfile  = 'x_s20_1_8_40_0lay_fgsmap.txt';
%   meanfile = 'x_s20_1_8_40_0lay_mean.txt';
%   hpdfile  = 'x_s20_1_8_40_0lay_hpds.txt';
%   plotfile = 'x_s20_1_8_40_0lay_marg_a.png';
%   mapfile  = 'x_s19_1_8_25_0lay_fgsmap.txt';
%   meanfile = 'x_s19_1_8_25_0lay_mean.txt';
%   hpdfile  = 'x_s19_1_8_25_0lay_hpds.txt';
%   plotfile = 'x_s19_1_8_25_0lay_marg_a.png';
%   mapfile  = 'x_s16_1_5_25_6lay_fgsmap.txt';
%   meanfile = 'x_s16_1_5_25_6lay_mean.txt';
%   hpdfile  = 'x_s16_1_5_25_6lay_hpds.txt';
%   plotfile = 'x_s16_1_5_25_6lay_marg_a.png';
%   mapfile  = 'x_s13_2_10_50_5layb_fgsmap.txt';
%   meanfile = 'x_s13_2_10_50_5layb_mean.txt';
%   hpdfile  = 'x_s13_2_10_50_5layb_hpds.txt';
%   plotfile = 'x_s13_2_10_50_5layb_marg_a.png';
%   mapfile  = 'x_s07_1_1_100_alf_2lay_fgsmap.txt';
%   meanfile = 'x_s07_1_1_100_alf_2lay_mean.txt';
%   hpdfile  = 'x_s07_1_1_100_alf_2lay_hpds.txt';
%   plotfile = 'x_s07_1_1_100_alf_2lay_marg_a.png';
%   plotsdfile = 'x_s07_1_1_100_alf_2lay_marg_y.png';
   mapfile  = 'x_s05_1_8_32_0lay_fgsmap.txt';
   meanfile = 'x_s05_1_8_32_0lay_mean.txt';
   hpdfile  = 'x_s05_1_8_32_0lay_hpds.txt';
   plotfile = 'x_s05_1_8_32_0lay_marg_a.png';
%   mapfile  = 'x_s04_1_3_16_0lay_fgsmap.txt';
%   meanfile = 'x_s04_1_3_16_0lay_mean.txt';
%   hpdfile  = 'x_s04_1_3_16_0lay_hpds.txt';
%   plotfile = 'x_s04_1_3_16_0lay_marg_a.png';
%   mapfile  = 'x_s02_1_3_25_6lay_fgsmap.txt';
%   meanfile = 'x_s02_1_3_25_6lay_mean.txt';
%   hpdfile  = 'x_s02_1_3_25_6lay_hpds.txt';
%   plotfile = 'x_s02_1_3_25_6lay_marg_a.png';
%   mapfile  = 'x_s01_1_3_20_5lay_fgsmap.txt';
%   meanfile = 'x_s01_1_3_20_5lay_mean.txt';
%   hpdfile  = 'x_s01_1_3_20_5lay_hpds.txt';
%   plotfile = 'x_s01_1_3_20_5lay_marg_a.png';
%   mapfile  = 'sim_T_1_4_17_0lay_fgsmap.txt';
%   meanfile = 'sim_T_1_4_17_0lay_mean.txt';
%   hpdfile  = 'sim_T_1_4_17_0lay_hpds.txt';
%   plotfile = 'sim_T_1_4_17_0lay_marg_a.png';
%   mapfile  = 'sim_C_1_7_40_6lay_fgsmap.txt';
%   meanfile = 'sim_C_1_7_40_6lay_mean.txt';
%   hpdfile  = 'sim_C_1_7_40_6lay_hpds.txt';
%   plotfile = 'sim_C_1_7_40_6lay_marg_a.png';
%   mapfile  = 'sim_B_1_10_50_1lay_fgsmap.txt';
%   meanfile = 'sim_B_1_10_50_1lay_mean.txt';
%   hpdfile  = 'sim_B_1_10_50_1lay_hpds.txt';
%   plotfile = 'sim_B_1_10_50_1lay_marg_a.png';
%   mapfile  = 'sim_A_1_3_16_2lay_fgsmap.txt';
%   meanfile = 'sim_A_1_3_16_2lay_mean.txt';
%   hpdfile  = 'sim_A_1_3_16_2lay_hpds.txt';
%   plotfile = 'sim_A_1_3_16_2lay_marg_a.png';
end

nsubfig = nx*ny;
xim = 0.03;
if(itrans == 1)
   yim = 0.3/ny;
else;
   yim = 0.1/ny;
end;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 30;
nbin2 = 10;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',10,'Height',1*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
%
% Sim B
%
%mtrue = [0.10, 1490, 1.45, 0.30, ...
%         0.30, 1550, 1.75, 0.25, ...
%         0.35, 1580, 1.95, 0.30, ...
%         3.50, 1540, 1.70, 0.10, ...
%               1600, 1.75, 0.50];

%load(ppd2);
%m2 = A(:,2:end);
%E2 = A(:,1);

load(data);
load(ppd);
if(isd == 0)
   m = A(:,2:end);
else
   m   = A(:,2:end-1);
   sd1 = A(:,end);
end;
if(itrans == 1)
   m = m(:,[1,4,2,7,5,3,6]);
end

E = A(:,1);

lo = F(1).minlim(1:end-3);
hi = F(1).maxlim(1:end-3);

lob = F(1).minlim(end-2:end);
hib = F(1).maxlim(end-2:end);

%% Site 21
%F(1).xticks = [1,2,3,4,5,6];
%F(2).xticks = [1500, 1600, 1700];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Site 16
%F(1).xticks = [1,2,3,4,5,6];
%F(2).xticks = [1500, 1600, 1700];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Site 13
%F(1).xticks = [1, 2,3];
%F(2).xticks = [1500, 1600, 1700];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Site 07
F(1).xticks = [0.50,1,1.5,2,2.5];
F(2).xticks = [1500, 1600, 1700, 1800];
F(3).xticks = [1.4, 1.6, 1.8, 2.0, 2.2];
F(4).xticks = [0, 0.5, 1.0 ];
F(5).xticks = [0.5, 1.0, 1.5, 2.0];

%% Site 05
%F(1).xticks = [0.25,0.50,0.75,1];
%F(2).xticks = [1500, 1600, 1700, 1800];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0, 2.2];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Site 04
%F(1).xticks = [0.5,1.00,1.5,2.0];
%F(3).xticks = [1450, 1500, 1550];
%F(2).xticks = [1.2, 1.4,1.6, 1.8];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Site 01
%F(1).xticks = [0.25,0.50,0.75,1];
%F(2).xticks = [1500, 1600, 1700, 1800];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0, 2.2];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Sim C:
%F(1).xticks = [0.5, 1, 1.5];
%F(2).xticks = [1500, 1600, 1700];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Sim B:
%F(1).xticks = [1,2,3];
%F(2).xticks = [1500, 1600, 1700];
%F(3).xticks = [1.4, 1.6, 1.8, 2.0];
%F(4).xticks = [0, 0.5, 1.0 ];

%% Sim A:
%F(1).xticks = [0.5,1,1.5];
%F(2).xticks = [1500, 1550, 1600];
%F(3).xticks = [1.4, 1.6, 1.8];
%F(4).xticks = [0, 0.5, 1.0 ];

%xparname = [{'Layer thickness (m)'},{'Density (g/ccm)'},{'Velocity (m/s)'},{'Atten. (dB/L)'},...
%            {'h'},{'Density (g/ccm)'},{'Velocity (m/s)'},{'nu'}];
%xparname = [{'Layer thickness (m)'},{'Velocity (m/s)'},{'Density (g/ccm)'},{'Atten. (dB/L)'},...
%            {'h'},{'Velocity (m/s)'},{'Density (g/ccm)'},{'nu'}];
xparname = [{'h1'},{'c1'},{'r1'},{'a1'};...
            {'h2'},{'c2'},{'r2'},{'a2'};...
            {'h3'},{'c3'},{'r3'},{'a3'};...
            {'h4'},{'c4'},{'r4'},{'a4'};...
            {'h5'},{'c5'},{'r5'},{'a5'};...
            {'h6'},{'c6'},{'r6'},{'a6'};...
            {'h7'},{'c7'},{'r7'},{'a7'};...
            {'h8'},{'c8'},{'r8'},{'a8'}];
xparname2 = [{'m1'},{'m2'},{'m3'},{'m4'}];

[cost_min i_min] = min(E);
xmap = m(i_min,:);
if(i_save == 1)
    save(mapfile,'xmap','-ASCII');
end

for ipar = 1:size(m,2)

    xmean(ipar) = mean(m(:,ipar));

end;
if(i_save == 1)
     save(meanfile,'xmean','-ASCII');
end;

[EE idxE] = min(E);
mstart = m(idxE,:);

npar = size(m,2);
nlay = (npar-3)/4;
%% Site 21:
%plo = [0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%            1448, 1.2, 0.0];
%phi = [2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%       2.0, 1710, 2.2, 1.0,...
%            1710, 2.2, 1.0];
%% Site 19:
%plo = [0.0, 1448, 1448, 1.2,...
%            1.2, 0.0, 0.0];
%phi = [2.5, 1900, 1900, 2.2,...
%            2.2, 1.0, 1.0];
%% Site 07:
plo = [0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0;...
       0.0, 1448, 1.2, 0.0];
phi = [2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0;...
       2.6, 1750, 2.2, 1.0];

if(isd == 1);
   plosd = [0.5];
   phisd = [2.];
end;

%% Site 04:
%plo = [0.0, 1.2, 1445, 0.0,...
%            1.2, 1445, 0.0];
%phi = [2.2, 1.8, 1555, 1.2,...
%            1.8, 1555, 1.2];
%% Site 16:
%plo = [0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%            1448, 1.2, 0.0];
%phi = [2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%       2.5, 1710, 2.2, 1.0,...
%            1710, 2.2, 1.0];
%% Site 13:
%plo = [0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%       0.0, 1448, 1.39, 0.0,...
%            1448, 1.39, 0.0];
%phi = [4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%       4.0, 1750, 2.1, 1.0,...
%            1750, 2.1, 1.0];
%
%% Site 05:
%plo = [0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%            1448, 1.2, 0.0];
%phi = [3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%       3.0, 1902, 2.2, 1.0,...
%            1902, 2.2, 1.0];
%% Site 01:
%plo = [0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%       0.0, 1448, 1.2, 0.0,...
%            1448, 1.2, 0.0];
%phi = [1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%       1.0, 1850, 2.4, 1.0,...
%            1850, 2.4, 1.0];
%% Sim C
%plo = [0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%       0.0, 1450, 1.2, 0.0,...
%            1450, 1.2, 0.0];
%phi = [2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%       2.0, 1750, 2.1, 1.0,...
%            1750, 2.1, 1.0];
%% Sim B
%plo = [0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%       0.0, 1449, 1.3, 0.0,...
%            1449, 1.3, 0.0];
%phi = [4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%       4.01, 1750, 2.2, 1.0,...
%             1750, 2.2, 1.0];
%% Sim A
%plo = [0.0, 1500, 1.3, 0.0,...
%       0.0, 1500, 1.3, 0.0,...
%       0.0, 1500, 1.3, 0.0,...
%            1500, 1.3, 0.0];
%phi = [2.00, 1650, 1.9, 1.0,...
%       2.00, 1650, 1.9, 1.0,...
%       2.00, 1650, 1.9, 1.0,...
%             1650, 1.9, 1.0];


nmods = size(m,1);
nfig = ceil(npar/nsubfig);

%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

for ipar = 1:npar

   [n1(:,ipar),lim(:,ipar)] = hist_tstar(m(:,ipar),E,nbin1,Tstar);
%   [n2(:,ipar),lim2(:,ipar)] = hist_tstar(m2(:,ipar),E2,nbin1,Tstar2);
   nf_int(ipar,:) = hpd(m(:,ipar),100,95);

end
if(isd == 1)
   [n1sd,lim1sd] = hist_tstar(sd1,E,nbin1,Tstar);
end;

%----------------------------------------------------------
%  Plot histograms:
%----------------------------------------------------------
%
% Make subplot mask:
%
idxsub = zeros(nlay+1,nx);
idxsub(1:end-1,[1,end-2:end]) = 1;
idxsub(end,end-2:end) = 1;

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc1=figure(1);
isubfig = 1;
ipar = 1;
for iy = 1:nlay+1
for ix = 1:nx

   if(idxsub(iy,ix) == 1)
      if(ipar > npar-3) 
         ylims_tmp(iy,ix+1) = max(n1(:,ipar));
      else
         ylims_tmp(iy,ix) = max(n1(:,ipar));
      end
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      hold on;box on;

      lim_2 = [min(m(:,ipar));min(m(:,ipar));
               lim(:,ipar); max(m(:,ipar)); max(m(:,ipar))];
      n1_2 = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
      fill(lim_2,n1_2,[0.8,0.8,0.8]);

%     lim2_2 = [min(m2(:,ipar));min(m2(:,ipar));
%               lim2(:,ipar); max(m2(:,ipar)); max(m2(:,ipar))];
%      n2_2 = [0; n2(1,ipar); n2(:,ipar); n2(end,ipar); 0];
%      stairs(lim2_2,n2_2);

      set(gca,'layer','top')
      ipar = ipar +1;
   end
   isubfig = isubfig + 1;
end
end

%ylims = max(max(ylims_tmp));
ylims = max(ylims_tmp);

ipar = 1;
isubfig = 1;
for iy = 1:nlay+1
for ix = 1:nx
   if(idxsub(iy,ix) == 1)
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      set(gca,'YTick',[],'FontSize',14);
      set(gca,'XTick',[F(ix).xticks],'FontSize',14);
      if(iy == nlay+1 | (iy == nlay & ix < nx-2))
         xlabel(xparname(iy,ix));
      else
         set(gca,'XTickLabel',[]);
      end 
      set(gca,'XLim',[plo(iy,ix) phi(iy,ix)],'FontSize',14);

      if(i_scaley == 1)
         if(ipar > npar-3) 
            set(gca,'YLim',[0 ylims(ix+1)+.05*ylims(ix+1)]);
         else
            set(gca,'YLim',[0 ylims(ix)+.05*ylims(ix)]);
         end 
      else
         if(iy > (npar-3)/4) 
            set(gca,'YLim',[0 ylims_tmp(iy,ix+1)+.05*ylims_tmp(iy,ix+1)]);
         else
            set(gca,'YLim',[0 ylims_tmp(iy,ix)+.05*ylims_tmp(iy,ix)]);
         end 
      end
      ipar = ipar +1;
   end
   isubfig = isubfig + 1;
end
end

if(isd == 1)
   gc2=figure(2);hold on;box on;
   set(gca,'XLim',[plosd phisd],'FontSize',14);
   size(sd1)
   size(lim1sd)
   lim1sd_2 = [min(sd1);min(sd1);
            lim1sd'; max(sd1); max(sd1)];
   n1sd_2 = [0; n1sd(1); n1sd'; n1sd(end); 0];
   fill(lim1sd_2,n1sd_2,[0.8,0.8,0.8]);

   set(gca,'layer','top')
end;

if(i_save == 1)
%   exportfig(gc1,plotfile,opts);
   saveas(gc1,plotfile,'png');
   saveas(gc2,plotsdfile,'png');
end

if(i_save == 1)
    save(hpdfile,'nf_int','-ASCII');
end

return;
