function [IMAP,ICOV,iar,i_varpar,irv,itr,iswd,iell,imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL,NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVsm,dVsp,dVpVsm, dVpVsp,sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile(parfile);  %for one-sided prior for dVs and dVpVs

fid = fopen(parfile);
%%  1
s = fgets(fid); %% IMAP:       Predict data for MAP file and exit.
IMAP = sscanf(s, '%d'); % convert to number
%%  2
s = fgets(fid); %% ICOV:       Selects likelihood function (0=implicit, 1=)
ICOV = sscanf(s, '%d'); % convert to number
%%  3
s = fgets(fid); %% IMAGSCALE:       Selects likelihood function (0=implicit, 1=)
IMAGSCALE = sscanf(s, '%d'); % convert to number
%%  4
s = fgets(fid); %% ENOS:       Selects likelihood function (0=implicit, 1=)
ENOS = sscanf(s, '%d'); % convert to number
%%  5
s = fgets(fid); %% IPOIPR:     Selects likelihood function (0=implicit, 1=)
IPOIPR = sscanf(s, '%d'); % convert to number
%%  6
s = fgets(fid); %% IAR:        Switches on autoregressive error model
iar = sscanf(s, '%d'); % convert to number
%%  7
s = fgets(fid); %% I_VARPAR:
i_varpar = sscanf(s, '%d'); % convert to number
%%  8
s = fgets(fid); %% IBD_SINGLE
ibd_single = sscanf(s, '%d'); % convert to number
%%  9
s = fgets(fid); %% I_RV:
irv = sscanf(s, '%d'); % convert to number
%%  10
s = fgets(fid); %% I_T: 
itr  = sscanf(s, '%d'); % convert to number
%%  11
s = fgets(fid); %% I_SWD:
iswd = sscanf(s, '%d'); % convert to number
%%  12
s = fgets(fid); %% I_ELL:
iell = sscanf(s, '%d'); % convert to number
%% 13
s = fgets(fid); %% I_MT:
imt = sscanf(s, '%d'); % convert to number
%% 14
s = fgets(fid); %% I_ZMT:
izmt = sscanf(s, '%d'); % convert to number
%%  15
s = fgets(fid); %% I_VREF:
ivref = sscanf(s, '%d'); % convert to number
%% 16
s = fgets(fid); %% I_VPVS:
ivpvs = sscanf(s, '%d'); % convert to number
%% 17
s = fgets(fid); %% ISMPPRIOR:  Switches on exchange moves (parallel tempering)
ISMPPRIOR = sscanf(s, '%d'); % convert to number
%% 18
s = fgets(fid); %% ISETSEED:  Switches on exchange moves (parallel tempering)
ISETSEED  = sscanf(s, '%d'); % convert to number
%% 19
s = fgets(fid); %% IEXCHANGE:  Switches on exchange moves (parallel tempering)
IEXCHANGE = sscanf(s, '%d'); % convert to number
%% 20
s = fgets(fid); %% IDIP:  Switches on exchange moves (parallel tempering)
idip = sscanf(s, '%d'); % convert to number
%% 21
s = fgets(fid); %% NDAT_SWD:    Max number of nodes.  
NDAT_SWD = sscanf(s, '%d'); % convert to number
%% 22
s = fgets(fid); %% NMODE:       Max number of nodes.  
NMODE = sscanf(s, '%d'); % convert to number
%% 23
s = fgets(fid); %% NDAT_ELL:    Max number of nodes.  
NDAT_ELL = sscanf(s, '%d'); % convert to number
%% 24
s = fgets(fid); %% NMODE_ELL:       Max number of nodes.  
NMODE_ELL = sscanf(s, '%d'); % convert to number
%% 25
s = fgets(fid); %% NDAT_MT:       Max number of nodes.  
NDAT_MT = sscanf(s, '%d'); % convert to number
%% 26
s = fgets(fid); %% NTIME:       Max number of nodes.  
NTIME = sscanf(s, '%d'); % convert to number
%% 27
s = fgets(fid); %% NSRC:       Max number of nodes.  
NSRC = sscanf(s, '%d'); % convert to number
%% 28
s = fgets(fid); %% NVMN:       Max number of nodes.  
NVMN = sscanf(s, '%d'); % convert to number
%% 29
s = fgets(fid); %% NVMX:       Max number of nodes.  
NVMX = sscanf(s, '%d'); % convert to number
%% 30
s = fgets(fid); %% ICHAINTHIN: Amount of chain thinning (keep low, more PT chains are better)
ICHAINTHIN = sscanf(s, '%d'); % convert to number
%% 31
s = fgets(fid); %% NKEEP:      Buffer length for saving sample.
NKEEP = sscanf(s, '%d'); % convert to number
%% 32
s = fgets(fid); %% NPTCHAINS1: No. chains at T = 1.
NPTCHAINS1 = sscanf(s, '%d'); % convert to number
%% 33
s = fgets(fid); %% dTlog:      PT chain spacing.
dTlog = sscanf(s, '%f'); % convert to number
%% 34
s = fgets(fid); %% lambda:      PT chain spacing.
lambda = sscanf(s, '%f'); % convert to number
%% 35
s = fgets(fid); %% hmx:      PT chain spacing.
hmx = sscanf(s, '%f'); % convert to number
%% 36
s = fgets(fid); %% hmin:      PT chain spacing.
hmin = sscanf(s, '%f'); % convert to number
%% 37
s = fgets(fid); %% armxH:      PT chain spacing.
armxH = sscanf(s, '%f'); % convert to number
%% 38
s = fgets(fid); %% armxV:      PT chain spacing.
armxV = sscanf(s, '%f'); % convert to number
%% 39
s = fgets(fid); %% armxSWD:      PT chain spacing.
armxSWD = sscanf(s, '%f'); % convert to number
%% 40
s = fgets(fid); %% armxSWD:      PT chain spacing.
armxELL = sscanf(s, '%f'); % convert to number
%% 41
s = fgets(fid); %% TCHCKPT:     Adapt the step size as function of acceptance.
TCHCKPT = sscanf(s, '%d'); % convert to number
%% 42
s = fgets(fid); %% shift:       Buffer length for adapting step size
shift = sscanf(s, '%f'); % convert to number
%% 43
s = fgets(fid); %% width:       Buffer length for adapting step size
width = sscanf(s, '%f'); % convert to number
%% 44
s = fgets(fid); %% wl:       Buffer length for adapting step size
wl = sscanf(s, '%f'); % convert to number
%% 45
s = fgets(fid); %% sampling_dt:       Buffer length for adapting step size
sampling_dt = sscanf(s, '%f'); % convert to number
%% 46
s = fgets(fid); %% dVsm:       Buffer length for adapting step size
dVsm = sscanf(s, '%f'); % convert to number
%% 47
s = fgets(fid); %% dVsp:       Buffer length for adapting step size
dVsp = sscanf(s, '%f'); % convert to number
%% 48
s = fgets(fid); %% dVpVsm:       Buffer length for adapting step size
dVpVsm = sscanf(s, '%f'); % convert to number
%% 49
s = fgets(fid); %% dVpVsp:       Buffer length for adapting step size
dVpVsp = sscanf(s, '%f'); % convert to number
%% 50
s = fgets(fid); %% sigmamin:       Buffer length for adapting step size
sigmamin = sscanf(s, '%f'); % convert to number
%% 51
s = fgets(fid); %% sigmamax:       Buffer length for adapting step size
sigmamax = sscanf(s, '%f'); % convert to number
%% 52
s = fgets(fid); %% sdmn:       Buffer length for adapting step size
sdmn = sscanf(s, '%f'); % convert to number
%% 53
s = fgets(fid); %% sdmx:       Buffer length for adapting step size
sdmx = sscanf(s, '%f'); % convert to number
%% 54
%s = fgets(fid); %% VpVsmin:       Buffer length for adapting step size
%VpVsmin = sscanf(s, '%f'); % convert to number
%% 55
%s = fgets(fid); %% VpVsmax:       Buffer length for adapting step size
%VpVsmax = sscanf(s, '%f'); % convert to number

%%
%% Read geometry file
%%
geomfile = 'sample.geom';
fid = fopen(geomfile);
itrace = 0;
while ~feof(fid);
  tline = fgets(fid);
  if(tline(1)=='#');
    continue;
  end;
  itrace = itrace + 1;
  tr(itrace,1:4) = sscanf(tline, '%f');

end;
baz = tr(:,1);
ntr = length(baz);

return;
