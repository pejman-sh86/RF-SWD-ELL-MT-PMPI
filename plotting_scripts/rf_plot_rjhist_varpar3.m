 !function [] = rf_plot_hist_varpar(filename);
 %%%%%%%%%%% for one sided prior for Vs and dVpVs
filename='simp04_sample.mat';
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

t1 = tic;

ifull   = 1; %% if 0, then use previously computed voronoi file
imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 2; %% Normalize profile marginals line by line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
isyn    = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ivoro   = 2; %% 1 -> Voronoi Cells; 2 -> Layer nodes
imap    = 0;
imead   = 0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imean   = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imax    = 0;
isave   = 1;
idatfit = 0;
idathist= 0;
iboudreau= 0;
isd     = 1;
add_to_vel_ref = 1; %%%%%%%%%%%%%%%%%%%%

NMISC = 0;
IEPS = 0;

filebase    = strrep(filename,'sample.mat','');
vorofile    = strrep(filename,'_sample.mat','_voro_sample.mat');
velreffile  = strrep(filename,'_sample.mat','_vel_ref.txt');
parfile     = strrep(filename,'_sample.mat','_parameter.dat')
datafile    = strrep(filename,'_sample.mat','_obs.dat');
datafileSWD = strrep(filename,'_sample.mat','_SWD.dat');
datafileELL = strrep(filename,'_sample.mat','_ELL.dat');
stffile     = strrep(filename,'_sample.mat','_stf_true.txt');
predfile    = strrep(filename,'_sample.mat','_pred.dat');
mapfiletrue = strrep(filename,'_sample.mat','_map_voro_true.dat');
profilefile = strcat(filebase,'profile.mat');
resfile     = strcat(filebase,'res.dat');
resarfile   = strcat(filebase,'resar.dat');
repensfile  = strcat(filebase,'rep_ens.dat');
%mapfile     = strrep(filename,'_sample.mat','_map_voro.dat'); %mine

[IMAP,ICOV,iar,i_varpar,irv,itr,iswd, iell, imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL, NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVsm,dVsp,dVpVsm,dVpVsp, sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile3(parfile);


NTIME2 = NSRC + NTIME - 1;
dt = sampling_dt;   
hmax   = hmx + 0.; %20.
hmax2 = hmax;
hmax2 = .1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zmin = 0.;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(irv == -1);irf = 1;else;irf = 0;end;
NRF    = ntr;
NSD    = 3*NRF + NMODE + NMODE_ELL +1;
IMAP
ICOV
iar
i_varpar
irv
itr
iswd
iell
imt
izmt
ivref
ivpvs
ISMPPRIOR
IEXCHANGE
idip
NDAT_SWD
NMODE
NDAT_ELL
NMODE_ELL
NDAT_MT
NTIME
NSRC
NVMN
NVMX
ICHAINTHIN
NKEEP
NPTCHAINS1
dTlog
hmx
hmin
armxH
armxV
TCHCKPT
shift
sampling_dt
dVsm
dVsp
dVpVsm
dVpVsp
sigmamin
sigmamax
sdmn
sdmx
ntr
baz

if (irv == -1 || irv == 1 || iswd == 1 || iell == 1)
    iseis = 1;
else
    iseis = 0;
end

if iseis == 1
    if ivpvs==1
    NPL= 3;
    else
        NPL = 2;
    end
    if(idip == 1);NPL = NPL+1;elseif(idip == 2);NPL = NPL+2;end;
else
    NPL = 1;
end
if (imt == 1); NPL = NPL+1; end;

if(idatfit == 1);
%if(irf == 0);
if irv >= 1
  predRfile = strcat(filebase,'predR.txt');
  predTfile = strcat(filebase,'predT.txt');
  predVfile = strcat(filebase,'predV.txt');
  predSfile = strcat(filebase,'predS.txt');
%else;
elseif irv == -1 %my change
  predRfile = strcat(filebase,'predRF.dat');
end;
if(iswd == 1);
  predSWDfile = strcat(filebase,'predSWD.dat');
  predSWDarfile = strcat(filebase,'predSWDar.dat');
end;
if(iell == 1);
  predELLfile = strcat(filebase,'predELL.dat');
  predELLarfile = strcat(filebase,'predELLar.dat');
end;
end;


plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';
plotfile1   = strcat(filebase,'varparhist.');
plotfile2   = strcat(filebase,'chainplots1.');
plotfile3   = strcat(filebase,'chainplots2.');
plotfile4   = strcat(filebase,'stddev.');
plotfile5   = strcat(filebase,'ARpar.');
plotfile6   = strcat(filebase,'kpar.');
plotfile7   = strcat(filebase,'numberpar.');
plotfile8   = strcat(filebase,'parcomplexity.');
plotfile9   = strcat(filebase,'logL.');
plotfile10  = strcat(filebase,'marg_prof.');
plotfile11  = strcat(filebase,'node_dens.');
plotfile12  = strcat(filebase,'data_fit.');
plotfile13  = strcat(filebase,'data_fit_swd.');
plotfile14  = strcat(filebase,'reshist.');
plotfile15  = strcat(filebase,'Axx.');
corefile    = 'core.mat';

%if(ivs == 1);
%  pmin = [ 3.0 1.60]';
%  pmax = [ 5.5 2.00]';
%  %pmin = [ 2.0 1.60]';
%  %pmax = [ 5.0 2.00]';
%else;
%  pmin = [ 3. 1.6  0.]';
%  pmax = [10. 2.0 30.]';
%end;
thinstep = 1;
if 1 == 2
if( isyn== 1 & irf == 0);
  stftru = dlmread(stffile);
  NSTFtru = length(stftru);
end;
end
if(irv==-1 || irv>=1) tmp = dlmread(datafile); end; %% mine if
if(iswd == 1);tmpSWD = dlmread(datafileSWD);end;
if(iell == 1);tmpELL = dlmread(datafileELL);end;
!if(irv >= 1);dobsR  = tmp(1:NRF,1:NTIME2);end;
if(irv==-1 || irv>=1);dobsR  = tmp(1:NRF,1:NTIME2);end; %% mine if 
if(itr == 1);dobsT  = tmp(2*NRF+1:3*NRF,1:NTIME2);end;
%tmp = dlmread(predfile);

if(idatfit == 1);
  %if(irf == 0);
  if irv >= 1 %mine
    dobsV   = tmp(NRF+1:2*NRF,1:NTIME2);
    predRraw   = load(predRfile);
    if(itr == 1);
    predTraw   = load(predTfile);end;
    predVraw   = load(predVfile);
    predSraw   = load(predSfile);
    NDsmp      = size(predRraw,1)/NRF;%% No data predictions
    predS = zeros(NDsmp,NRF,NSRC);
    predR = zeros(NDsmp,NRF,NTIME);
    predV = zeros(NDsmp,NRF,NTIME);
    itmp = 0;
    for ind=1:NDsmp;
      for iaz=1:NRF;
        itmp = itmp + 1;
        predS(ind,iaz,1:NSRC) = predSraw(itmp,1:NSRC);
        predR(ind,iaz,1:NTIME) = predRraw(itmp,1:NTIME);
        resR(ind,iaz,1:NTIME) = predRraw(itmp,1:NTIME)-dobsR(iaz,1:NTIME);
        if(itr == 1);
        predT(ind,iaz,1:NTIME) = predTraw(itmp,1:NTIME);
        resT(ind,iaz,1:NTIME) = predTraw(itmp,1:NTIME)-dobsT(iaz,1:NTIME);end;
        if(irf == 0);
        predV(ind,iaz,1:NTIME) = predVraw(itmp,1:NTIME);
        resV(ind,iaz,1:NTIME) = predVraw(itmp,1:NTIME)-dobsV(iaz,1:NTIME);end;
      end;
    end;
  %else
  elseif irv == -1 %% mine
    predRraw   = load(predRfile);
    NDsmp      = size(predRraw,1)/NRF;%% No data predictions
    itmp = 0;
    for ind=1:NDsmp;
      for iaz=1:NRF;
        itmp = itmp + 1;
        predR(ind,iaz,:) = predRraw(itmp,:);
        resR(ind,iaz,:) = predRraw(itmp,:)-dobsR(iaz,1:NTIME);
      end;
    end;
  end;
  if (irv == -1);NPRED   = size(predR,1); end; % mine if
%stop
if(iswd == 1);
  dobsSWD = tmpSWD(:,2);
  swdphase= tmpSWD(:,1);
  dobsSWD = dobsSWD';
  swdphase= swdphase';
  predSWD = load(predSWDfile);
  if(iar == 1);predSWDar = load(predSWDarfile);end;
  for i=1:size(predSWD,1);
    if(iar == 0);
      resSWD(i,:) = predSWD(i,:)-dobsSWD;
    else;
      resSWD(i,:) = predSWD(i,:)-dobsSWD-predSWDar(i,:);
    end;
  end;
end;
end;
%%
%% Ref velocity model:
%%
if(ivref == 1);
  tmp = dlmread(velreffile);
  NREF = tmp(1,1);
  ntmp = tmp(1,2);
  NPREM = ntmp-NREF;
  vel_ref = tmp(2:NREF+1,:);
%  vel_ref(NREF+1,:) = [hmax,tmp(end,2:end)];
  vel_prem = tmp(NREF+1:NREF+1+NPREM,:);
end;

%%
%% Source inversion simulation
%%
if(isyn == 1);
if(ivoro >= 1);
  %%
  %% True parameters for simulation
  tmpmap = load(mapfiletrue);
  ktru = tmpmap(1);
  for ivo = 1:ktru;
    vorotru(ivo,:) = tmpmap(2+(ivo-1)*NPL:2+(ivo*NPL)-1);
  end;
  voroidxtru = ones(size(vorotru));
  voroidxtru(vorotru<-99.) = 0;
  sdpartru(1:10) = [5.e-3,1.e-2,1.e-4,1.e-4,0.0171,0.0171,0.0171,0.0171,0.0171,0.0171]; %% true H value
  if irv == -1 || irv == 1 || iswd == 1 || iell == 1 %mine
  if(ivref == 1);
  for ivo = 1:ktru;
    if(voroidxtru(ivo,2) == 1);
      %vorotru(ivo,1),vel_ref
      [vref,vpvsref]=rf_getref(vorotru(ivo,1),vel_ref);
      vorotru(ivo,2) = vref + vorotru(ivo,2);
    end
    if ivpvs == 1 %mine
    if(voroidxtru(ivo,3) == 1);
      [vref,vpvsref]=rf_getref(vorotru(ivo,1),vel_ref);
      vorotru(ivo,3) = vpvsref + vorotru(ivo,3);
      %disp([vorotru(ivo,1),vref,vorotru(ivo,2)]);
    end
    end
  end;
  end; %mine
  end %mine
%   if (imt==1); vorotru(:,end) = 10.0 .^ (-vorotru(:,end));end;
if (imt==1); vorotru(:,end) = -vorotru(:,end);end;
  if(ivoro == 1);
    [mtru,nuniquetru] = rf_voro_to_lay(vorotru,voroidxtru,ktru,NPL,ktru);
  elseif(ivoro == 2);
    [mtru,nuniquetru] = rf_laynode_to_lay(vorotru,voroidxtru,ktru,NPL,ktru);
  end;
end;
end;

%% must replace with my prior
%% Prior bounds:
%%
% if 1 == 2
% pmin = [min(vel_ref(:,2))-dVs, min(vel_ref(:,3))-dVpVs]';
% pmax = [max(vel_ref(:,2))+dVs, max(vel_ref(:,3))+dVpVs]';
% if(idip == 1);
%   pmin(3) = [ 5];
%   pmax(3) = [45];
% elseif(idip == 2);
%   pmin(3:4) = [ 5,  5];
%   pmax(3:4) = [45, 50];
% end;
% end
% still needs correcting for dipping layres and NPL > 4 
% my comment: better to write it in terms of iseis and imt rather than NPL
% if NPL == 2 
%     if imt == 1
%         pmin = sigmamin;
%         pmax = sigmamax;
%     else %elseif ivpvs ~= 1
%        pmin = min(vel_ref(:,2))-dVs;
%        pmax = max(vel_ref(:,2))+dVs;
%    end
%elseif NPL == 3
%    if imt == 1
%    pmin = [min(vel_ref(:,2))-dVs, sigmamin]';
%    pmax = [max(vel_ref(:,2))+dVs; sigmamax]';
%    else %elseif ivpvs ==1
%      pmin = [min(vel_ref(:,2))-dVs, min(vel_ref(:,3))-dVpVs]';
%      pmax = [max(vel_ref(:,2))+dVs, max(vel_ref(:,3))+dVpVs]';
%    end
%elseif NPL == 4
%    pmin = [min(vel_ref(:,2))-dVs, min(vel_ref(:,3))-dVpVs]';
%    pmax = [max(vel_ref(:,2))+dVs, max(vel_ref(:,3))+dVpVs]';
%    pmin(3) = sigmamin;
%    pmax(3) = sigmamax;
%end

pmin = [];
pmax = [];
if iseis == 1
    pmin = [pmin; min(vel_ref(:,2))+dVsm];
    pmax = [pmax; max(vel_ref(:,2))+dVsp];
    if ivpvs == 1
       pmin = [pmin; min(vel_ref(:,3))+dVpVsm];
       pmax = [pmax; max(vel_ref(:,3))+dVpVsp];
    end
end
if imt == 1
%     pmin = [pmin; sigmamin];
%     pmax = [pmax; sigmamax];
%     pmin = [pmin; 10.0^(-sigmamax)];
%     pmax = [pmax; 10.0^(-sigmamin)];
pmin = [pmin; -sigmamax];
    pmax = [pmax; -sigmamin];
end

if(ifull == 1);
  %%
  %% Convert file from layer nodes to layers:
  %%
  t1_a = tic;
  disp('Converting voro to lay');
  if(ivoro == 1);
    rf_convert_sample_voro_to_lay(vorofile,NPL,NVMX,NMISC,NRF)
  elseif(ivoro == 2);
    rf_convert_sample_laynode_to_lay(vorofile,NPL,NVMX,NMISC,NRF,NMODE,NMODE_ELL,imt,iseis,ivref,ivpvs,hmax)
  end;
  t2 = toc(t1_a);
  disp('Done converting voro to lay. Time:');disp(t2);
end;
load(vorofile);
B=A;
load(filename);
if(ivoro == 0);B=A;end;
disp('Done loading sample.');

logLmin = min(A(:,1))-10;
logLmax = max(A(:,1))+10;
if(iar == 1)
   order = 1;
   armin = -0.6*ones(1,NRF);
   armax = 1.*ones(1,NRF);
else
   order = 1;
end;
NPROF = size(A,1);

k = A(:,4);
%% Compute prior volume
logP = zeros(size(k));
for i=1:size(k,1);
  vol = 1.;
  for ik=1:k(i);
    vol = vol * (hmax-(ik-1.)*hmin)*0.489*(pmax(end)-pmin(end));
  end;
  vol = vol * 0.489 * (pmax(end)-pmin(end));  %% This is half-space volume
  logP(i) = -log(vol);
end;

NFP = (k*NPL)+NPL-1;
m = A(:,5:end-5);

if(ifull == 1);
t1_a = tic;
%% Voronoi nodes sample file
if(ivoro >= 1);
  kv = B(:,4);
  gv = zeros(sum(kv),1);
  kvpar = zeros(size(kv,1),NPL);
  mv = B(:,5:end-5);
  jj1 = 1;jj2 = 1;jj3 = 1;jj4 = 1;jj5 = 1;
  %if imt == 1; isigma = NPL; jjs = 1; end; % my changes
  igv = 1;
  for i=1:size(B,1);
    for j=1:kv(i);
     %if irv == -1 || irv == 1 || iswd == 1   
      if(NPL > 1);
      if(mv(i,(j-1)*NPL+2) > -99.);
        Vs(jj1,:)   = mv(i,[(j-1)*NPL+1,(j-1)*NPL+2]);
        kvpar(i,1) = kvpar(i,1) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj1 = jj1 + 1;
      end;end;
      if(NPL > 2);
      if(mv(i,(j-1)*NPL+3) > -99.);
        VpVs(jj2,:)  = mv(i,[(j-1)*NPL+1,(j-1)*NPL+3]);
        kvpar(i,2) = kvpar(i,2) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj2 = jj2 + 1;
      end;end;
      if(NPL > 3);
      if(mv(i,(j-1)*NPL+4) > -99.);
        dip(jj3,:) = mv(i,[(j-1)*NPL+1,(j-1)*NPL+4]);
        kvpar(i,3) = kvpar(i,3) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj3 = jj3 + 1;
      end;end;
      if(NPL > 4);
      if(mv(i,(j-1)*NPL+5) > -99.);
        strike(jj4,:)   = mv(i,[(j-1)*NPL+1,(j-1)*NPL+5]);
        kvpar(i,4) = kvpar(i,4) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj4 = jj4 + 1;
      end;end;
      if(NPL > 5);
      if(mv(i,(j-1)*NPL+6) > -99.);
        Ks(jj5,:) = mv(i,[(j-1)*NPL+1,(j-1)*NPL+6]);
        kvpar(i,5) = kvpar(i,5) + 1;
        if(j>1);gv(igv) = gv(igv) + 1;end;
        jj5 = jj5 + 1;
      end;end;
     %end
%      if imt == 1
%         if(mv(i,(j-1)*NPL+isigma) > -99.);
%         sigma(jjs,:)   = mv(i,[(j-1)*NPL+1,(j-1)*NPL+isigma]);
%         kvpar(i,NPL) = kvpar(i,NPL) + 1;
%         if(j>1);gv(igv) = gv(igv) + 1;end;
%         jjs = jjs + 1;
%       end;
%      end
      %% Position parameter always present:
      kvpar(i,NPL) = kvpar(i,NPL) + 1;
      igv = igv + 1; 
    end;
  end;
  gv(find(gv==0))=[];
  nparv = sum(kvpar,2);
  t2 = toc(t1_a);
  disp('Done counting nodes. Time:');disp(t2);
  %%
  %% Save node file
  
 %%
  if(NPL == 2);
    save tmp_nodes.mat kvpar gv nparv Vs;
  elseif(NPL == 3);
    save tmp_nodes.mat kvpar gv nparv Vs VpVs;
  elseif(NPL == 4);
    save tmp_nodes.mat kvpar gv nparv Vs VpVs dip
  end;
  %save tmp_nodes.mat kvpar gv nparv Vs VpVs dip strike Ks;
end;
else
  %%
  %% Load node file
  %%
  kv = B(:,4);
  load tmp_nodes.mat;
end; %% end ifull
if(ivoro >= 1);
  %%
  %% Histograms of No. nodes, No. unique layers, No. parameters,
  %% No. parameters per node
  %%
  figw = 10;
  figh = 3;
  fig1 = figure();
  set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
  nx = 4;
  ny = 1;
  xim = 0.04;
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
  xlabel('No. nodes');
  set(gca,'XLim',[0 NVMX],'YLim',[0 0.4],'TickDir','out');
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
  %set(gca,'XLim',[min(A(:,4))-1 max(A(:,4))+1],'TickDir','out');
  set(gca,'XLim',[0 15],'YLim',[0 0.4],'TickDir','out');
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
  %set(gca,'XLim',[min(A(:,3))-1 max(A(:,3))+1],'TickDir','out');
  set(gca,'XLim',[10 40],'YLim',[0 0.4],'TickDir','out');
  if(isyn == 1);
    plot([sum(sum(voroidxtru)) sum(sum(voroidxtru))],[0 1],'-w');
    plot([sum(sum(voroidxtru)) sum(sum(voroidxtru))],[0 1],'--k');
  end;

  h4 = subplot('Position',[loc(1,4) loc(2,4) spw sph]);
  hold on; box on;set(gca,'FontSize',14);
  [n1,lim]=hist(gv,[1:NPL]);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
  n1 = n1/sum(n1);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n1,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n1,'k');
  clear n1 lim;
  xlabel('No. parameters per node');
  set(gca,'XLim',[0 NPL],'YLim',[0 0.7],'TickDir','out');
  set(gca,'XTick',[0:1:4]);

end;
%%
%% MAP for max(P(k))
%%
kmin = min(k);
kmax = max(k);
if(kmax > 0);
  for i = 1:kmax-kmin+1;
    ik = kmin+i-1;
    idx = find(A(:,4) == ik);
    nmod_k(i) = size(idx,1);
    [a,b] = max(A(idx,1));
    mapk(i).par=A(idx(b),5:4+A(idx(b),4)*NPL+(NPL-1));
    clear idx;
  end
  [a,b] = max(nmod_k);
  kpeak = b+kmin;
  idx = find(A(:,4) == kpeak);
  [a,b] = max(A(idx,1));map=A(idx(b),4:end-6);
  clear idx;
else
  [a,b] = max(A(:,1));map=A(b,4:end-6);
end;
%save(mapfile,'map','-ascii');

if(iar == 1)
   (size(B,2)-6-order*NRF+1:size(B,2)-6)
   alpha = B(:,end-6-(order*NRF)+1:end-6);
end

logL = A(:,1);
if(imap == 1);
   [logLmap,jmap] = max(logL);
   kmap = k(jmap);
   NFPmap = NFP(jmap);
   mmap = m(jmap,1:NFPmap);
end;

if(kmax > 0);
  h   = zeros(size(A,1),max(A(:,4)));
  dep = zeros(size(A,1),max(A(:,4)));
  for i = 1:size(A,1);
    if(A(:,4)>0);
      idxh = [0:k(i)-1]*NPL+1;
      %h(i,1:A(i,4)) = m(i,idxh(end));
        dep(i,1:k(i)) = m(i,idxh);
        h(i,1:k(i)) = diff([0,dep(i,1:k(i))]);
      clear idxh;
    end;
  end;
  dep = dep(:);dep(find(dep==0))=[];
  h = h(:);h(find(h==0))=[];
  [a_enos,b_enos]=hist(h,200);a_enos=a_enos/trapz(b_enos,a_enos);
  a_enos=[0,a_enos,0];
  b_enos=[b_enos(1),b_enos,b_enos(end)];
  fig2 = figure('visible','on');
  [xx,yy]=stairs(b_enos,a_enos,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
end;
if(isd == 1)
  %% The order of std devs in the sample is 
  %% sd_H, sd_V, sd_T, sd_SWD, sd_MT
  sdlabel = [{'Radial'},{'Vertical'},{'Transverse'},{'SWD'}, {'ELL'}, {'MT'}];
  sd(:,1:NSD) = A(:,end-5-(3*NRF)-NMODE-NMODE_ELL-NSD:end-6-(3*NRF)-NMODE-NMODE_ELL);
  disp('Done getting sigma. Indeces:');
  disp([size(A,2)-5-(3*NRF)-NMODE-NMODE_ELL-NSD:size(A,2)-6-(3*NRF)-NMODE-NMODE_ELL]);
  
  if(itr == 0);
    sd(:,3) = [];
    NSD = NSD - 1;
    sdlabel(3) = [];
  end;
  if(irv == -1); % my change
    sd(:,2) = [];
    NSD = NSD - 1;
    sdlabel(2) = [];
  end;
  if(irv == 0); % my change
    sd(:,1:2) = [];
    NSD = NSD - 2;
    sdlabel(1:2) = [];
  end;
  if(iswd == 0);
    sd(:,end-2) = [];
    NSD = NSD - 1;
    sdlabel(end-2) = [];
  end;
  if(iell == 0);
    sd(:,end-1) = [];
    NSD = NSD - 1;
    sdlabel(end-1) = [];
  end;
  if(imt == 0); %my change
    sd(:,end) = [];
    NSD = NSD - 1;
    sdlabel(end) = [];
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% CHAIN PLOTS
%%
fig2=figure;
NPT = max(A(:,end-1));
nx = NPT;
ny = nx;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ith = 1;
ipt = 1;
jj=0;
for j=1:NPT;
   jj=jj+1;
   h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
   hold on; box off;
   set(gca,'FontSize',14);

   clear idx1 idx2;
   if(find(A(:,end)==i));
   idx1 = find(A(:,end)==i);
   idx2 = find(A(idx1,end-1)==j);
   if(ivoro == 0);
     [AX,H1,H2] = plotyy([1:length(idx2)],...
                A(idx1(idx2),1),[1:length(idx2)],...
                A(idx1(idx2),4));
   else
     [AX,H1,H2] = plotyy([1:length(idx2)],...
                B(idx1(idx2),1),[1:length(idx2)],...
                B(idx1(idx2),4));
   end;
   set(AX(1),'YTick',[0:40:1000]);
   if(rem(i-1,nx)==0);
      set(get(AX(1),'Ylabel'),'String','logL') 
      set(AX(1),'YTickLabel',[0:40:1000]);
   else;
      set(AX(1),'YTickLabel',[]);
   end;
   set(AX(2),'YTick',[0:2:100]);
   if(rem(i,nx)==0);
      set(get(AX(2),'Ylabel'),'String','No. interfaces') 
      set(AX(2),'YTickLabel',[0:2:100]);
   else;
      set(AX(2),'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   set(AX(2),'XTickLabel',[],'YLim',[kmin-1 kmax+1]);
   set(gca,'YLim',[logLmin logLmax]);
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Posterior CHAIN PLOTS
%%
fig3 = figure;
NPT = max(A(:,end-1));
nx = NPT;
ny = nx;
xim = 0.01;
yim = 0.06;
xymarg = [0.1 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

ith = 1;
ipt = 1;
jj=0;
for j=1:NPT;
   jj=jj+1;
   h1 = subplot('Position',[loc(1,jj) loc(2,jj) spw sph]);
   hold on; box off;
   set(gca,'FontSize',14);

   clear idx1 idx2;
   if(find(A(:,end)==i));
   idx1 = find(A(:,end)==i);
   idx2 = find(A(idx1,end-1)==j);
   plot([1:length(idx2)],A(idx1(idx2),1)+logP(idx1(idx2)));
   set(gca,'YTick',[0:40:1000]);
   if(rem(i-1,nx)==0);
      set(gca,'YTickLabel',[0:40:1000]);
   else;
      set(gca,'YTickLabel',[]);
   end;
   if(i>NTH-((nx*ny)-NTH+1));
      xlabel('rjMCMC step');
   end;
   %set(gca,'YLim',[logLmin+max(logP) logLmax+min(logP)]);
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SIGMA PLOT
%%
%isyn = 1;
Npanel = 5;
NSD2 = 3;
font_size = 25;
% sdlabel = {'RF', 'SWD', 'H/V', 'MT'};
sdlabel = {'MT', 'HV', 'SWD'};
% sdlabel = {'SWD'};
% sdlabel = {'SWD', 'RF','MT'};
% nbins = [50 75 40];  %real
nbins = [50 50 50];
% sdpartru(1:10) = [1.e-4,1.e-2,1.e-2,1.e-4,0.0171,0.0171,0.0171,0.0171,0.0171,0.0171];
sdpartru = [2.e-2,5.e-2,3.e-2];
% sdpartru = 3.e-2;
ilabel = 1;
if(isd == 1)
   sdmax = [4.1e-2, .9e-1, 4.1e-2];%[1.e-3, 1.e-1, 1.e-1];%[4., 2., 1.2];%[0.03, 0.12, 0.03];%
   sdmin = [1.e-2, 2.e-2, .5e-2];%[1.e-5, 1.e-3, 1.e-3];%[1.e-0, 1.e-0, 1.e-0];%
   Pmax = .1; %0.05 for sim 0.23 real
   fig4=figure;
   %nx = NSD/3;
   %ny = 3;
   nx = Npanel;
   ny = 3;
   xim = 0.01;
   yim = 0.08;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
%    spw = spw - 0.005;
%    spw2 = spw;
%    spw = spw * 3./4.;
%    loc(1,5) = loc(1,5) - spw2/4.;
%    loc(1,6) = loc(1,6) - 2.*spw2/4.;
%    XTick = { [7.e-5, 1.e-4, 1.5e-4], [7.e-3, 1.e-2, 1.5e-2], [7.e-3, 1.e-2, 1.5e-2] }; %sim
%    exponent = [-4, -2, -2]; %sim
XTick = { [2. 3.]*1.e-2, [ 3. 5.]*1.e-2, [1. 3.]*1e-2 }; %sim
% XTick = { [1., 2., 3.], [1., 1.5], [1., 1.1] }; %real
   ilabel = 1;
   sd2 = zeros(size(sd,1), NSD2);

%    make the order based on the present dataset
   sd2(:,1:NSD) = sd(:,end:-1:1); %%check sdpartru, sdmax, andsdmin to be in reverse order, too

%    make the order fixed: always MT-SWD-RF
%     if(imt==1)
%         sd2(:,1:NSD) = sd(:,end:-1:1);
%     elseif(iswd==1)
%         sd2(:,2:NSD+1) = sd(:,end:-1:1);
%     elseif(irv==-1)
%         sd2(:,3:NSD+2) = sd(:,end:-1:1);
%     end

   for i=1:Npanel;
     h1 = subplot('Position',[loc(1,i+Npanel) loc(2,i+Npanel) spw sph]);
     hold on; box on;
     set(gca,'FontSize',font_size);
     if(i > 1);set(gca,'YTickLabel',[]);end;
     if(i>NSD); set(gca,'XTickLabel',[]);end

     if i<=NSD
     [n1,lim]=hist(sd2(:,i),nbins(i));n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
     n1 = n1/sum(n1);
     lim = lim - (lim(3)-lim(2))/2;
%      [xx,yy]=stairs(n1,lim,'k');
     [xx,yy]=stairs(lim,n1,'k');
     patch(xx,yy,[0.8,0.8,0.8]);
     stairs(lim,n1,'k');
     clear n1 lim;
%      if(i == 1);ylabel('Data error standard deviation');end;
%      xlabel('Probability');
     if(i == 1);ylabel('Probability');end;
     %xlabel('std');
%      %if(i==1);isdmx = 1;end;
%      %if(i==NRF+1);isdmx = 2;end;
%      %if(i==2*NRF+1);isdmx = 3;end;
     isdmx = i; %%mine
%      %set(gca,'XLim',[0.0 sdmax(isdmx)],'YTick',[0.0:0.01:0.1])
%      %set(gca,'YLim',[0.0 0.05],'XTick',[0.0:0.01:0.1])
     set(gca,'XLim',[sdmin(isdmx), sdmax(isdmx)]) %,'YTick',[0.0:0.01:0.1])
     set(gca,'YLim',[0.0 Pmax])
     set(gca,'XTick',XTick{i});
%      ax = gca;
%      ax.XAxis.Exponent = exponent(i);
%      if(i > 1);set(gca,'YTickLabel',[]);end;
     box on;
%      %text(.01,.0475,sdlabel(ilabel),'FontSize',12,'Color',[0,0,0]);  %mine
     if(isyn == 1);
       plot([sdpartru(i) sdpartru(i)],[0 Pmax],'-w');
       plot([sdpartru(i) sdpartru(i)],[0 Pmax],'--k');
     end;
     %ilabel = floor(i/NRF)+1; %%mine
     label = strcat( '{',sdlabel(ilabel) );
     label = strcat( label, '}' );
     xlabel( strcat('s_', label) )
     ilabel = ilabel+1;
   end;
   end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% ALPHA PLOT
%%
if(iar == 1)
   figw = 12;
   figh = 8;
   fig5=figure('visible','on');
   set(fig5,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = NRF;
   ny = 3;
   xim = 0.01;
   yim = 0.01;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   for i=1:NRF
     h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);
     h2 = subplot('Position',[loc(1,i+NRF) loc(2,i+NRF) spw sph]);
     hold on; box off;
     set(gca,'FontSize',12);

     subplot(h1);hold on;box on;
     [n,lim]=hist(alpha(:,i),60);n = [0, n, 0];lim = [lim(1) lim lim(end)];
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

     idxar = ones(size(alpha(:,i)));
     idx = find(alpha(:,i) < -0.5);
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
fig6=figure;
subplot(2,1,1);hold on;box on;
set(gca,'FontSize',14);
if(ivoro == 0);
  plot([1:thinstep:length(k)],k(1:thinstep:end),'k')
else;
  plot([1:thinstep:length(kv)],kv(1:thinstep:end),'k')
end;
ylabel('No. nodes');
xlabel('rjMCMC step');
set(gca,'XLim',[0 length(k)],'TickDir','out');
set(gca,'YLim',[min(kv)-1 max(kv)+1],'TickDir','out');

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
%plot([1:20], poisson_dist([1:20], 10))
if isyn == 1
plot([ktru,ktru], [0,1],'--k') %mine
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% No. parameters PLOT
%%
if(ivoro >= 1);
fig7=figure;
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
fig8=figure;
font_size = 25;
font_size2 = 25;
figw = 12;
figh = 6;
set(fig8,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
NPL2 = 6;%4;  %%number of panels + 1
nx = NPL2-1;%%NPL-1;
ny = 3;
xim = 0.01;
yim = 0.08;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
%sph = spw;
% spw = spw - 0.005;
% spw2 = spw;
% spw = spw * 3./4.;
% loc(1,5) = loc(1,5) - spw2/4.;
% loc(1,6) = loc(1,6) - 2.*spw2/4.;
XTick = [5:5:10];  
Pmax = 1.; %real
% Pmax =1;
NVMX2 = 10;

%make resistivity first if it is present
if (imt==1)
    kvpar2 = zeros(size(kvpar));
    kvpar2(:,1) = kvpar(:,end-1);
    kvpar2(:,2:end-1) = kvpar(:,1:end-2);
    kvpar2(:,end) = kvpar(:,end);
else
    kvpar2 = kvpar;
end

% make resistivity always first
% kvpar2 = zeros(size(kvpar,1),NPL2);
% if (imt==1)
%     kvpar2(:,1)=kvpar(:,end-1);
%     kvpar2(:,2:NPL-1) = kvpar(:,1:end-2);
%     kvpar2(:,NPL2) = kvpar(:,end);
% else
%     kvpar2(:,2:NPL+1) = kvpar(:,1:end);
% end

% h1=[];
h2=[];

for ipar=1:NPL2-1;
%   hh1 = subplot('Position',[loc(1,ipar) loc(2,ipar) spw sph]); % I have subtracted 0.01
%   hold on; box on;
%   set(gca,'FontSize',20);
%   set(gca,'YLim',[0 NVMX],'TickDir','out');
%   set(gca,'XLim',[0 length(k)],'TickDir','out');
%   if(ipar==1);ylabel('No. nodes');
%   else;set(gca,'YTickLabel',[]);end;
% %   xlabel('rjMCMC step');
%   if ipar<=1
  hh2 = subplot('Position',[loc(1,ipar+NPL2-1) loc(2,ipar+NPL2-1) spw sph]);
  set(gca,'XLim',[0 NVMX2],'TickDir','out');
  set(gca,'XTick',XTick)
  set(gca,'YLim',[0. Pmax],'TickDir','in');
  box on;
  hold on; %box off;
  set(gca,'FontSize',font_size);
  if ipar>1; set(gca,'YTickLabel',[]); end;
%   h1 = [h1, hh1];
  h2 = [h2, hh2];
%   end
end

  for ipar=1:NPL-1%NPL2-1  %NPL2: resistivity always first NPL:resistivity first if present

%   subplot(h1(ipar));
%   plot([1:thinstep:length(kvpar2(:,ipar))],kvpar2(1:thinstep:end,ipar),'k')

  subplot(h2(ipar));
  [n,lim]=hist(kvpar2(:,ipar),[0:NVMX]);
  n = [0, n, 0];
  lim = [lim(1) lim lim(end)];
  n = n/sum(n);
  lim = lim - (lim(3)-lim(2))/2;
  [xx,yy]=stairs(lim,n,'k');
  patch(xx,yy,[0.8,0.8,0.8]);
  stairs(lim,n,'k');
  clear n lim;
  %if(ipar == 1);xlabel('No. nodes with V_P');end;
  %if(ipar == 1);xlabel('No. nodes with V_S');end;
  if(ipar==1)
      if(imt==1)
          xlabel('No. nodes resistivity','FontSize',font_size2);
%           xlabel('k_{\rho}');
      else
          xlabel('No. nodes with V_P');
          xlabel('No. nodes V_S','FontSize',font_size2);
%           xlabel('k_{V_S}');
      end
  end
  if(ipar == 2)
      if (imt==1)
          xlabel('No. nodes V_S','FontSize',font_size2);
%            xlabel('k_{V_S}');
      else
          if (ivpvs==1)
              xlabel('No. nodes V_P/V_S','FontSize',font_size2);
%               xlabel('k_{V_P/V_S}');
          end
      end
% %       if (ivpvs==1)
% %           xlabel('No. nodes V_S','FontSize',font_size2);
% %       else
% %           xlabel('No. nodes log-conductivity','FontSize',font_size2);
% %       end
  end
  if(ipar == 3);xlabel('No. nodes with dip');end;
%   if(ipar == 3);xlabel('No. nodes log-conductivity','FontSize',font_size2);end;
  if(ipar == 3);xlabel('No. nodes V_P/V_S','FontSize',font_size2);end;
%   if(ipar == 3);xlabel('k_{V_P/V_S}');;end;
  if(ipar == 4);xlabel('No. nodes with strike');end;
  if(ipar == 5);xlabel('No. nodes with \gamma_S');end;
  if(ipar == 1);ylabel('Probability');
  else;set(gca,'YTickLabel',[]);end;
  if(isyn == 1);

%       make resistivity first if present
      if (imt==1)
          voroidxtru2 = zeros(size(voroidxtru));
          voroidxtru2(:,2) = voroidxtru(:,end);
          voroidxtru2(:,3:end) = voroidxtru(:,2:end-1);
          voroidxtru2(:,1) = voroidxtru(:,1); 
      else
          voroidxtru2 = voroidxtru; 
      end

%     make resistivity always first
%       voroidxtru2 = zeros(size(voroidxtru,1),NPL2);
%       if (imt==1)
%           voroidxtru2(:,2) = voroidxtru(:,end);
%           voroidxtru2(:,3:NPL) = voroidxtru(:,2:end-1);
%       else
%           voroidxtru2(:,3:NPL+1) = voroidxtru(:,2:end);
%       end
%       voroidxtru2(:,1) = voroidxtru(:,1);

    plot([sum(voroidxtru2(:,ipar+1)) sum(voroidxtru2(:,ipar+1))],[0 1],'-w');  % I have changed voroidxtru to voroidxtru2
    plot([sum(voroidxtru2(:,ipar+1)) sum(voroidxtru2(:,ipar+1))],[0 1],'--k');
%     plot([3 3],[0 1],'-w');  % I have changed voroidxtru to voroidxtru2
%     plot([3 3],[0 1],'--k');
%       plot([sum(voroidxtru2(:,ipar+1)) sum(voroidxtru2(:,ipar+1))],[0 1],'--r');
  end;
%   set(gca,'XLim',[0 NVMX2],'TickDir','out');
%   set(gca,'YLim',[ 0. 1.0],'TickDir','in');
%   box on;
 end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% logL PLOT
%%
fig9=figure;
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',14);

%for i=1:max(A(:,end));
%   idx = find(A(:,end)==i);
%   plot(A(idx(1:thinstep:end),1),'k');
%   clear idx;
   plot([1:thinstep:length(logL)],logL(1:thinstep:end),'k');
%end;

ylabel('log Likelihood');
xlabel('rjMCMC step');
%set(gca,'XLim',[0 length(logL)])
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
   NZ = 1000; %400; %4000, 2000
   NZI= 1000; %400; %4000, 2000
   nsmooth= ceil(NZ/80.);
   NC = 200; %200; %1000, 400
   NR = 200;
   NA = 200;
   zlim = cumsum(hmx/NZI*ones(1,NZI));
   clim = pmin(1)+cumsum((pmax(1)-pmin(1))/NC*ones(1,NC));
   if NPL > 2
   rlim = pmin(2)+cumsum((pmax(2)-pmin(2))/NR*ones(1,NR));
   end
   if(NPL>3);
   alim = pmin(3)+cumsum((pmax(3)-pmin(3))/NA*ones(1,NA));
   end;

   dz  = hmax/(NZ-1);
   z   = cumsum(dz*ones(1,NZ))-dz;
   dzi = hmax/(NZI-1);
   zi  = cumsum(dzi*ones(1,NZI))-dzi;

   h = zeros(1,NZI);
   nlo = zeros(NZ,1);
   nhi = zeros(NZ,1);
   c = zeros(NPROF,NZ);
   if NPL>2; r = zeros(NPROF,NZ); end;
   %s = zeros(NPROF,NZ);
   if NPL > 3;a = zeros(NPROF,NZ);end;
   disp('Sample size: '),disp(size(m))
   for iprof = 1:NPROF

     if(rem(iprof,5000)==0)
        fprintf(1,'%8i',iprof)
     end
     clear idxh idxc idxr idxa prof;
     %% Find index for current model
     if(k(iprof) > 0)
        idxh = (([1:k(iprof)]-1)*NPL)+1;
        idxc = (([1:k(iprof)]-1)*NPL)+2;
        if NPL>2;idxr = (([1:k(iprof)]-1)*NPL)+3;end;
        if NPL>3;idxa = (([1:k(iprof)]-1)*NPL)+4;end;
        idxh = [idxh idxh(end)];
        idxc = [idxc idxc(end)+NPL-1];
        if NPL>2; idxr = [idxr idxr(end)+NPL-1];end;
        if NPL > 3; idxa = [idxa idxa(end)+NPL-1];end;
     else
        idxh = [];
        idxc = [1];
        if NPL>2; idxr = [2]; end;
        if NPL>3; idxa = [3]; end;
     end

     %% Compute the profile for current model
     if(k(iprof) > 0)
       prof(1:k(iprof),1) = cumsum(m(iprof,idxh(1:end-1)),2);
       prof(1:k(iprof),1) = m(iprof,idxh(1:end-1));
       prof(k(iprof)+1,1) = prof(k(iprof),1)+m(iprof,idxh(end));
       prof(:,2) = m(iprof,idxc);
       if NPL>2; prof(:,3) = m(iprof,idxr); end;
       if NPL>3 prof(:,4) = m(iprof,idxa);end;

       for ilay=2:k(iprof)+1  %% k is # layers of current model
          idxzi = round(prof(ilay-1,1)/dzi);
          if(idxzi == 0)idxzi = 1;end;
          h(idxzi) = h(idxzi) + 1;
       end;
       c(iprof,:) = prof(1,2);
       if NPL>2; r(iprof,:) = prof(1,3); end;
       if NPL>3; a(iprof,:) = prof(1,4); end;
       for ilay=2:k(iprof)+1  %% k is # layers of current model
          idxz = round(prof(ilay-1,1)/dz);
          if(idxz == 0)idxz = 1;end;
%          h(idxz)     = h(idxz) + 1;
          c(iprof,idxz:end) = prof(ilay,2);
          if NPL>2; r(iprof,idxz:end) = prof(ilay,3); end;
          if NPL>3; a(iprof,idxz:end) = prof(ilay,4); end;
       end;

     else
       prof(:,2) = m(iprof,idxc);
       if NPL>2; prof(:,3) = m(iprof,idxr); end;
       if NPL>3; prof(:,4) = m(iprof,idxa); end;
       c(iprof,:) = prof(1,2);
       if NPL>2; r(iprof,:) = prof(1,3); end;
       if NPL>3; a(iprof,:) = prof(1,4); end;
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
      [nfc(iz,:)] = hpd(c(:,iz),100,95);
      meac(iz) = median(c(:,iz));
      if NPL>2
        [Nr(iz,:),binsr] = hist(r(:,iz),rlim);
        [nfr(iz,:)] = hpd(r(:,iz),100,95);
        mear(iz) = median(r(:,iz));
      end
      if(NPL>3);
        [Ndip(iz,:),binsa] = hist(a(:,iz),alim);
        [nfa(iz,:)] = hpd(a(:,iz),100,95);
        meaa(iz) = median(a(:,iz));
      end;
   end;
   [Ndep,binsdep] = hist(dep,zlim);
   %
   % Normalize Histograms
   %
   if(inorm == 0)
     for iz=1:NZ
       Nc(iz,:) = Nc(iz,:)/NPROF;
       if NPL>2; Nr(iz,:) = Nr(iz,:)/NPROF; end;
       if(NPL>3);Ndip(iz,:) = Ndip(iz,:)/NPROF;end;
     end;
   elseif(inorm == 1)
      for iz=1:NZ
         Nc(iz,:) = Nc(iz,:)/trapz(binsc,Nc(iz,:));
         if NPL>2; Nr(iz,:) = Nr(iz,:)/trapz(binsr,Nr(iz,:)); end;
         if(NPL>3);Ndip(iz,:) = Ndip(iz,:)/trapz(binsa,Ndip(iz,:));end;
      end;
       Ndep = Ndep/trapz(binsdep,Ndep);
       Ndep = [0,Ndep,0];binsdep = [binsdep(1),binsdep,binsdep(end)];
   elseif(inorm == 2)
     for iz=1:NZ
       Nc(iz,:) = Nc(iz,:)/max(Nc(iz,:));
       if NPL>2; Nr(iz,:) = Nr(iz,:)/max(Nr(iz,:)); end;
       if(NPL>3);Ndip(iz,:) = Ndip(iz,:)/max(Ndip(iz,:));end;
     end;
      Ndep = Ndep/trapz(binsdep,Ndep);
      Ndep = [0,Ndep,0];binsdep = [binsdep(1),binsdep,binsdep(end)];
   end;
   disp('Done histograms.');

c_mean = mean(c);
c_mead = median(c);
if NPL>2
r_mean = mean(r);
r_mead = median(r);
end
if(NPL>3);
  a_mean = mean(a);
  a_mead = median(a);
  %[ntmp,idxamax] = max(Ndip,[],2); %my change %
end;

[ntmp,idxcmax] = max(Nc,[],2);
if NPL>2 [ntmp,idxrmax] = max(Nr,[],2); end;
if NPL>3 [ntmp,idxamax] = max(Ndip,[],2);end; %mine
for iz=1:NZ
   c_max(iz) = clim(idxcmax(iz));
   if NPL>2; r_max(iz) = rlim(idxrmax(iz)); end;
   if(NPL>3);a_max(iz) = alim(idxamax(iz));end;
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
   if NPL>2; r_max_sm(i) = sqrt(mean(r_max(NAVE1:NAVE2).^2)); end;
   if(NPL>3);a_max_sm(i) = sqrt(mean(a_max(NAVE1:NAVE2).^2));end;
end

end; % end imarg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% PLOT PROFILE MARGINALS
%%
fig10 = figure;hold on; box on;
% set(fig10, 'renderer', 'opengl')
cmap = colormap( bone(1000) );
cmap = flip(cmap,1);
cmap(1:1,:) = ones(1,3);%[1 1 1];
colormap(cmap);

% imean=0;
plot_interfaces = 1;

font_size = 20;

lwt = .5;
lwp = 2;
plot_prior_Vs = 0;
plot_prior_VpVs = 0;
plot_prior_res = 0;

xxlim = [0. 85.];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.55 for sim 0.11 and .8 (.175) for real
% pmin2 = pmin;
% pmax2 = pmax;
xlims_min = [pmin(1) pmin(1) 1.];  % sim %%[Vs VpVs log_/rho]
xlims_max = [1.3 pmax(1) pmax(end)];
% xlims_min = [3.6 pmin(1) 2.6]; %real %%[Vs VpVs log_/rho]
% xlims_max = [4.8  pmax(1) 3.6];
if (iseis==1 && ivpvs==1); xlims_min(2) = pmin(2); xlims_max(2) = pmax(2);end
pmin2 = [];
pmax2 = [];
if iseis == 1
    pmin2 = [pmin2; xlims_min(1)];
    pmax2 = [pmax2; xlims_max(1)];
    if ivpvs == 1
       pmin2 = [pmin2; xlims_min(2)];
       pmax2 = [pmax2; xlims_max(2)];
    end
end
if imt == 1
pmin2 = [pmin2; xlims_min(end)];
pmax2 = [pmax2; xlims_max(end)];
end

XTick_res = [1:2.:4.];
XTick_Vs = [1:1.:6.];
XTick_Vs = [.2 1. 2.];
XTick_VpVs = [1.8:.8:3.];
%isyn = 0;
opts = struct('bounds','tight','LockAxes',1, ...
              'Width',8,'Height',4.8,'Color','cmyk',...
              'Renderer','painters','Format','png',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
%[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

% loc = [0.08,  0.18, 0.4533, 0.5633;...
%        0.14, 0.14, 0.14, 0.14];
%spw1 = 0.09;
spw1 = 0.05;
%spw2 = 0.2633;
spw2 = 0.1;
%spw3 = 0.10;
spw3 = 0.05;
%sph = 0.82;
sph = .51;
sep = 0.01;
% % loc = [0.08,  0.01+0.18, 0.03+0.4530-0.02, 0.05+0.7266-0.04;...
% %        0.14, 0.14, 0.14, 0.14];
% spw = [spw1, spw2, spw2, spw3];
% loc = zeros(2,4);
% loc(2,:) = .14;
% %loc(1,1) = .08;
% loc(1,1) = 0.35;
% for ii = 2:4
%     loc(1,ii) = loc(1,ii-1) + spw(ii-1) + sep;
% end

%%%%%%%%%%%%%%%%%%%%%%
spw = [spw1 spw2, spw2, spw3, spw2, spw2, spw3, spw2];
loc = zeros(2,8);
loc(2,:) = .14;
loc(1,1) = .08;
for ii = 2:8
    loc(1,ii) = loc(1,ii-1) + spw(ii-1) + sep;
end
%%%%%%%%%%%%%%%%%%%%%

if add_to_vel_ref
    vel_ref2 = zeros( size(vel_ref,1)+1, size(vel_ref,2) );
    vel_ref2(1:size(vel_ref,1),:) = vel_ref;
    vel_ref2(end,:) = vel_ref(end,:);
    vel_ref2(end,1) = hmx;
else
    vel_ref2 = vel_ref;
end

% fig10 = figure;hold on; box on;
%set(fig10, 'renderer', 'opengl')
%h4 is for interfaces
%h1 is always for the first physical parameter in voro: sigma for (iseis=0) and Vs for (iseis=1)
%h2 is always for the second physical parameter in voro: Vp/Vs for
% (ivpvs=1) and sigma for (ivpvs~=1)
%h3 is always for the third physical parameter in voro: sigma
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]); %sigma for iseis=0, Vs for iseis=1
first_plot = 1;
hold on; box off;
h2 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]); hold on; box on; %%!!!!!!!!
set(h2,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(h2,'XDir','reverse','TickDir','out','XTickLabel',[]);
h3 = subplot('Position',[loc(1,4) loc(2,4) spw3 sph]); hold on; box on; %%!!!!!!!!
set(h3,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(h3,'XDir','reverse','TickDir','out','XTickLabel',[]);
if NPL>2
    if(ivpvs==1) 
        h2 = subplot('Position',[loc(1,3) loc(2,3) spw3 sph]); hold on; box off;%here h2 is for Vp/Vs
        h3 = subplot('Position',[loc(1,3)+spw3+sep loc(2,3) spw2 sph]); hold on; box on;%%!!!!!!!!!!!!
        set(h3,'YDir','reverse','TickDir','out','YTickLabel',[]);
        set(h3,'XDir','reverse','TickDir','out','XTickLabel',[]);
    else
        h1 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]); hold on; box off; %here h1 is for Vs
        h2 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]); hold on; box off; %here h2 is for sigma
        first_plot = 2;
        h3 = subplot('Position',[loc(1,4) loc(2,4) spw3 sph]); hold on; box on;  %%!!!!!!
        set(h3,'YDir','reverse','TickDir','out','YTickLabel',[]);
        set(h3,'XDir','reverse','TickDir','out','XTickLabel',[]);

    end
end
if(NPL>3);
  h1 = subplot('Position',[loc(1,3) loc(2,3) spw2 sph]); hold on; box off; %here h1 is for Vs
  h2 = subplot('Position',[loc(1,4) loc(2,4) spw3 sph]); hold on; box off; %here h2 is for Vp/Vs
  h3 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]); hold on; box off; %here h3 is for sigma
  first_plot = 3;
end;
if (plot_interfaces); h4 = subplot('Position',[loc(1,1) loc(2,1) spw1 sph]);hold on; box off;first_plot=4;end
%%%%%%%%%%%%%%%
h5 = subplot('Position',[loc(1,5) loc(2,5) spw2 sph]); hold on; box on;  %%!!!!!!
set(h5,'YDir','reverse','TickDir','in','YTickLabel',[]);
set(h5,'XDir','reverse','TickDir','in','XTickLabel',[]);
h6 = subplot('Position',[loc(1,6) loc(2,6) spw2 sph]); hold on; box on;  %%!!!!!!
set(h6,'YDir','reverse','TickDir','in','YTickLabel',[]);
set(h6,'XDir','reverse','TickDir','in','XTickLabel',[]);
h7 = subplot('Position',[loc(1,7) loc(2,7) spw3 sph]); hold on; box on;  %%!!!!!!
set(h7,'YDir','reverse','TickDir','in','YTickLabel',[]);
set(h7,'XDir','reverse','TickDir','in','XTickLabel',[]);
h8 = subplot('Position',[loc(1,8) loc(2,8) spw2 sph]); hold on; box on;  %%!!!!!!
set(h8,'YDir','reverse','TickDir','in','YTickLabel',[]);
set(h8,'XDir','reverse','TickDir','in','XTickLabel',[]);
%%%%%%%%%%%%%%%

if plot_interfaces
subplot(h4);hold on;box on;
if(imarg == 1)
  %h = h/sum(h);
  %[xx,yy]=stairs(h,zi,'-k');
  %xx = [0;xx;0];
  %yy = [yy(1);yy;yy(end)];
  %patch(xx,yy,[0.8,0.8,0.8]);
  %[xx,yy]=stairs(h,zi,'-k');

  [xx,yy]=stairs(Ndep,binsdep,'-k');
  patch(xx,yy,[0.8,0.8,0.8]);
  %[xx,yy]=stairs(h,z,'-k');
  set(gca,'Fontsize',font_size,'YLim',[zmin hmax2],'XLim',xxlim);
  set(gca,'YDir','reverse','TickDir','out');
%   xlabel('Interface probability');
  xlabel({'Interface', 'probability'});
  ylabel('Depth (km)');
  
  if isyn == 1
      for in = 2:ktru
          %if in == 3; continue; end
      plot(xxlim, [vorotru(in,1) vorotru(in,1)], '-w','LineWidth', lwt)
      plot(xxlim, [vorotru(in,1) vorotru(in,1)], '--k','LineWidth', lwt)
      end
      box on;
  end
end;
end

subplot(h1);
set(gca,'FontSize',font_size);
if(imarg == 1)
   pcolor(clim,z,Nc);shading flat;
   if(imead == 1);plot(c_mead,z,'.r','Linewidth',1);end;
   if(imean == 1);plot(c_mean,z,'.r','Linewidth',.5);end;
end;
%surf(clim,z,Nc);shading flat;
set(h1,'layer','top');
set(gca,'Fontsize',font_size,'XLim',[pmin2(1) pmax2(1)],'YLim',[zmin hmax2]);
set(gca,'YDir','reverse','TickDir','out');
if first_plot~=1
    set(gca,'YTickLabel',[]);
else
     ylabel('Depth (km)');
end
if(isyn == 1)
   h1=plprof(mtru,hmax,NPL,'-w',1,0);
   set(h1,'LineWidth',lwt);
   h1=plprof(mtru,hmax,NPL,'--k',1,0);
   set(h1,'LineWidth',lwt);
end;
if irv == -1 || irv == 1  || iswd == 1
%if 1 == 2    
if(ivref == 1)
    if (plot_prior_Vs)
   %h1=plot(vel_ref(:,2),vel_ref(:,1),'-w','LineWidth',2);
   %h1=plot(vel_ref(:,2),vel_ref(:,1),'--k','LineWidth',2);

   h1=plot(vel_ref2(:,2)+dVsp,vel_ref2(:,1),'-w','LineWidth',lwp);
   h1=plot(vel_ref2(:,2)+dVsp,vel_ref2(:,1),'--k','LineWidth',lwp);

   h1=plot(vel_ref2(:,2)+dVsm,vel_ref2(:,1),'-w','LineWidth',lwp);
   h1=plot(vel_ref2(:,2)+dVsm,vel_ref2(:,1),'--k','LineWidth',lwp);
    end
end;
%end
else
   if (plot_prior_res)
   h1=plot(pmin(end)*ones(size(z)),z,'-w','LineWidth',lwp);
   h1=plot(pmin(end)*ones(size(z)),z,'--k','LineWidth',lwp);

   h1=plot(pmax(end)*ones(size(z)),z,'-w','LineWidth',lwp);
   h1=plot(pmax(end)*ones(size(z)),z,'--k','LineWidth',lwp);
   end
end
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
if iseis == 0
%     set(gca,'XTick',[-10:1:0]);
%     xlabel('log_{10}\sigma (S/m)')  
set(gca,'XTick',XTick_res);
    xlabel('log_{10}\rho (\Omega m)')  
else
    set(gca,'XTick',XTick_Vs);
    xlabel('V_S (km/s)')
end
xtickangle(0)
box on;
if NPL>2
subplot(h2)
set(gca,'FontSize',font_size);
if(imarg == 1)
  pcolor(rlim,z,Nr);shading flat;
  if(imead == 1);plot(r_mead,z,'.r','Linewidth',1);end;
  if(imean == 1);plot(r_mean,z,'.r','Linewidth',.5);end;
end;
set(h2,'layer','top')
set(gca,'Fontsize',font_size,'XLim',[pmin2(2) pmax2(2)],'YLim',[zmin hmax2]);
set(gca,'YDir','reverse','TickDir','out');
if first_plot~=2
    set(gca,'YTickLabel',[]);
else
     ylabel('Depth (km)');
end
if ivpvs == 1
    set(gca,'XTick',XTick_VpVs);
    xlabel('V_P/V_S ratio');
else
%     set(gca,'XTick',[-10:1:0]);
%     xlabel('log_{10}\sigma (S/m)') 
    set(gca,'XTick',XTick_res);
    xlabel('log_{10}\rho (\Omega m)') 
end
xtickangle(0)
%if irv == -1 || irv == 1  || iswd == 1
%if 1 == 2
if(ivref == 1 && ivpvs == 1)
    if (plot_prior_VpVs)
   %h1=plot(vel_ref(:,2),vel_ref(:,1),'-w','LineWidth',2);
   %h1=plot(vel_ref(:,2),vel_ref(:,1),'--k','LineWidth',2);
   
   h2=plot(vel_ref2(:,3)+dVpVsp,vel_ref2(:,1),'-w','LineWidth',lwp);
   h2=plot(vel_ref2(:,3)+dVpVsp,vel_ref2(:,1),'--k','LineWidth',lwp);

   h2=plot(vel_ref2(:,3)+dVpVsm,vel_ref2(:,1),'-w','LineWidth',lwp);
   h2=plot(vel_ref2(:,3)+dVpVsm,vel_ref2(:,1),'--k','LineWidth',lwp);
    end
else
   if (plot_prior_res)
   h2=plot(pmin(end)*ones(size(z)),z,'-w','LineWidth',lwp);
   h2=plot(pmin(end)*ones(size(z)),z,'--k','LineWidth',lwp);

   h2=plot(pmax(end)*ones(size(z)),z,'-w','LineWidth',lwp);
   h2=plot(pmax(end)*ones(size(z)),z,'--k','LineWidth',lwp);
   end
end
%end
%end
if(isyn == 1)
   h2=plprof(mtru,hmax,NPL,'-w',2,0);
   set(h2,'LineWidth',lwt);
   h2=plprof(mtru,hmax,NPL,'--k',2,0);
   set(h2,'LineWidth',lwt);
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
end
if(NPL>3);
%   subplot(h3)
%   set(gca,'FontSize',14);
%   if(imarg == 1)
%     pcolor(alim,z,Ndip);shading flat;
%     if(imead == 1);plot(a_mead,z,'.k','Linewidth',2);end;
%     if(imean == 1);plot(a_mean,z,'.k','Linewidth',2);end;
%   end;
%   %surf(alim,z,Ndip);shading flat;
%   set(h3,'layer','top')
%   set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax]);
%   set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
%   xlabel('Dip (deg.)');
%   if(isyn == 1)
%     h1=plot(vorotru(:,4),vorotru(:,1),'.k');
%     h1=plot(vorotru(:,4),vorotru(:,1),'ow');
%   end;
%   if(imap == 1)
%     plprof(mmap,hmax,NPL,'--k',3);
%   end;
%   if(imax == 1)
%     for i=1:length(z);
%       [amx,j] = max(Ndip(i,:));
%       a_max(i) = alim(j);
%     end;
%     plot(a_max,z,'w');
%     save('max_model.mat','z','c_max','r_max','a_max');
%   end;
% %  cmap = colormap(jet);
% %  colormap(cmap);
%   box on;
subplot(h3);
set(gca,'FontSize',font_size);
if(imarg == 1)
   pcolor(alim,z,Ndip);shading flat;
   if(imead == 1);plot(a_mead,z,'.r','Linewidth',1.);end;
   if(imean == 1);plot(a_mean,z,'.r','Linewidth',.5);end;
end;
%surf(clim,z,Nc);shading flat;
set(h3,'layer','top');
set(gca,'Fontsize',font_size,'XLim',[pmin2(3) pmax2(3)],'YLim',[zmin hmax2]);
set(gca,'YDir','reverse','TickDir','out');
if first_plot~=3
    set(gca,'YTickLabel',[]);
else
     ylabel('Depth (km)');
end
if(isyn == 1)
   h3=plprof(mtru,hmax,NPL,'-w',3,0);
   set(h3,'LineWidth',lwt);
   h3=plprof(mtru,hmax,NPL,'--k',3,0);
   set(h3,'LineWidth',lwt);
end;
% if irv == -1 || irv == 1  || iswd == 1
% if(ivref == 1)
   % h1=plot(vel_ref(:,2),vel_ref(:,1),'-w','LineWidth',2);
   % h1=plot(vel_ref(:,2),vel_ref(:,1),'--k','LineWidth',2);
 
    %h1=plot(vel_ref(:,2)+dVs,vel_ref(:,1),'-w','LineWidth',2);
    %h1=plot(vel_ref(:,2)+dVs,vel_ref(:,1),'--k','LineWidth',2);
 
    %h1=plot(vel_ref(:,2)-dVs,vel_ref(:,1),'-w','LineWidth',2);
    %h1=plot(vel_ref(:,2)-dVs,vel_ref(:,1),'--k','LineWidth',2);
 %end;
 %end
 if (plot_prior_res)
 h3=plot(pmin(end)*ones(size(z)),z,'-w','LineWidth',lwp); %%%%!!!!!!!!!!!!!!!!!!!
 h3=plot(pmin(end)*ones(size(z)),z,'--k','LineWidth',lwp);

 h3=plot(pmax(end)*ones(size(z)),z,'-w','LineWidth',lwp);
 h3=plot(pmax(end)*ones(size(z)),z,'--k','LineWidth',lwp);
 end
if(imap == 1)
   plprof(mmap,hmax,NPL,'--k',1,0);
end;
if(imax == 1)
  for i=1:length(z);
     [amx,j] = max(Ndip(i,:));
     a_max(i) = alim(j);
  end;
  plot(a_max,z,'w');
end;
% set(gca,'XTick',[-10:1:0]);
% set(gca,'XScale','log');
% xlabel('log_{10}\sigma (S/m)')
set(gca,'XTick',XTick_res);
%set(gca,'XScale','log');
xlabel('log_{10}{\rho} (\Omega m)')
xtickangle(0)
box on;
end;
%suptitle('RF - SWD - MT')
%%
%% Interface PLOT
%%
iinterface = 0;
npanely = 3;
if iinterface==1
fig20=figure;
deplim = [0. 11. 51. 220.];
color = [1. 0. 0.;...
         0. 1. 0.;...
         0. 0. 1.];
% deplim = [0. 220.];
Pmax = [.8 .25 .065];
font_size = 25;
figw = 12;
figh = 6;
set(fig20,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = length(deplim) - 1;
ny = NSD2;
xim = 0.01;
yim = 0.05;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
% spw = spw - 0.005;
spw2 = spw;
spw = spw * 2./3.;
for ipanel=1:nx
    for ipanely=2:ny
        loc(1,ipanely+(ipanel-1)*ny) = loc(1,ipanely+(ipanel-1)*ny) - (ipanely-1)*(spw2-spw);
    end
end
% XTick = [0:5:10];

for ipanel=1:nx
    for ipanely = 1:ny
    
      h1 = subplot('Position',[loc(1,ipanely+(ipanel-1)*ny) loc(2,ipanely+(ipanel-1)*ny) spw sph]);
      hold on; box on;
      if ipanely==npanely
          [xx,yy]=stairs(binsdep,Ndep,'-k');
          patch(xx,yy,[.8 .8 .8]);
%           patch(xx,yy,color(ipanel,:));
          for iiface = 1:ipanel-1
              xlims = ( binsdep>=deplim(iiface) & binsdep<=deplim(iiface+1) );
              binsdep2 = binsdep(xlims);
              Ndep2 = Ndep(xlims);
              Ndep2 = [0,Ndep2,0];
              binsdep2 = [binsdep2(1),binsdep2,binsdep2(end)];
              [xx2,yy2]=stairs(binsdep2,Ndep2,'-k');
%               patch(xx2,yy2,color(iiface,:));
               patch(xx2,yy2,[.4 .4 .4]);
          end
        %   if isyn == 1
        %       for in = 2:ktru
        %           %if in == 3; continue; end
        %       plot(xxlim, [vorotru(in,1) vorotru(in,1)], '-w','LineWidth', lwt)
        %       plot(xxlim, [vorotru(in,1) vorotru(in,1)], '--k','LineWidth', lwt)
        %       end
        %       box on;
        %   end
      end
    
      set(gca,'Fontsize',font_size)
      set(gca,'XLim',[deplim(1) deplim(ipanel+1)])
%       set(gca,'XLim',[deplim(ipanel) deplim(ipanel+1)])
    %   set(gca,'XLim',deplim)
%       if (ipanel<nx); set(gca,'XTickLabel',[]);end
      set(gca,'YLim',[0., Pmax(ipanel)]);
      if (ipanely>1); set(gca,'YTickLabel',[]);end
      set(gca,'TickDir','in');
      if (ipanely==1); ylabel({'Interface', 'probability'});end;
      if (ipanel==nx); xlabel('Depth (km)'); end;

    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Node density of voronoi cells:
%%

if(ivoro >= 1);
fig11 = figure;hold on; box on;
set(fig11, 'renderer', 'painters')
h1 = subplot('Position',[loc(1,2) loc(2,2) spw2 sph]);
hold on; box off;
if NPL>2
h2 = subplot('Position',[loc(1,3) loc(2,3) spw3 sph]);
hold on; box off;
end
if(NPL>3);
  h3 = subplot('Position',[loc(1,4) loc(2,4) spw3 sph]);
  hold on; box off;
end;

subplot(h1);
set(gca,'FontSize',14,'layer','top','LineWidth',1);
if(ivref == 0)
  set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(Vs(:,2),Vs(:,1),[pmin(1) pmax(1) 0 hmax],true,[1001 1001]);
  set(gca,'XTick',[1:1:10]);
  xlabel('V_S (m/s)');
else
  set(gca,'Fontsize',14,'XLim',[dVsm dVsp],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(Vs(:,2),Vs(:,1),[dVsm dVsp 0 hmax],true,[1001 1001]);
  set(gca,'XTick',[dVsm:.2:dVsp]);
  xlabel('dV_S (m/s)');
end;
cmap = colormap( jet(1000) );
cmap(1,:) = [1 1 1];
colormap(cmap);
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out');
set(gca,'CLim',[0 1.5],'FontSize',14)
ylabel('Depth (km)');
box on;
if NPL>2
subplot(h2)
set(gca,'FontSize',14,'layer','top','LineWidth',1)
set(gca,'Fontsize',14,'XLim',[pmin(2) pmax(2)],'YLim',[0 hmax2]);
if(ivref == 0)
  set(gca,'Fontsize',14,'XLim',[pmin(1) pmax(1)],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(VpVs(:,2),VpVs(:,1),[pmin(2) pmax(2) 0 hmax],true,[1001 1001]);
  set(gca,'XTick',[1:1:10]);
  xlabel('V_P/V_S ratio');
else
  set(gca,'Fontsize',14,'XLim',[dVsm dVsp],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(VpVs(:,2),VpVs(:,1),[dVpVsm dVpVsp 0 hmax],true,[1001 1001]);
  set(gca,'XTick',[dVpVsm:2*dVpVsp/5:dVpVsp]);
  xlabel('dVpVs');
end;
%colormap( 1-hot );
%colormap( 1-gray(256) );
set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
set(gca,'CLim',[0 1.5],'FontSize',14)
box on;
end
if(NPL>3);
  subplot(h3)
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  set(gca,'Fontsize',14,'XLim',[pmin(3) pmax(3)],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(dip(:,2),dip(:,1),[pmin(3) pmax(3) 0 hmax],true,[1001 1001]);
  %colormap( 1-hot );
  %colormap( 1-gray(256) );
  set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
  set(gca,'CLim',[0 1.5],'FontSize',14)
  xlabel('Dip (deg.)');
  set(gca,'XTick',[0:5:45]);
  box on;
end;

if(NPL>4);
  subplot(h5)
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  set(gca,'Fontsize',14,'XLim',[pmin(4) pmax(4)],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(strike(:,2),strike(:,1),[pmin(4) pmax(4) 0 hmax],true,[1001 1001]);
  %colormap( 1-hot );
  %colormap( 1-gray(256) );
  set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
  set(gca,'CLim',[0 1.5],'FontSize',14)
  xlabel('Strike angle (deg.)');
  set(gca,'XTick',[-15:1:0]);
  box on;
end;

if(NPL>5);
  subplot(h6)
  set(gca,'FontSize',14,'layer','top','LineWidth',1)
  set(gca,'Fontsize',14,'XLim',[pmin(5) pmax(5)],'YLim',[0 hmax2]);
  [hdens]=cloudPlot(Ks(:,2),Ks(:,1),[pmin(5) pmax(5) 0 hmax],true,[1001 1001]);
  %colormap( 1-hot );
  %colormap( 1-gray(256) );
  set(gca,'YDir','reverse','TickDir','out','YTickLabel',[]);
  set(gca,'CLim',[0 1.5],'FontSize',14)
  xlabel('log(Ks)');
  set(gca,'XTick',[0:5:25]);
  box on;
end;

end;

if(idatfit == 1);
  %%
  %%
  figw = 10;
  figh = 8;
  fig12=figure();
  set(fig12,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
  %nx = 4;
  %ny = NRF;
  nx = 1;
  ny = 3;
  xim = 0.00;
  yim = 0.08;
  xymarg = [0.1 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  if(idathist == 1);
    NV = 3000;
    %%
    %% Vertical
    %%
    isp = 1;
    for iaz=1:NRF
      %if(irf == 0);
      if irv >= 1
      vmin = min(min(min(predV)))-0.1;
      vmax = max(max(max(predV)))+0.1;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for idat=1:NTIME
        Nd(:,idat) = histc(predV(:,iaz,idat),vlim);
      end;
      %% Normalize Histograms
      for idat=1:NTIME
        Nd(:,idat) = Nd(:,idat)/max(Nd(:,idat));
      end;
      h2 = subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);
      isp = isp + 1;
      hold on; 
      imagesc(dt*[0:NTIME-1],vlim,Nd);shading flat;
      clear Nd,vlim;
      p2=plot(dt*[0:NTIME2-1],dobsV(iaz,:),'-w','LineWidth',1);
      p2=plot(dt*[0:NTIME2-1],dobsV(iaz,:),'--k','LineWidth',1);
      cmap = colormap(jet(1000));
      cmap(1,:) = [1 1 1];
      colormap(cmap);
      xlabel('Time (s)');ylabel('V component');
      %legend([p2],'data');
      %set(gca,'XLim',[0 dt*NTIME],'TickDir','out');
      set(gca,'XLim',[0 35],'TickDir','out');
      set(gca,'YLim',[min(min(min(predV)))-0.05,max(max(max(predV)))+0.05],'TickDir','out');
      box on;set(gca,'FontSize',14);
      end;

      %%
      %% Radial
      %%
      vmin = min(min(min(min(predR))))-0.1;
      vmax = max(max(max(max(predR))))+0.1;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for idat=1:NTIME
        Nd(:,idat) = histc(predR(:,iaz,idat),vlim);
      end;
      %% Normalize Histograms
      for idat=1:NTIME
        Nd(:,idat) = Nd(:,idat)/max(Nd(:,idat));
      end;
      h1 = subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);
      isp = isp + 1;
      hold on; 
      imagesc(dt*[0:NTIME-1],vlim,Nd);shading flat;
      clear Nd,vlim;
      p2=plot(dt*[0:NTIME2-1],dobsR(iaz,:),'-w','LineWidth',1);
      p2=plot(dt*[0:NTIME2-1],dobsR(iaz,:),'--k','LineWidth',1);
      cmap = colormap(jet(1000));
      cmap(1,:) = [1 1 1];
      colormap(cmap);
      xlabel('Time (s)');ylabel('R component');
      %legend([p2],'data');
      %set(gca,'XLim',[0 dt*NTIME],'TickDir','out');
      set(gca,'XLim',[0 35],'TickDir','out');
      set(gca,'YLim',[min(min(min(predR)))-0.05,max(max(max(predR)))+0.05],'TickDir','out');
      box on;set(gca,'FontSize',14);

      %%
      %% Transverse
      %%
      if(itr == 1);
      vmin = min(min(min(predT)))-0.1;
      vmax = max(max(max(predT)))+0.1;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for idat=1:NTIME
        Nd(:,idat) = histc(predT(:,iaz,idat),vlim);
      end;
      %% Normalize Histograms
      for idat=1:NTIME
        Nd(:,idat) = Nd(:,idat)/max(Nd(:,idat));
      end;
      h1 = subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);
      isp = isp + 1;
      hold on; 
      imagesc(dt*[0:NTIME-1],vlim,Nd);shading flat;
      clear Nd,vlim;
      p2=plot(dt*[0:NTIME2-1],dobsT(iaz,:),'-w','LineWidth',1);
      p2=plot(dt*[0:NTIME2-1],dobsT(iaz,:),'--k','LineWidth',1);
      %xlabel('Time (s)');ylabel('T component');
      %legend([p2],'data');
      %set(gca,'XLim',[0 dt*NTIME],'TickDir','out');
      set(gca,'XLim',[0 35],'TickDir','out');
      set(gca,'YLim',[min(min(min(predT)))-0.05,max(max(max(predT)))+0.05],'TickDir','out');
      box on;set(gca,'FontSize',14);
      end;

      %if(irf == 0);
      if irv >= 1
      %%
      %% Source
      %%
      vmin = min(min(min(predS)))-0.1;
      vmax = max(max(max(predS)))+0.1;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for idat=1:NSRC
        Nd(:,idat) = histc(predS(:,iaz,idat),vlim);
      end;
      %% Normalize Histograms
      for idat=1:NSRC
        Nd(:,idat) = Nd(:,idat)/max(Nd(:,idat));
      end;
      %h3 = subplot('Position',[loc(1,3) loc(2,3) 1.5*spw*(NSRC/NTIME2) sph]);
      h3 = subplot('Position',[loc(1,isp) loc(2,isp) 0.5*spw sph]);
      isp = isp + 1;
      hold on;
      imagesc(dt*[0:NSRC-1],vlim,Nd);shading flat;
      clear Nd,vlim;
      if(isyn == 1);
        p2=plot(dt*[0:NSTFtru-1]+1,stftru,'-w','LineWidth',1);
        p2=plot(dt*[0:NSTFtru-1]+1,stftru,'--k','LineWidth',1);
      end;
      xlabel('Time (s)');ylabel('Source estimate');
      %legend([p2],'data');
      set(gca,'XLim',[0 dt*NSRC],'TickDir','out');
      set(gca,'YLim',[min(min(min(predS)))-0.05,max(max(max(predS)))+0.05],'TickDir','out');
      box on;set(gca,'FontSize',14);
      end;
    end;

    if(iswd == 1);
    %%
    %% SWD
    %%
    vmin = min(min(predSWD))-0.1;
    vmax = max(max(predSWD))+0.1;
    vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
    for idat=1:NDAT_SWD
      Nd(:,idat) = histc(predSWD(:,idat),vlim);
    end;
    %% Normalize Histograms
    for idat=1:NDAT_SWD
      Nd(:,idat) = Nd(:,idat)/max(Nd(:,idat));
    end;
    Nd2=[Nd,zeros(size(Nd,1),1)];
    Nd2=[Nd2;zeros(1,size(Nd2,2))];
    vlim2=[vlim,vlim(end)];
    swdphase2=[swdphase,swdphase(1,end)+(swdphase(1,end)-swdphase(1,end-1))];
    dph = diff(swdphase2)/2;
    swdphase2(1:21)=swdphase(1:21)-dph;
    swdphase2(22)=swdphase2(21)+2*dph(end);
    fig13 = figure();
    hold on; box on;set(gca,'FontSize',14);
    pcolor(swdphase2,vlim2,Nd2);shading flat;
    clear Nd,vlim;
    p2=plot(swdphase,dobsSWD,'-ow','LineWidth',1);
    p2=plot(swdphase,dobsSWD,'--k','LineWidth',1);
    %if(isyn == 1);
    %  p2=plot(swdphase,swdtru,'-w','LineWidth',1);
    %  p2=plot(swdphase,swdtru,'--k','LineWidth',1);
    %end;
    xlabel('Perid (s)');ylabel('Phase velocity');
    set(gca,'XLim',[swdphase(1) swdphase(end)],'TickDir','out');
    set(gca,'YLim',[min(min(predSWD))-0.05,max(max(predSWD))+0.05],'TickDir','out');
    end;

  else;
    h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
    hold on; box on;set(gca,'FontSize',14);
    for ipred=1:NPRED;
      plot(dt*[0:NTIME-1],predR(ipred,:),'LineWidth',0.5)
    end;
    plot(dt*[0:NTIME2-1],dobsR,'xr','Markersize',4);
    ylabel('H component');
    set(gca,'XLim',[0 dt*NTIME2],'TickDir','out','XTickLabel',[]);

    h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
    hold on; box on;set(gca,'FontSize',14);
    for ipred=1:NPRED;
      p1=plot(dt*[0:NTIME-1],predV(ipred,:),'LineWidth',0.5);
    end;
    p2=plot(dt*[0:NTIME2-1],dobsV,'xr','Markersize',4);
    xlabel('Time (s)');ylabel('V component');
    legend([p2,p1],'data','estimate');
    set(gca,'XLim',[0 dt*NTIME2],'TickDir','out');

    h3 = subplot('Position',[loc(1,3) loc(2,3) spw*(NSRC/NTIME2) sph]);
    hold on; box on;set(gca,'FontSize',14);
    for ipred=1:NPRED;
      plot(dt*[0:NSRC-1],predS(ipred,1:NSRC),'LineWidth',0.5);
    end;
    plot(dt*[0:NSTFtru-1],stftru(1:NSTFtru),'-r','LineWidth',1);
    xlabel('Time (s)');ylabel('Source-time function');
    set(gca,'XLim',[0 dt*NSRC],'TickDir','out');
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% PLOT res hist 
  %%
  fig14=figure();
  figw = 9;
  figh = 5;
  set(fig14,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
  nx = 3;
  ny = NRF;
  xim = 0.04;
  yim = 0.04;
  xymarg = [0.1 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

  x = -8.:.5:8.;
  xx = -8.:.01:8.;
  isp = 1;
  for iaz=1:NRF
    for i=1:size(predR,1);
      %% remove mean
      if(irv==-1);resR(i,iaz,:) = resR(i,iaz,:)-mean(resR(i,iaz,:)); end; %mine if
      %if(irf == 0);
      if irv >= 1
      resV(i,iaz,:) = resV(i,iaz,:)-mean(resV(i,iaz,:));end;
      if(itr == 1);
      resT(i,iaz,:) = resT(i,iaz,:)-mean(resT(i,iaz,:));end;
      %% detrend
      if(irv==-1);resR(i,iaz,:) = detrend(squeeze(resR(i,iaz,:))); end; %mine if
      %if(irf == 0);
      if irv >= 1
      resV(i,iaz,:) = detrend(squeeze(resV(i,iaz,:)));end;
      if(itr == 1);
      resT(i,iaz,:) = detrend(squeeze(resT(i,iaz,:)));end;
      %% Radial
      if irv == -1 %mine
      [n1H(i,iaz,:),xoutR] = hist(squeeze(resR(i,iaz,:))./std(squeeze(resR(i,iaz,:))),x);
      narea = sum(n1H(i,iaz,:)) * (xoutR(2)-xoutR(1));
      n1H(i,iaz,:) = n1H(i,iaz,:)/narea;
      end
      %% Vertical
      %if(irf == 0);
      if irv >= 1
      [n1V(i,iaz,:),xoutV] = hist(squeeze(resV(i,iaz,:))./std(squeeze(resV(i,iaz,:))),x);
      narea = sum(n1V(i,iaz,:)) * (xoutV(2)-xoutV(1));
      n1V(i,iaz,:) = n1V(i,iaz,:)/narea;end;
      %% Transverse:
      if(itr == 1);
      [n1T(i,iaz,:),xoutT] = hist(squeeze(resT(i,iaz,:))./std(squeeze(resT(i,iaz,:))),x);
      narea = sum(n1T(i,iaz,:)) * (xoutT(2)-xoutT(1));
      n1T(i,iaz,:) = n1T(i,iaz,:)/narea;end;
    end;
  end;
  if(iswd == 1);
    for i=1:size(predR,1);
      resSWD(i,:) = resSWD(i,:)-mean(resSWD(i,:));
      resSWD(i,:) = detrend(resSWD(i,:));
      %% SWD
      [n1S(i,:),xoutS] = hist(resSWD(i,:)./std(resSWD(i,:)),x);
      narea = sum(n1S(i,:)) * (xoutS(2)-xoutS(1));
      n1S(i,:) = n1S(i,:)/narea;
    end;
  end;

  nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
  %plot(xx,nd,'-k');
  %% Plot all histograms as stairs (stairs takes 2D matrix): 
  %xout = xout - (xout(2)-xout(1))/2;
  %stairs(xout,n1','k');
  %axis([-6 6 0 0.6]);
  %set(gca,'XTick',[-4 0 4],'FontSize',14);
  %set(gca,'XTickLabel',[],'FontSize',14);
  yyH = 0:max(max(max(n1H)))/50:max(max(max(n1H)));
  %if(irf == 0);
  if irv >= 1
  yyV = 0:max(max(max(n1V)))/50:max(max(max(n1V)));end;
  if(itr  == 1);
  yyT = 0:max(max(max(n1T)))/50:max(max(max(n1T)));end;
  isp = 1;
  for iaz=1:NRF
    for i=1:length(x);
      %% Radial
      if irv == -1 %mine
      [n2R(:,iaz,i),xout2R] = hist(squeeze(n1H(:,iaz,i)),yyH);
      narea = sum(squeeze(n2R(:,iaz,i))) * (xout2R(2)-xout2R(1));
      n2R(:,iaz,i) = n2R(:,iaz,i)/narea;
      end
      %% Vertical
      %if(irf == 0);
      if irv >= 1
      [n2V(:,iaz,i),xout2V] = hist(squeeze(n1V(:,iaz,i)),yyV);
      narea = sum(n2V(:,iaz,i)) * (xout2V(2)-xout2V(1));
      n2V(:,iaz,i) = n2V(:,iaz,i)/narea;end;
      %% Transverse
      if(itr == 1);
      [n2T(:,iaz,i),xout2T] = hist(squeeze(n1T(:,iaz,i)),yyT);
      narea = sum(n2T(:,iaz,i)) * (xout2T(2)-xout2T(1));
      n2T(:,iaz,i) = n2T(:,iaz,i)/narea;end;
    end;
    %% Radial
    if irv == -1 %mine
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    set(gca,'FontSize',14);
    xout2R = xout2R - (xout2R(2)-xout2R(1))/2;
    x2 = x - (x(2)-x(1))/2;
    pcolor(x2,xout2R,squeeze(n2R(:,iaz,:)));shading flat;
    plot(xx,nd,'-w');plot(xx,nd,'--k');
    axis([-4 4 0 max(max(max(n1H)))]);
    xlabel('Standard deviation')
    end
    %% Vertical
    %if(irf == 0);
    if irv >= 1
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    set(gca,'FontSize',14);
    xout2V = xout2V - (xout2V(2)-xout2V(1))/2;
    pcolor(x2,xout2V,squeeze(n2V(:,iaz,:)));shading flat;
    plot(xx,nd,'-w');plot(xx,nd,'--k');
    axis([-4 4 0 max(max(max(n1V)))]);
    xlabel('Standard deviation');end;
    %% Transverse
    if(itr == 1);
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    set(gca,'FontSize',14);
    xout2T = xout2T - (xout2T(2)-xout2T(1))/2;
    pcolor(x2,xout2T,squeeze(n2T(:,iaz,:)));shading flat;
    plot(xx,nd,'-w');plot(xx,nd,'--k');
    axis([-4 4 0 max(max(max(n1V)))]);
    xlabel('Standard deviation');end;
  end;

  fig15=figure();
  if(iswd == 1);
  subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
  yyS = 0:max(max(n1S))/13:max(max(n1S));end;
  for i=1:length(x);
    %% SWD
    if(iswd == 1);
    [n2S(:,i),xout2S] = hist(n1S(:,i),yyS);
    narea = sum(n2S(:,i)) * (xout2S(2)-xout2S(1));
    n2S(:,i) = n2S(:,i)/narea;end;
  end;
  %% SWD
  if(iswd == 1);
  %subplot('Position',[loc(1,3) loc(2,3) spw sph]);
  hold on;box on;
  set(gca,'FontSize',14);
  xout2S = xout2S - (xout2S(2)-xout2S(1))/2;
  pcolor(x2,xout2S,n2S);shading flat;
  plot(xx,nd,'-w');plot(xx,nd,'--k');
  axis([-4 4 0 max(max(n1S))]);
  xlabel('Standard deviation');end;

  %save tmp.mat
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% PLOT res Axx
  %%
  if 1==2 %mine
  fig16=figure();
  figw = 9;
  figh = 5;
  set(fig15,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
  nx = 3;
  ny = NRF;
  xim = 0.04;
  yim = 0.04;
  xymarg = [0.1 0.04 0.04 0.14];
  [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
  for iaz=1:NRF
    for i=1:size(predR,1);
      %% Radial
      [AxxH(i,iaz,:),lagsR] = xcov(resR(i,iaz,:),'coeff');
      %% Vertical
      %if(irf == 0);
      if irv >= 1
      [AxxV(i,iaz,:),lagsV] = xcov(resV(i,iaz,:),'coeff');end;
      %% Transverse
      if(itr == 1);
      [AxxT(i,iaz,:),lagsT] = xcov(resT(i,iaz,:),'coeff');end;
    end;
  end;
  yyH = min(min(min(AxxH))):(1-min(min(min(AxxH))))/50:1.;
  %if(irf == 0);
  if irv >= 1
  yyV = min(min(min(AxxV))):(1-min(min(min(AxxV))))/50:1.;end;
  if(itr == 1);
  yyT = min(min(min(AxxT))):(1-min(min(min(AxxT))))/50:1.;end;
  isp = 1;
  for iaz=1:NRF
    for i=1:length(lagsR);
      %% Radial
      [n2R(:,iaz,i),xout2R] = hist(AxxH(:,iaz,i),yyH);
      narea = sum(n2R(:,iaz,i)) * (xout2R(2)-xout2R(1));
      n2R(:,iaz,i) = n2R(:,iaz,i)/narea;
      %% Vertical
      %if(irf == 0);
      if irv >= 1
      [n2V(:,iaz,i),xout2V] = hist(AxxV(:,iaz,i),yyV);
      narea = sum(n2V(:,iaz,i)) * (xout2V(2)-xout2V(1));
      n2V(:,iaz,i) = n2V(:,iaz,i)/narea;end;
      %% Transverse
      if(itr == 1);
      [n2T(:,iaz,i),xout2T] = hist(AxxT(:,iaz,i),yyT);
      narea = sum(n2T(:,iaz,i)) * (xout2T(2)-xout2T(1));
      n2T(:,iaz,i) = n2T(:,iaz,i)/narea;end;
    end;
    %% Theoretical Axx for Gaussian noise:
    AxxT1 = xcov(randn(1,size(resR,3)),'coeff');
    %% Radial
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    set(gca,'FontSize',14);
    xout2R = xout2R - (xout2R(2)-xout2R(1))/2;
    x2 = x - (x(2)-x(1))/2;
    pcolor(sampling_dt*lagsR,xout2R,squeeze(n2R(:,iaz,:)));shading flat;
    plot(sampling_dt*lagsR,AxxT1,'-w');plot(sampling_dt*lagsR,AxxT1,'--k');
    axis([sampling_dt*lagsR(1) sampling_dt*lagsR(end) yyH(1) 1]);
    xlabel('lag (s)')
    %% Vertical
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    %if(irf == 0);
    if irv >= 1
    set(gca,'FontSize',14);
    xout2V = xout2V - (xout2V(2)-xout2V(1))/2;
    pcolor(sampling_dt*lagsV,xout2V,squeeze(n2V(:,iaz,:)));shading flat;
    plot(sampling_dt*lagsV,AxxT1,'-w');plot(sampling_dt*lagsV,AxxT1,'--k');
    axis([sampling_dt*lagsV(1) sampling_dt*lagsV(end) yyV(1) 1]);
    xlabel('lag (s)');end;
    %% Transverse
    subplot('Position',[loc(1,isp) loc(2,isp) spw sph]);hold on;box on;
    isp = isp + 1;
    if(itr == 1);
    set(gca,'FontSize',14);
    xout2T = xout2T - (xout2T(2)-xout2T(1))/2;
    pcolor(sampling_dt*lagsT,xout2T,squeeze(n2T(:,iaz,:)));shading flat;
    plot(sampling_dt*lagsT,AxxT1,'-w');plot(sampling_dt*lagsT,AxxT1,'--k');
    axis([sampling_dt*lagsT(1) sampling_dt*lagsT(end) yyT(1) 1]);
    xlabel('lag (s)');end;
  end;

  if(iswd == 1);
    for i=1:size(predR,1);
      %% SWD
      if(iswd == 1);
      [AxxS(i,:),lagsS] = xcov(resSWD(i,:),'coeff');end;
    end;
    AxxT2 = xcov(randn(1,size(resSWD,2)),'coeff');
    fig16=figure();
    yyS = min(min(AxxS)):(1-min(min(AxxS)))/50:1.;
    for i=1:length(lagsS);
      %% SWD
      [n2S(:,i),xout2S] = hist(AxxS(:,i),yyS);
      narea = sum(n2S(:,i)) * (xout2S(2)-xout2S(1));
      n2S(:,i) = n2S(:,i)/narea;
    end;
    %% SWD
    subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
    set(gca,'FontSize',14);
    xout2S = xout2S - (xout2S(2)-xout2S(1))/2;
    pcolor(lagsS,xout2S,n2S);shading flat;
    plot(lagsS,AxxT2,'-w');plot(lagsS,AxxT2,'--k');
    axis([lagsS(1) lagsS(end) yyS(1) 1]);
    xlabel('lag (s)');
  end;
  end %mine
end;


t2 = toc(t1);
disp('time before saving:');disp(t2);
if(isave == 1)
   saveas(fig1,strcat(plotfile1,plotext2),'png');
   if(IEPS == 1);
   print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');end;

   saveas(fig2,strcat(plotfile2,plotext2),'png');
   if(IEPS == 1);
   print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');end;

   saveas(fig3,strcat(plotfile3,plotext2),'png');
   if(IEPS == 1);
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');end;

   saveas(fig4,strcat(plotfile4,plotext2),'png');
   if(IEPS == 1);
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');end;

   if(iar == 1);
     saveas(fig5,strcat(plotfile5,plotext2),'png');
     if(IEPS == 1);
     print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');end;
   end;

   saveas(fig6,strcat(plotfile6,plotext2),'png');
   if(IEPS == 1);
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');end;

   saveas(fig7,strcat(plotfile7,plotext2),'png');
   if(IEPS == 1);
   print(fig7,'-painters','-r250',strcat(plotfile7,plotext3),'-depsc');end;

   saveas(fig8,strcat(plotfile8,plotext2),'png');
   if(IEPS == 1);
   print(fig8,'-painters','-r250',strcat(plotfile8,plotext3),'-depsc');end;

   saveas(fig9,strcat(plotfile9,plotext2),'png');
   if(IEPS == 1);
   print(fig9,'-painters','-r250',strcat(plotfile9,plotext3),'-depsc');end;

   saveas(fig10,strcat(plotfile10,plotext2),'png');
   if(IEPS == 1);
   print(fig10,'-painters','-r250',strcat(plotfile10,plotext3),'-depsc');end;

   saveas(fig11,strcat(plotfile11,plotext2),'png');
   if(IEPS == 1);
   print(fig11,'-painters','-r250',strcat(plotfile11,plotext3),'-depsc');end;

  if(idatfit == 1);
    saveas(fig12,strcat(plotfile12,plotext2),'png');
    if(IEPS == 1);
    print(fig12,'-painters','-r250',strcat(plotfile12,plotext3),'-depsc');end;
    if(iswd == 1);
      saveas(fig13,strcat(plotfile13,plotext2),'png');
      if(IEPS == 1);
      print(fig13,'-painters','-r250',strcat(plotfile13,plotext3),'-depsc');end;
    end;
    saveas(fig14,strcat(plotfile14,plotext2),'png');
    if(IEPS == 1);
    print(fig14,'-painters','-r250',strcat(plotfile14,plotext3),'-depsc');end;
    saveas(fig15,strcat(plotfile15,plotext2),'png');
    if(IEPS == 1);
    print(fig15,'-painters','-r250',strcat(plotfile15,plotext3),'-depsc');end;
  end;
end;
t2 = toc(t1);
disp('total time:');disp(t2);
!return;
