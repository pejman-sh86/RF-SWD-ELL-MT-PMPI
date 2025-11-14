%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Regression trans-D BIC optimization for starting models derived from 
%% linear finite fault inversion results. 
%% Jan Dettmer May 2015	dettmer.jan@gmail.com
%%
function []=ffi_tsunami_sample_startmodel(filebase);

%%
%% Input variables
%%
imapstart = 1; %% Start at pre-computed MAP model
imapstep  = 0; %% Start with previous step size
ibic  = 0;     %% 0->logL; 1->BIC
inoise = 1;    %% Add noise to simulated model?

kmin =  1;
kmax = 60;
amp_min = -3.;
amp_max = 10.;
%%
%% Annealing variables:
%%
NT1   =   1;   %% No. chains at T=1
NT    =  10;   %% Total No. PT chains
bmin  =  -2.;  %% Lower limit of tempering range
bmax  =   0.;  %% upper limit of tempering range (should always be zero)
iexchange = 0;if(NT>1);iexchange = 1;end;
NMCMC = 1e8;   %% Max No. MCMC steps
NTHIN =   5;   %% Amount of chain thinning
NKEEP =  20;   %% Save sample every NKEEP MCMC steps

myCluster = parcluster('local');
myCluster.NumWorkers = NT; % 'Modified' property now TRUE
saveProfile(myCluster);    % 'local' profile now updated,
                           % 'Modified' property now FALSE

parfile  = strcat(filebase,'_parameter.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEPtmp,IADAPT,NBUF,...
 MAXDISP,MINDISP]=ffi_tsunami_read_parfile(parfile);

NSTN    = h5readatt(datfile,'/Observed_data','N_sta');
deltt   = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc = h5read(datfile,'/Observed_data/Hypo_Loc_Cart');
mat_excl= h5read(datfile,'/Sensitivity_kernel/mat_excl');
mat_excl = fliplr(mat_excl);
dhyp    = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP    = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx     = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx     = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin   = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax   = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');
NTW     = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
%NTW    = 1
NTW  = cast(NTW,'like',1);
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);
NPV = 2 + NTW;
%%
%% Set up grid
%%
deltr  = double(Rmx)/double(NRAN);
deltz  = double(Zmx)/double(NDEP);
deltrn = double(deltr)/double(Rmx);
deltzn = double(deltz)/double(Zmx);
x_evn = [deltrn/2.:deltrn:1.-deltrn/2.];
z_evn = [deltzn/2.:deltzn:1.-deltzn/2.];
x_ev = x_evn*Rmx;
z_ev = z_evn*Zmx;
z_ev = z_ev';
%%
%% No. data and exclusion array
%%
mat_excl2 = repmat(mat_excl,[NTW,1]);
NDAT = sum(sum(mat_excl2));
rex = mat_excl2(:);
%%
%% Find left and right exclusion indices:
%%
%% Here we need to find the last occurence of zero.
idx_exl = zeros(NDEP,1);
for iz = 1:NDEP;
  for ir = 1:NRAN/2;
    if(mat_excl(iz,ir) == 0)
      idx_exl(iz) = ir;
end;end;end;
idx_exr = (NRAN+1)*ones(NDEP,1);
%% Here we need to find the first occurence of zero (hence EXIT).
for iz = 1:NDEP
  for ir = NRAN/2:NRAN
    if(mat_excl(iz,ir) == 0)
      idx_exr(iz) = ir;
      break;
end;end;end;

%%
%% Priors:
%%
minlim(1:2)  = [ 0.,  0.];
maxlim(1:2)  = [Rmx, Zmx];
minlim(3:2+NTW) = -3.;
maxlim(3:2+NTW) = 10.;
maxpert = maxlim-minlim;
%%
%% Parameter structure
%%
param.NPV = NPV;
param.NRAN = NRAN;
param.NDEP = NDEP;
param.NTW = NTW;
param.NDAT = NDAT;
param.x_ev = x_ev;
param.z_ev = z_ev;
param.mat_excl = mat_excl;
param.minlim = minlim;
param.maxlim = maxlim;
param.maxpert = maxpert;
param.ibic = ibic; 
param.rex = rex;
param.idx_exl = idx_exl;
param.idx_exr = idx_exr;
param.kmin = kmin;
param.kmax = kmax;

%%
%% "Observed" image to fit
%%
sl1tru = load('true_sl1.txt');
sl1tru2 = zeros(NDEP,NRAN,NTW);
for itw = 1:NTW;
  sl1tru2(:,:,itw) = sl1tru((itw-1)*NDEP+1:itw*NDEP,:);
end;
obs = sl1tru2(:);  %% This reshapes into single row vector
if(inoise);
  obs = obs + 0.1*randn(size(obs));
end;

%%
%% Annealing schedule:
%%
if(bmax>bmin);
  beta = ones(1,NT);
  logbeta = [bmax-(bmax-bmin)/(NT-NT1-1):-(bmax-((bmax-bmin)/(NT-NT1-1))-bmin)/(NT-NT1-1):bmin];
  beta(NT1+1:NT) = 10.^logbeta;
else;
  beta = ones(1,NT);
end;

%%
%% Starting model
%%
for ic = 1:NT;
  obj(ic).voro       = zeros(kmax,2+NTW);
  obj(ic).step_sz    = zeros(kmax,2+NTW);
  obj(ic).iaccept    = zeros(kmax,2+NTW);
  obj(ic).ipropose   = zeros(kmax,2+NTW);
  obj(ic).iacceptbd  = 0;
  obj(ic).iproposebd = 0;
  if(imapstart == 0);
    idxtmp = randperm(kmax);
    obj(ic).k = idxtmp(1);
    %% MCMC start model
    ioutside = 1;
    while ioutside == 1;
      ioutside = 0;
      obj(ic).voro = repmat(param.minlim,[obj(ic).k 1]) + ...
                     repmat(param.maxpert,[obj(ic).k 1]).*rand(obj(ic).k,2+NTW);
      [ioutside]=checkbounds(obj(ic),ic,param);
    end;
  else
    load('mapstart.mat');
    %% MCMC start model
    %obj(ic) = mapobj;
    obj(ic).k = mapobj.k;
    obj(ic).voro = mapobj.voro;
  end;
  if(imapstep == 0);
    %% MCMC step size 
    obj(ic).step_sz(:,1:2) = repmat(maxpert(1:2)/500.,[kmax,1]);
    obj(ic).step_sz(:,3:2+NTW) = repmat(maxpert(3:2+NTW)/250.,[kmax,1]);
  else
    load('mapstart.mat');
    %% MCMC step size 
    obj(ic).step_sz = mapobj.step_sz;
  end;

  obj(ic).obs = obs;
  [obj(ic).BIC,obj(ic).sl1]=get_bic(obj(ic).voro(1:obj(ic).k,:),obj(ic).obs,param);
  [ioutside]=checkbounds(obj(ic),ic,param);
  %% Tempering level needs to be set after everything else...
  obj(ic).beta = beta(ic);
end;

if(ioutside);error('Starting model outside bounds, abort.');end;
save('sample.mat','beta','sl1tru2','minlim','maxlim','mat_excl','NDAT');

%%
%% Time the forward model:
%%
t1 = tic;
for i=1:10;
  [sl1]=voro_interp(obj(1).voro,param);
end;
time = toc(t1);
disp(['Time voro_interp:',num2str(time)]);

%%
%% Simulated annealing with BIC:
%%
ibuf = 1;
nex_prop = 0;
nex_accept  = 0;
poolobj = gcp('nocreate');
delete(poolobj);
parpool(NT);
%parfor ic = 1:4;
%
%  objtmp = obj(ic);
%  [obj(ic)] = explore_mhbd(obj(ic),obj(ic).beta,param);
%  [obj(ic)] = explore_mh(obj(ic),obj(ic).beta,param);
%
%end;
%stop

idbacc = zeros(NT);
tic;
for imcmc = 1:NMCMC;
  %for ithin = 1:NTHIN;
    parfor ic = 1:NT;
    for ithin = 1:NTHIN;
      if(kmin ~= kmax);[obj(ic)] = explore_mhbd(obj(ic),obj(ic).beta,param);end;
      [obj(ic)] = explore_mh(obj(ic),obj(ic).beta,param);
    end;
    end;
    %% Pass structure array to temp_exchange:
    if(iexchange == 1);
      [obj,nex_accept,nex_prop] = temp_exchange(obj,nex_accept,nex_prop,param);
    end;
  %end;
  %% Adapt step size based on acceptance history
  %% (does not properly keep track of b/d changes)
  if(mod(imcmc,NKEEP)==0);
    for ic = 1:NT;
      for kwhich = 1:obj(ic).k;
        for iwhich = 1:NPV;
          if(obj(ic).iaccept(kwhich,iwhich)/obj(ic).ipropose(kwhich,iwhich) < 0.2);
            obj(ic).step_sz(kwhich,iwhich) = obj(ic).step_sz(kwhich,iwhich)/1.2;
          end;
          if(obj(ic).iaccept(kwhich,iwhich)/obj(ic).ipropose(kwhich,iwhich) > 0.4);
            obj(ic).step_sz(kwhich,iwhich) = obj(ic).step_sz(kwhich,iwhich)*1.2;
          end;
        end;
      end;
    end;
  end;
  for ic = 1:NT1;
    buf(ibuf) = obj(ic);
    ibuf = ibuf + 1;
  end;
  if(mod(ibuf,NKEEP)==0);
    for ic = 1:NT;
      ibdacc(ic) = obj(ic).iacceptbd;
    end;
    tmcmc = toc;
    disp([imcmc,1.,obj(1).k,obj(1).BIC,obj(1).iacceptbd/obj(1).iproposebd,...
          nex_accept/nex_prop,tmcmc]);
    disp(ibdacc)
    save('sample.mat','buf','-append');
    tic;
  end;
end;

return;  %%% End of main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Metropolis-Hastings birth/death step
%%
function [obj]=explore_mhbd(obj,beta,param);

%% Propose b/d with equal probability:
obj.iproposebd = obj.iproposebd + 1;
i_bd = 0;
ran_unik = rand;
if(obj.k == param.kmax);      %% If k == kmax, no birth allowed
  if(ran_unik>=0.3333);i_bd = 2;end;
elseif(obj.k == param.kmin);  %% If k == kmin, no death allowed
  if(ran_unik>=0.3333);i_bd = 1;end;
else
  if(ran_unik<=0.3333);i_bd = 1;end;
  if(ran_unik>=0.6666);i_bd = 2;end;
end
%% Important to keep obj unaltered here to ensure proper rejection.
if(i_bd == 0);objnew=obj;end;
if(i_bd == 1);[objnew]=birth(obj,param);end;
if(i_bd == 2);[objnew]=death(obj,param);end;

%% If a b or d was performed, eval MH, otherwise skip to normal MH step.
if(obj.k ~= objnew.k);
  [objnew.BIC,objnew.sl1] = get_bic(objnew.voro(1:objnew.k,:),objnew.obs,param);
  %% Check MH probability (b/d is from prior, so only likelihood ratio needed)
  logLratio = (objnew.BIC - obj.BIC)*beta;
  if(rand < exp(logLratio));
    %% accept
    %disp('accept')
    objnew.iacceptbd = objnew.iacceptbd + 1;
    obj = objnew;
  end;
end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Metropolis-Hastings step without dimension change.
%%
function [obj]=explore_mh(obj,beta,param);

for kwhich = 1:obj.k;
for iwhich = 1:param.NPV;
  obj.ipropose(kwhich,iwhich) = obj.ipropose(kwhich,iwhich) + 1;
  %% Important to have obj unchanged here to ensure proper rejection.
  [objnew] = perturb(obj,iwhich,kwhich,param);
  [ioutside]=checkbounds2(objnew,kwhich,param);
   
  if(ioutside == 0);
    [objnew.BIC,objnew.sl1] = get_bic(objnew.voro(1:objnew.k,:),objnew.obs,param);
    %% Check MH probability
    logLratio = (objnew.BIC - obj.BIC)*beta;
    if(rand < exp(logLratio));
      %% accept
      %disp('accept')
      objnew.iaccept(kwhich,iwhich) = objnew.iaccept(kwhich,iwhich) + 1;
      obj = objnew;
    end;
  else;
    %% reject
    ioutside = 0;
  end;
end;
end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Perturb voro value (Gaussian centered on current)
%%
function [objnew]=perturb(obj,iwhich,kwhich,param);
objnew = obj;
objnew.voro(kwhich,iwhich) = obj.voro(kwhich,iwhich) + obj.step_sz(kwhich,iwhich)*randn();
%if(obj.voro(kwhich,iwhich) < minlim(iwhich));ioutside = 1;end;
%if(obj.voro(kwhich,iwhich) > maxlim(iwhich));ioutside = 1;end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Birth node (from prior)
%%
function [objnew]=birth(obj,param);

objnew = obj;
objnew.k = obj.k + 1;
objnew.voro = zeros(objnew.k,param.NPV);
objnew.voro(1:obj.k,:) = obj.voro;
ioutside = 1;
while ioutside == 1;
  objnew.voro(objnew.k,1) = rand*param.maxlim(1);
  objnew.voro(objnew.k,2) = rand*param.maxlim(2);
  [tmp,idxr]=min(abs(round(objnew.voro(objnew.k,1)-param.x_ev)));
  [tmp,idxz]=min(abs(round(objnew.voro(objnew.k,2)-param.z_ev)));
  if(idxr <= param.idx_exl(idxz) | idxr >= param.idx_exr(idxz));
    ioutside = 1;
  else;
    ioutside = 0;
  end;
end;
%% Sample new parameters from prior (results in |J|=1 and 
%% MH acceptance simplifies to likelihood ratio).
for i=1:param.NTW;
  objnew.voro(objnew.k,2+i) = param.minlim(2+i) + rand*param.maxpert(2+i);
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Death node (from prior)
%%
function [objnew]=death(obj,param);

objnew = obj;
objnew.k = obj.k - 1;
idxrand = randperm(obj.k);
kdel = idxrand(1);
objnew.voro(kdel,:) = [];
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Exchange move (parallel tempering)
%%
function [obj,nex_accept,nex_prop]=temp_exchange(obj,nex_accept,nex_prop,param);

%% Propose 2*NT swaps in chain cloud (arbitrary choice)...
for ic = 1:2*length(obj);
  nex_prop = nex_prop + 1;  %% Proposal counter
  idxtmp = randperm(length(obj));
  ic1 = idxtmp(1);
  ic2 = idxtmp(2);
  betaratio = obj(ic2).beta-obj(ic1).beta;
  logratio  = betaratio*(obj(ic1).BIC-obj(ic2).BIC);
  if(rand <= exp(logratio));
    %% ACCEPT SWAP
    objtmp1    = obj(ic1);
    objtmp2    = obj(ic2);
    obj(ic1)   = objtmp2;
    obj(ic2)   = objtmp1;
    %% Temperature does not swap
    obj(ic1).beta = objtmp1.beta;
    obj(ic2).beta = objtmp2.beta;
    %% BD acceptance counters do not swap
    obj(ic1).iacceptbd  = objtmp1.iacceptbd;
    obj(ic2).iacceptbd  = objtmp2.iacceptbd;
    obj(ic1).iproposebd = objtmp1.iproposebd;
    obj(ic2).iproposebd = objtmp2.iproposebd;
    %% Voro acceptance counters do not swap
    obj(ic1).iaccept  = objtmp1.iaccept;
    obj(ic2).iaccept  = objtmp2.iaccept;
    obj(ic1).ipropose = objtmp1.ipropose;
    obj(ic2).ipropose = objtmp2.ipropose;
    %% Step sizes do not swap
    obj(ic1).step_sz = objtmp1.step_sz;
    obj(ic2).step_sz = objtmp2.step_sz;
    
    nex_accept    = nex_accept + 1;  %% Acceptance counter
  end;
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Perturb voro value (Gaussian centered on current)
%%
function [ioutside]=checkbounds(obj,ic,param);

ioutside = 0;
for kwhich = 1:obj.k;
  [tmp,idxr]=min(abs(round(obj.voro(kwhich,1)-param.x_ev)));
  [tmp,idxz]=min(abs(round(obj.voro(kwhich,2)-param.z_ev)));
  if(idxr <= param.idx_exl(idxz) | idxr >= param.idx_exr(idxz));
    ioutside = 1;
    if(ic == 1);
      disp(['bounds out for k = ',num2str([kwhich,obj.voro(kwhich,1:2)])]);
    end;
  end;
  for iwhich = 1:param.NPV;
    if(obj.voro(kwhich,iwhich) <= param.minlim(iwhich));
      ioutside = 1;
      if(ic == 1);
        disp(['min bounds out for ',num2str([kwhich,iwhich,obj.voro(kwhich,:)])]);
      end;
    end;
    if(obj.voro(kwhich,iwhich) >= param.maxlim(iwhich));
      ioutside = 1;
      if(ic == 1);
        disp(['max bounds out for ',num2str([kwhich,iwhich,obj.voro(kwhich,:)])]);
      end;
    end;
  end;
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Perturb voro value (Gaussian centered on current)
%%
function [ioutside]=checkbounds2(obj,kwhich,param);

ioutside = 0;
[tmp,idxr]=min(abs(round(obj.voro(kwhich,1)-param.x_ev)));
[tmp,idxz]=min(abs(round(obj.voro(kwhich,2)-param.z_ev)));
if(idxr <= param.idx_exl(idxz) | idxr >= param.idx_exr(idxz));
  ioutside = 1;
  %disp(['bounds out',num2str([kwhich,obj.voro(kwhich,1:2)])]);
end;
for iwhich = 1:param.NPV;
  if(obj.voro(kwhich,iwhich) <= param.minlim(iwhich));ioutside = 1;end;
  if(obj.voro(kwhich,iwhich) >= param.maxlim(iwhich));ioutside = 1;end;
end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute BIC value
%%
function [BIC,sl1]=get_bic(voro,obs,param);

k = size(voro,1);
MPAR = k*param.NPV;

[sl1]=voro_interp(voro,param);
pred = sl1(:);  %% This reshapes into single row vector
res = (obs-pred).*param.rex;
logL = -param.NDAT/2.* log(sum(res.*res));
sdpar = sqrt(1./param.NDAT*sum(res.*res));
if(param.ibic == 1);
  BIC = 2.*logL - MPAR*log(param.NDAT);
else;
  BIC = logL;
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Voronoi nearest neighbour interp (Euklidean norm)
%%
function [sl1]=voro_interp(voro,param);

k=size(voro,1);
d = zeros(k,1);
dx = zeros(k,1);
dz = zeros(k,1);
sl1 = zeros(param.NDEP,param.NRAN,param.NTW);

x  = voro(1:k,1);
z  = voro(1:k,2);

%% Plot MAP
for iz = 1:param.NDEP;
  dz = param.z_ev(iz)-z;
  dzsq = dz.*dz;
  for ir = 1:param.NRAN;
    if(param.mat_excl(iz,ir) == 1);
      dx = param.x_ev(ir)-x;
      d  = sqrt(dx.*dx + dzsq);
      [dmin,iv] = min(d);
      sl1(iz,ir,1:param.NTW) = voro(iv,3:2+param.NTW);
    end;
  end;
end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Voronoi nearest neighbour interp (Euklidean norm)
%%
function [sl1]=voro_interp_array(voro,param);

k=size(voro,1);
d = zeros(k,1);
dx = zeros(k,1);
dz = zeros(k,1);

x  = voro(1:k,1);
z  = voro(1:k,2);

%% this is k by NDEP
dzsq = (repmat(z,[1,NDEP])-repmat(param.z_ev,[k,1])).^2;
%% this is k by NRAN
dxsq = (repmat(x,[1,NRAN])-repmat(param.x_ev,[k,1])).^2;
%% Plot MAP
for iz = 1:NDEP;
  %dz = z_ev(iz)-z;
  %dzsq = dz.*dz;
  for ir = 1:NRAN;
    if(param.mat_excl(iz,ir) == 1);
      %dx = x_ev(ir)-x;
      d  = sqrt(dxsq(:,ir) + dzsq(:,iz));
      [dmin,iv] = min(d);
      sl1(iz,ir,1:param.NTW) = voro(iv,3:2+param.NTW);
    end;
  end;
end;

return;
%%
%%
%%EOF
