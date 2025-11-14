%
% Simulate Hollands data unsing parameters and model from
% parameters.dat
%
function [] = simulate_hol_grad();

dense  = 0;
ilin   = 0;
ispher = 0;
inoise = 0;

%load('syn_covbmf');
if inoise == 1
  load('cov_mat_est3');
  B=F;
end
load('hol_site20b');

if dense == 0
%    simfile = 'rep_1.mat';
%    simfile = 'site19_mapdata.mat';
    simfile = 'site20b_rep1.mat';
else
    simfile = 'site20b_rep1_dense.mat';
end

%m = F(1).msim;
m = [1. 1.2877492      1.6459754     0.30756672      1475.8274           1450     0.45016572];
%m = [1.5 1.3469251      1.4470517     0.71610216      1473.2794      1473.0577     0.29768832];
F(1).lay_thick = 0.02; % discretization for layer thickness - controls
                      % upper frequency limit
F(1).znorm=F(1).lay_thick/2:F(1).lay_thick:1-F(1).lay_thick/2;
F(1).sz=length(F(1).znorm);

%F(1).freq = F(1).freq(6);
%F(1).ang = F(6).ang;

for ifreq = 1:length(F(1).freq)
  nang = length(F(ifreq).ang);
  F(ifreq).dat = zeros(1,nang);
  if inoise == 1
    F(ifreq).csave = B(ifreq).csave;
  end
  if dense == 1
      F(ifreq).ang = [1:0.01:89];
  end
end

if inoise == 1;clear B;end;

if ilin == 0 
    if ispher == 0
      [E,Fm] = forward_hol_grad(m,F);
    else
      [E,Fm] = forward_hol_grad_spher(m,F);
    end 
else
    [E,Fm] = forward_hol_lingrad(m,F);
end

for ifreq = 1:length(F(1).freq)

  F(ifreq).dat = Fm(ifreq).dat;

  if inoise == 1
%    randn('state',sum(100*clock))
%    s = randn('state')
    noise = randn(length(F(ifreq).dat),1);
    L = chol(F(ifreq).csave);
    cnoise = L'*noise;

    F(ifreq).datnn = F(ifreq).dat;
    F(ifreq).dat = cnoise + Fm(ifreq).dat;
  
    F(ifreq).sd = std(cnoise);
  end

end

B = F(1:length(F(1).freq));
clear F;
F=B;

save(simfile,'F','outfile','converg_cov','converg_marg');

return;
