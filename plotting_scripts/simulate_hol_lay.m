%
% Simulate Holland data unsing parameters and model from
% parameters.dat
% Frequency average for layered models
%
function [] = simulate_hol_lay();

%load ../sim_A_1dB_1.mat;
load oastl_bench.mat;
simfile = ('sph_40ave.mat');
inoise = 0;
ispher = 1;
sd = 1;

%m = [ 5.  1540.0 1.4 0.1 ...
%      5.  1520.0 1.5 0.1 ...
%          1560.0 1.6  0.1]';
%m = [1.5 1511 1.55 0.10 ...
%         1511 1.55 0.10]';
%m = [1550.00 1.55  0.1]';
%m = [2.0  1520.00 1.400  0.05 ...
%     18.0 1540.00 1.550  0.1 ...
%     80.0 1600.00 1.600 0.15 ...
%          2400.00 2.200 0.01]';
%m = [1540 1.40 0.05]';
%m = [1.5 1540 1.40 0.05 ...
%         1520 1.55 0.10]';
%m = [1.5 1540 1.40 0.05 ...
%     3.0 1520 1.55 0.10 ... 
%     1.0 1560 1.75 0.15 ...
%         1600 1.90 0.10]';
%m = [1.5 1540 1.40 0.05 ... 
%     3.0 1520 1.55 0.10 ...
%     1.0 1560 1.75 0.15 ... 
%     7.0 1600 1.90 0.10 ...
%     1.0 1700 2.05 0.10 ...
%         1600 2.10 0.10]';
m = [1.5 1540 1.40 0.05 ...
     3.0 1520 1.55 0.1 ...
     1.0 1560 1.75 0.15 ...
     7.0 1600 1.90 0.1 ...
     1.0 1700 2.05 0.1 ...
     4.0 1600 2.10 0.1 ...
         1650 2.20 0.01]';

%load sim_A_1dB_4_map;
%m = xmap;

F(1).nmod = length(m);
%F(1).freq = [50 63 80 100 125 160 200 250 315];
%F(1).freq = [400 500 630 800 1000];
%F(1).freq = [400 500 630 800 1000 1250 1600 2000];
F(1).freq = [400 600 900 1350 2000];
%F(1).freq = [2000];
%save bla F;
for ifreq = 1:length(F(1).freq)
  nang = length(F(1).ang);
  F(ifreq).dat = zeros(1,nang)';
end

if ispher == 0
    [E,Fm] = forward_hol_lay(m,F);
else
    [E,Fm] = forward_hol_lay_spher(m,F);
end

for ifreq = 1:length(F(1).freq)

  F(ifreq).dat = Fm(ifreq).dat;
  if ispher == 1
    F(ifreq).plane = Fm(ifreq).plane;
  end

  if inoise == 1
%    randn('state',sum(100*clock))
%    s = randn('state')
    noise = randn(length(F(ifreq).dat),1);
%    L = chol(F(ifreq).csave);
%    cnoise = L'*noise;
    cnoise = sd*noise;

    F(ifreq).datnn = F(ifreq).dat;
    F(ifreq).dat = cnoise + Fm(ifreq).dat;
  
    F(ifreq).sd = std(cnoise);
  end

end
%for ifreq = 1:length(F(1).ifreq)
%  B(ifreq).ang = F(ifreq).ang;
%end

B = F(1:length(F(1).freq));
clear F;
F=B;

save(simfile,'F');
%save(simfile,'B');

return;
