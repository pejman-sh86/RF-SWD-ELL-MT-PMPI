%
% Simulate Harrisons data unsing parameters and model from
% parameters.dat
%
function [] = simulate_har();

tab = readtextfile('measured_parameters_3_750-1265.dat');
[outfile,F,converg_cov,converg_marg] =...
convert_parfile(tab);
%extension = strcat('rep_data',num2str(nlay),num2str(f));
simfile = strrep(outfile,'sample.dat','rep_data_3_750-1265.mat');
nang = length(F(1).ang);
%m = F(1).msim;
%m = [0.63914583      1622.5336      1.6601592 0.78298314      2.9083414      1580.6068      1.0967968     0.37488555 1850      2.0664243     0.76113981]';
%m = [0.329139      1581.3119      1.4393652     0.49350746     0.55101746 1527.9471      1.2320978   0.43199968      1716.2655      1.7616306 0.54225654]';
m = [0.29030184       1579.469      1.3696368 0.99992177     0.58739455      1564.7047      1.0284575      0.3533903 1730.8198      1.5588327     0.68915472]';

[Bft] = Beams([.5:.5:16],F(1).freq,1500,1024);

for ifreq = 1:length(F(1).freq)
  F(ifreq).dat = zeros(1,nang);
end

[E,B] = forward_har_new(m,F,Bft);

save(simfile,'B');

return;
