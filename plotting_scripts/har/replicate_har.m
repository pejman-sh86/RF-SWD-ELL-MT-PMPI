%
% Simulate Harrisons data unsing parameters and model from
% parameters.dat
%
function [B] = replicate_har(m,F);

%tab = readtextfile('measured_parameters_2_562.dat');
%[outfile,F,converg_cov,converg_marg] =...
%convert_parfile(tab);
%extension = strcat('rep_data',num2str(nlay),num2str(f));
%simfile = strrep(outfile,'sample.dat','rep_data_2_562.mat');
nang = length(F(1).ang);
%m = F(1).msim;
%m = []';

[Bft] = Beams([.5:.5:16],F(1).freq,1500,1024);

for ifreq = 1:length(F(1).freq)
  F(ifreq).dat = zeros(1,nang);
end

[E,B] = forward_har_new(m,F,Bft);

return;
