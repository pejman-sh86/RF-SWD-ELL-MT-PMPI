%
% Rearrange Holland data using
% parameters.dat
%
function [] = rearrange_hol();

ext1 = 'site20';
ext2 = '_b';
%plotfile1 = strrep(sample,'sample.dat','marginals.eps');

%dataload = strcat('xde4.mat');
%dataload = strcat('xbl_mn20_2_',ext1,'_b.mat');
%dataload = strcat('xbl_',ext1,'_noise.mat');
dataload = strcat('xbl_',ext1,'_ph3_neg.mat');
%dataload = strcat('xbl_',ext1,'_ph2_pos.mat');
datawrite = strcat('hol_',ext1,ext2,'.mat');
parfile = strcat('parameters_hol',ext2,'.dat')

A = load(dataload);
strrep(dataload,'.mat','');
tab = readtextfile(parfile);
varname = genvarname(strrep(dataload,'.mat',''))
A = getfield(A,varname);
%A = getfield(A,'xde4');

[outfile,F,converg_cov,converg_marg] =...
convert_parfile_hol(tab);
nfreq = length(F(1).ifreq);
nang = length(F(1).iang);

for ifreq = 1:nfreq
  jj = 1;
  for j = 1:nang
    if(isnan(A.bl(F(1).iang(j),ifreq-1+F(1).ifreq(1))) ~= 1)
      dat(jj) = A.bl(F(1).iang(j),ifreq-1+F(1).ifreq(1));
      ang(jj) = A.ang(F(1).iang(j));
      jj = jj+1;
    end
  end
  [F(ifreq).ang,IDX] = sort(ang);
  F(ifreq).dat = dat(IDX)';
end
F(1).freq = A.pref(F(1).ifreq,2);

F=rmfield(F,'ifreq');F=rmfield(F,'iang');
F=rmfield(F,'astart');F=rmfield(F,'aend');
F=rmfield(F,'ave');

save(datawrite,'F','outfile','converg_cov','converg_marg');

return;

