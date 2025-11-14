function [outfile,F,converg_cov,converg_marg] =...
          convert_parfile_hol(tab);

n1 = 19;
outfile    = tab(1,16:end);
inoise      = str2num(tab(2,16:end));
F(1).ifreq  = [str2num(tab(3,16:end)):str2num(tab(4,16:end)):...
               str2num(tab(5,16:end))];
F(1).iang   = [str2num(tab(6,16:end)):str2num(tab(7,16:end)):...
               str2num(tab(8,16:end))];
F(1).astart = str2num(tab(9,16:end));
F(1).aend   = str2num(tab(10,16:end));
F(1).nmod   = str2num(tab(11,16:end));
F(1).sd     = str2num(tab(12,16:end));
converg_cov = str2num(tab(13,16:end));
converg_marg= str2num(tab(14,16:end));
F(1).ave    = str2num(tab(15,16:end));
F(1).cw     = str2num(tab(16,16:end));
F(1).rw     = str2num(tab(17,16:end));
n2 = n1 + F(1).nmod + 1;
n3 = n2 + F(1).nmod + 1;
% model parameters: h, c1,rho1,alpha1, c2,rho2,alpha2
for i = 1:F(1).nmod
  msim(i)   = str2num(tab(n1+i-1,16:end));
end
F(1).msim = msim';
for i = 1:F(1).nmod
  minlim(i) = str2num(tab(n2+i-1,16:end));
end
F(1).minlim = minlim';
for i = 1:F(1).nmod
  maxlim(i) = str2num(tab(n3+i-1,16:end));
end
F(1).maxlim = maxlim';


return;
