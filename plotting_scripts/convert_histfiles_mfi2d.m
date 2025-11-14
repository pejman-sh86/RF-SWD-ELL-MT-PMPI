function []=convert_histfiles_mfi2d();

NHIST = 150;
NDEP  = 1001;
NRAN  = 401;
c_ens = load('rjmh_pecan_ensavec.dat');
r_ens = load('rjmh_pecan_ensaver.dat');
a_ens = load('rjmh_pecan_ensavea.dat');

disp('loading c');
c_hst = load('rjmh_pecan_histc.dat');
disp('loading r');
r_hst = load('rjmh_pecan_histr.dat');
disp('loading a');
a_hst = load('rjmh_pecan_hista.dat');
bins  = load('bins.dat');

for i=1:NHIST;
  c_hst2(1:NDEP,1:NRAN,i) = c_hst((i-1)*NDEP+1:i*NDEP,:);
  r_hst2(1:NDEP,1:NRAN,i) = r_hst((i-1)*NDEP+1:i*NDEP,:);
  a_hst2(1:NDEP,1:NRAN,i) = a_hst((i-1)*NDEP+1:i*NDEP,:);
end;

clear c_hst r_hst a_hst i;
save hists.mat;

return;
