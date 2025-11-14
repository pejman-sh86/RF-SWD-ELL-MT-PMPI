
bands = [975. 1100. 1250. 2100. 2400. 2700.];
nf = length(bands);
pstgl  = 1;
pendgl = 170;
fileext  = '.txt';
fileext2  = '.png';

env = dlmread('track_environment_z.dat');
NPING = size(env,1);
NPL = 4;

ii = 1;
for iping = 1:NPING;
  for ilay = 1:env(iping,1);
    ipar = 3+(ilay-1)*NPL;
    cenv(ii) = env(iping,ipar);
    renv(ii) = env(iping,ipar+1);
    ii = ii+1;
  end;
end;
minlim = [1450. 1.2];
maxlim = [1750. 2.2];

rho = [minlim(2):(maxlim(2)-minlim(2))/1000.:maxlim(2)];

cl=(1.54-0.907*rho+0.3695*rho.^1.88)*1.5004*1000.;
ch=(1.6-0.907*rho+0.3695*rho.^2.01)*1.5014*1000.;

cl(find(cl<minlim(1))) = minlim(1);
ch(find(ch>maxlim(1))) = maxlim(1);

figure();
plot(rho,ch,'b');
hold on;
plot(rho,cl,'b');
plot(renv,cenv,'xr')