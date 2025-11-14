
filebase='tmp_'
iar = 1;
repfile   = strcat(filebase,'replica.dat');
repscfile = strcat(filebase,'replicascaled.dat');
resfile   = strcat(filebase,'residuals.dat');
resrawfile= strcat(filebase,'residualsraw.dat');
resarfile = strcat(filebase,'residualsar.dat');
obsfile   = strcat(filebase,'observed.dat');
Afile     = strcat(filebase,'repscaleA.dat');
rep   = load(repfile);
repsc = load(repfile);
res   = load(resfile);
resraw= load(resrawfile);
resar = load(resarfile);
obs   = load(obsfile);
A     = load(Afile);

NBANDS = 11;
NHYD   = 48;

NREP = length(rep)/NBANDS/2;

zobs = complex(obs(1:NBANDS,:),obs(NBANDS+1:2*NBANDS,:));
zobs = zobs.';

for ifr=1:NBANDS;
   obsmxr(ifr) = max(abs(real(zobs(:,ifr))));
   obsmxi(ifr) = max(abs(imag(zobs(:,ifr))));
end;

for irep=1:NREP;

   ireal = (irep-1)*NBANDS*2+1;
   iimag = (irep-1)*NBANDS*2+1+NBANDS;
   zrep(irep,:,:) = complex(rep(ireal:ireal+NBANDS-1,:),...
                    rep(iimag:iimag+NBANDS-1,:));
   zrepsc(irep,:,:) = complex(repsc(ireal:ireal+NBANDS-1,:),...
                      repsc(iimag:iimag+NBANDS-1,:));
   zres(irep,:,:) = complex(res(ireal:ireal+NBANDS-1,:),...
                    res(iimag:iimag+NBANDS-1,:));
   zresraw(irep,:,:) = complex(resraw(ireal:ireal+NBANDS-1,:),...
                       resraw(iimag:iimag+NBANDS-1,:));
   zA(irep,:) = complex(A((irep-1)*2+1,:),...
                A((irep-1)*2+2,:));
   if(iar == 1)
      zar(irep,:,:) = complex(resar(ireal:ireal+NBANDS-1,:),...
                      resar(iimag:iimag+NBANDS-1,:));
   end;

end;

for irep=1:NREP;
for ifr=1:NBANDS;
   ztmp = squeeze(zrep(irep,ifr,:));
   zrepscaled(irep,ifr,:) = ztmp*((ztmp'*zobs(:,ifr))/...
                      (ztmp'*ztmp));
   zAA(irep,ifr) = ((ztmp'*zobs(:,ifr))/(ztmp'*ztmp));
%  zrepscaled=zrep*((zrep'*zobs)/(zrep'*zrep));

end;
end;

matfile = strcat(filebase,'repres.mat');
save(matfile,'zobs','zrep','zrepsc','zres','zresraw','zar','obsmxr','obsmxi','zrepscaled','NREP','zA','zAA');

