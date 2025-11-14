
filebase = 'HON';
include_rv = 0;
include_swd = 0;
include_ell = 0;
include_mt = 1;
parfile = strcat(filebase,'_parameter.dat');
icovsfile = strcat(filebase,'_icovmat.dat');
Cdifile = strcat(filebase,'_Cdi.dat');
Cdiswdfile = strcat(filebase,'_CdiSWD.dat');
Cdiellfile = strcat(filebase,'_CdiELL.dat');
Cdimtfile = strcat(filebase,'_CdiMT.dat');

%%
[IMAP,ICOV,iar,i_varpar,irv,itr,iswd, iell, imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL, NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVs,dVpVs, sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile(parfile);
NRF = ntr;

%%
covs = importdata(icovsfile); 

%%
if include_rv == 0; NRF = 0; end % be sure when include_rv == -1, NRF ~= 0
if include_swd == 0; NMODE = 0; end % be sure when include_swd == 1, NMODE ~= 0
if include_ell == 0; NMODE_ELL = 0; end % be sure when include_ell == 1, NMODE_ELL ~= 0

%%
if include_rv == -1
    Cdi = covs(1:NTIME, 1:NTIME);
    save(Cdifile, 'Cdi', '-ascii')
end

if include_swd == 1
    Cdiswd = covs(NRF*NTIME+1:NRF*NTIME+NDAT_SWD, 1:NDAT_SWD);
    save(Cdiswdfile, 'Cdiswd', '-ascii')
end

if include_ell == 1
    Cdiell = covs(NRF*NTIME+NMODE*NDAT_SWD+1:NRF*NTIME+NMODE*NDAT_SWD+NDAT_ELL, 1:NDAT_ELL);
    save(Cdiellfile, 'Cdiell', '-ascii')
end

if include_mt == 1
    Cdimt = covs(NRF*NTIME+NMODE*NDAT_SWD+NMODE_ELL*NDAT_ELL+1:NRF*NTIME+NMODE*NDAT_SWD+NMODE_ELL*NDAT_ELL+NDAT_MT, 1:2*NDAT_MT);
    save(Cdimtfile, 'Cdimt', '-ascii')
end
