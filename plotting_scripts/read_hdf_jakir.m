
nfrac = 4.;

datfile  = 'tohoku_data.hdf5';

NSTN   = h5readatt(datfile,'/Observed_data','N_sta');
deltt  = h5read(datfile,'/Observed_data/Sample_rate');
hyp_loc= h5read(datfile,'/Observed_data/Hypo_Loc');
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx    = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx    = h5readatt(datfile,'/Sensitivity_kernel','max_y');

dat = h5read(datfile,'/Observed_data/displacements');
dat = dat';
NTSMP = h5read(datfile,'/Observed_data/Ntraces');
NDAT = length(dat);

A  = h5read(datfile,'/Sensitivity_kernel/A');
x  = h5read(datfile,'/Sensitivity_kernel/Linear_solution');
pred  = h5read(datfile,'/Sensitivity_kernel/synthetic_displacements');

figure();
plot(dat,'-b');hold on;
plot(pred,'--r');

res=dat-pred';
figure(); 
plot(res)

for istn=1:NSTN;
  NAVE = round(double(NTSMP(istn))/nfrac);
  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;
  dat2(istn).res = res(istart:iend);
  dat2(istn).axx = xcorr(dat2(istn).res);
  sd(istn) = std(dat2(istn).res);

  for idat=istart:iend;
      if(idat <= istart+NAVE);
        NAVE1 = istart;
      else
        NAVE1 = idat-NAVE;
      end;
      if(idat >= iend-NAVE-1);
        NAVE2 = iend;
      else
        NAVE2 = idat+NAVE;
      end;
      dat2(istn).sd_nonstat(idat) = sqrt(mean(res(NAVE1:NAVE2).^2.));
      %disp([idat,NAVE1,NAVE2,sd_nonstat(idat)]);
    end;

end;

%% Average noise level at each station:
figure();
plot(sd);

figure();
for istn=1:NSTN;
end;



