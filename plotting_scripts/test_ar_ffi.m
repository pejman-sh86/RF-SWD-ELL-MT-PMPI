function XX= test2_arma(filebase);

parfile  = strcat(filebase,'_parameter.dat');
repfile  = strcat(filebase,'_replica.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF]=ffi_read_parfile(parfile);

NSTN  = h5readatt(datfile,'/Observed_data','N_sta');
deltt = h5readatt(datfile,'/Observed_data','Sample_rate');
hyp_loc = h5read(datfile,'/Observed_data/Hypo_Loc');
dhyp  = h5readatt(datfile,'/Sensitivity_kernel','Hyp_interval');
NRAN  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP  = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
Rmx   = h5readatt(datfile,'/Sensitivity_kernel','max_x');
Zmx   = h5readatt(datfile,'/Sensitivity_kernel','max_y');
Vrmin = h5readatt(datfile,'/Sensitivity_kernel','V_r_min');
Vrmax = h5readatt(datfile,'/Sensitivity_kernel','V_r_max');

dat = h5read(datfile,'/Observed_data/displacements');
dat = dat';
NTSMP = h5read(datfile,'/Observed_data/Ntraces');
NDAT = length(dat)
rep=dlmread(repfile);

res = dat-rep(3,:);

for istn = 1:NSTN;

  iend = sum(NTSMP(1:istn));
  istart = iend-NTSMP(istn)+1;

  x=res(istart:iend);
  N=length(x);
  order = 2


  fi = zeros(order,order);
  for k = order:order;
    clear XX a alpha est_x

    a = aryule(x,k)

    alpha = -a(2:end);
    est_x = filter([0 -a(2:end)],1,x);    % Estimated signal

    salpha(k) = sum(alpha);
    saalpha(k) = sum(abs(alpha));
    fi(k,1:k) = alpha;

   XX(1)=0;
   for i=2:N;
      XX(i) = 0;
      if(k>=i)
         for j=1:i-1;
            XX(i) = XX(i) + alpha(j) * x(i-j);
         end;
      else
         for j=1:k;
            XX(i) = XX(i) + alpha(j) * x(i-j);
         end;
      end;
   end;

   figure(2);box on;
   plot(x,'-b');hold on;
%   plot(est_x,'-r');
   plot(XX,':r');
   plot([0 N],[0 0],':k');
   set(gca,'XLim',[0 N]);
   title('Original Signal vs. LPC Estimate');
   xlabel('Datum');ylabel('Residual error (m/s)');
   legend('Original Signal','ARYULE Forward Estimate')

   figure(3);hold on;box on;
   plot(x-XX,':b');plot([0 N],[0 0],':k');
   set(gca,'XLim',[0 N]);
   xlabel('Datum');ylabel('Residual error (m/s)');
   legend('ARYULE Forward Estimate')

   figure(4);hold on;box on;
   [acs1,lags1] = xcorr(x,'coeff');   % ACS of prediction error
   [acs2,lags2] = xcorr(x-XX,'coeff');   % ACS of forw prediction error
   plot(lags1,acs1,'-b');
   plot(lags2,acs2,'-r');
   plot([-N N],[0 0],':k');
   set(gca,'XLim',[-N N]);
   rms(k) = sqrt(mean((x-XX).^2));
   xlabel('Lag');ylabel('ACF');
   legend('raw','ARYULE Forward Estimate')

   figure(5);hold on;box on;
   [acs1,lags1] = xcorr(x,'coeff');   % ACS of prediction error
   [acs2,lags2] = xcorr(XX,'coeff');   % ACS of forw prediction error
   plot(lags1,acs1,'-b');
   plot(lags2,acs2,'-r');
   plot([-N N],[0 0],':k');
   set(gca,'XLim',[-N N]);
   rms(k) = sqrt(mean((x-XX).^2));
   xlabel('Lag');ylabel('ACF');
   legend('raw','ARYULE Forward Estimate')

   figure(6);hold on;box on;
   plot(fi(k,1:k));xlabel('Coeff IDX');ylabel('Alpha');
   plot([0 order],[0 0],':k');

end;
end;
figure(6);box on;
plot(rms);xlabel('Order');ylabel('RMS error');
figure(7);box on;
plot(salpha);xlabel('Order');ylabel('sum(alpha)');
figure(8);box on;
plot(saalpha);xlabel('Order');ylabel('sum(abs(alpha))');

return;
