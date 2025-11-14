function XX= test2_arma;

%%
%% Load delta residual data (power law residuals)
%%
NBANDS = 11;
NHYD   = 48;
idiff = 0;

rep=load('rjmh_mfi_replica.dat');
obs=load('rjmh_mfi_observed.dat');

zrep = complex(rep(1:528),rep(529:end));
zobs = complex(obs(1:528),obs(529:end));

for ifr=1:NBANDS;

   zzrep(:,ifr) = zrep(((ifr-1)*NHYD)+1:ifr*NHYD);
   zzobs(:,ifr) = zobs(((ifr-1)*NHYD)+1:ifr*NHYD);

   zrepscaled(:,ifr) = zzrep(:,ifr)*((zzrep(:,ifr)'*zzobs(:,ifr))/...
                      (zzrep(:,ifr)'*zzrep(:,ifr)));

end;
resr = real(zrepscaled-zzobs);
resi = imag(zrepscaled-zzobs);

%res = resr;
res = resi;

for ifr = 1:NBANDS;
   if(idiff == 1)
      resnew = res(:,ifr);
      resnew(2:end) = mean(diff(res(:,ifr)))+res(1:end-1,ifr);
      figure(1);
      subplot(4,3,ifr);hold on;box on;
      plot(res(:,ifr),'-b');hold on;
      plot(resnew,'--r');
      plot(res(:,ifr)-resnew,'--k');
      res(:,ifr) = res(:,ifr)-resnew;
   end

   x=res(:,ifr);
   N=length(x);
   order = 3;

   fi = zeros(order,order);
   for k = order:order;
      clear XX a alpha est_x
      XX = zeros(NHYD,1);
      XX2= zeros(NHYD,1);
      XX3= zeros(NHYD,1);

      a = aryule(x,k);

      aa = aryule(x,k)
      bb = aryule(flipud(x),k)

      alpha = -a(2:end);
      est_x = filter([0 -a(2:end)],1,x);    % Estimated signal

      salpha(k) = sum(alpha);
      saalpha(k) = sum(abs(alpha));
      fi(k,1:k) = alpha;

%      XX(1)=x(1);
      XX(1)=0;

      x2 = flipud(x);
%      XX2(1)=x2(1);
      XX2(1)=0;

      for i=2:N;
         XX(i) = 0;
         XX2(i) = 0;
         if(k>=i)
            for j=1:i-1;
               XX(i) = XX(i) + alpha(j) * x(i-j);
               XX2(i) = XX2(i) + alpha(j) * x2(i-j);
            end;
         else
            for j=1:k;
               XX(i) = XX(i) + alpha(j) * x(i-j);
               XX2(i) = XX2(i) + alpha(j) * x2(i-j);
            end;
         end;
      end;
      XX2 = flipud(XX2);
      XX3 = (XX+XX2)/2.;
%      XX3 = XX2;

      figure(2);box on;
      subplot(4,3,ifr);hold on;box on;
      plot(x,'-b');hold on;
%      plot(est_x,'-r');
      plot(XX,':r');
      plot(XX2,':k');
      plot(XX3,'-k');
      plot([0 N],[0 0],':k');
      set(gca,'XLim',[0 N]);
      title('Original Signal vs. LPC Estimate');
      xlabel('Datum');ylabel('Residual error (m/s)');
      if(ifr == 1);legend('Original Signal','ARYULE Forward Estimate','ARYULE Back Estimate','ARYULE Forw-Back Estimate');end;

      figure(3);hold on;box on;
      subplot(4,3,ifr);hold on;box on;
      plot(x-XX,':b');plot([0 N],[0 0],':k');
      plot(x-XX2,':r');
      plot(x-XX3,'k');
      set(gca,'XLim',[0 N]);
      xlabel('Datum');ylabel('Residual error (m/s)');
      if(ifr == 1);legend('ARYULE Forward Estimate','','ARYULE Back Estimate','ARYULE Forw-Back Estimate');end;

      figure(4);hold on;box on;
      subplot(4,3,ifr);hold on;box on;
      [acs1,lags1] = xcorr(x,'coeff');   % ACS of prediction error
      [acs2,lags2] = xcorr(x-XX,'coeff');   % ACS of forw prediction error
      [acs3,lags3] = xcorr(x-XX3,'coeff');   % ACS of forw-back prediction error
      plot(lags1,acs1,'-b');
      plot(lags2,acs2,'-r');
      plot(lags3,acs3,'-k');
      plot([-N N],[0 0],':k');
      set(gca,'XLim',[-N N]);
      rms(k) = sqrt(mean((x-XX).^2));
      xlabel('Lag');ylabel('ACF');
      if(ifr == 1);legend('raw','ARYULE Forward Estimate','ARYULE Forw-Back Estimate');end;

      figure(5);hold on;box on;
      subplot(4,3,ifr);hold on;box on;
      plot(fi(k,1:k));xlabel('Coeff IDX');ylabel('Alpha');
      plot([0 order],[0 0],':k');

   end; % end order for
   figure(6);box on;
   subplot(4,3,ifr);hold on;box on;
   plot(rms);xlabel('Order');ylabel('RMS error');
   figure(7);box on;
   subplot(4,3,ifr);hold on;box on;
   plot(salpha);xlabel('Order');ylabel('sum(alpha)');
   figure(8);box on;
   subplot(4,3,ifr);hold on;box on;
   plot(saalpha);xlabel('Order');ylabel('sum(abs(alpha))');

end % end freq for
return;
