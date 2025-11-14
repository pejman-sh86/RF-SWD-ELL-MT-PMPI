function XX= test2_arma;

%%
%% Load delta residual data (power law residuals)
%%
%NDATA = 5;
%NDATB = 30;
NDATA = 1;
NDATB = 44;
idiff = 0;

dat=load('fort.3');
dat = dat(2:end,:);
dat = dat(NDATA:NDATB,:);
rep=load('power_replica.dat');
%rep=load('test_rep.txt');
rep = rep(NDATA:NDATB,:);
res = dat(:,2)-rep(:,2);

if(idiff == 1)
   resnew = res;
   resnew(2:NDATB) = mean(diff(res))+res((1:end-1));
   figure(1);
   plot(res,'-b');hold on;
   plot(resnew,'--r');
   plot(res-resnew,'--k');
   res = res-resnew;
end

x=res;
N=length(x);
order = 14;

fi = zeros(order,order);
for k = order:order;
   clear XX a alpha est_x

   a = aryule(x,k);

   aa = aryule(x,k)
   bb = aryule(flipud(x),k)

   alpha = -a(2:end);
   est_x = filter([0 -a(2:end)],1,x);    % Estimated signal

   salpha(k) = sum(alpha);
   saalpha(k) = sum(abs(alpha));
   fi(k,1:k) = alpha;

%   XX(1)=x(1);
   XX(1)=0;

   x2 = flipud(x);
%   XX2(1)=x2(1);
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
   XX = XX';
   XX2 = flipud(XX2');
   size(XX2)
   XX3 = (XX+XX2)/2.;
%   XX3 = XX2;

   figure(2);box on;
   plot(x,'-b');hold on;
%   plot(est_x,'-r');
   plot(XX,':r');
   plot(XX2,':k');
   plot(XX3,'-k');
   plot([0 N],[0 0],':k');
   set(gca,'XLim',[0 N]);
   title('Original Signal vs. LPC Estimate');
   xlabel('Datum');ylabel('Residual error (m/s)');
   legend('Original Signal','ARYULE Forward Estimate','ARYULE Back Estimate','ARYULE Forw-Back Estimate')

   figure(3);hold on;box on;
   plot(x-XX,':b');plot([0 N],[0 0],':k');
   plot(x-XX2,':r');
   plot(x-XX3,'k');
   set(gca,'XLim',[0 N]);
   xlabel('Datum');ylabel('Residual error (m/s)');
   legend('ARYULE Forward Estimate','','ARYULE Back Estimate','ARYULE Forw-Back Estimate')

   figure(4);hold on;box on;
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
   legend('raw','ARYULE Forward Estimate','ARYULE Forw-Back Estimate')

   figure(5);hold on;box on;
   [acs1,lags1] = xcorr(x,'coeff');   % ACS of prediction error
   [acs2,lags2] = xcorr(XX,'coeff');   % ACS of forw prediction error
   [acs3,lags3] = xcorr(XX3,'coeff');   % ACS of forw-back prediction error
   plot(lags1,acs1,'-b');
   plot(lags2,acs2,'-r');
   plot(lags3,acs3,'-k');
   plot([-N N],[0 0],':k');
   set(gca,'XLim',[-N N]);
   rms(k) = sqrt(mean((x-XX).^2));
   xlabel('Lag');ylabel('ACF');
   legend('raw','ARYULE Forward Estimate','ARYULE Forw-Back Estimate')

   figure(6);hold on;box on;
   plot(fi(k,1:k));xlabel('Coeff IDX');ylabel('Alpha');
   plot([0 order],[0 0],':k');

end;
figure(6);box on;
plot(rms);xlabel('Order');ylabel('RMS error');
figure(7);box on;
plot(salpha);xlabel('Order');ylabel('sum(alpha)');
figure(8);box on;
plot(saalpha);xlabel('Order');ylabel('sum(abs(alpha))');

return;
