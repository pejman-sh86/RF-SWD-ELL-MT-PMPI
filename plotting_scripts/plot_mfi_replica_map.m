function [] = plot_mfi_replica();

inorm1  = 1; %% Normalize profile marginals line by line
inorm2  = 1; %% Normalize profile marginals line by line
iar = 1;

NBANDS=  3;  %% No. freq bands
NSR   =  3;  %% No. source ranges
NHYD  = 32;  %% No. hydrophones
NDAT = NBANDS*NSR;

rep   = load('rjmh_pecan_replica.dat');
obs   = load('rjmh_pecan_observed.dat');

z     = rep(1,:);
rep   = rep(2:end,:);
zobs  = complex(obs(1:2:end,:),obs(2:2:end,:));
zrep  = complex(rep(1:2:end,:),rep(2:2:end,:));
zobs  = zobs.';
zrep  = zrep.';

for ifr=1:NDAT;
   obsmxr(ifr) = max(abs(real(zobs(:,ifr))));
   obsmxi(ifr) = max(abs(imag(zobs(:,ifr))));
end;

for ifr=1:NDAT;
   zrepscaled(:,ifr) = zrep(:,ifr)*((zrep(:,ifr)'*zobs(:,ifr))/...
                      (zrep(:,ifr)'*zrep(:,ifr)));
end;

%%
%%  Data plots
%%
   fig1=figure(1);hold on;box on;
   title('Data misfit real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NDAT;
      subplot(4,3,ifr);hold on;box on;
      plot([1:NHYD],real(zobs(:,ifr))/obsmxr(ifr),'-k');
      plot([1:NHYD],real(zrepscaled(:,ifr))/obsmxr(ifr),'--b');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
   end;

   fig2=figure(2);hold on;box on;
   title('Data misfit imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NDAT;
      subplot(4,3,ifr);hold on;box on;
      plot([1:NHYD],imag(zobs(:,ifr))/obsmxi(ifr),'-k');
      plot([1:NHYD],imag(zrepscaled(:,ifr))/obsmxi(ifr),'--b');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
   end;

   fig3=figure(3);hold on;box on;
   title('Data misfit real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NDAT;
      subplot(4,3 ,ifr);hold on;box on;
      plot([1:NHYD],(real(zobs(:,ifr))-real(zrepscaled(:,ifr)))/obsmxr(ifr),'-k');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
   end;


return;
