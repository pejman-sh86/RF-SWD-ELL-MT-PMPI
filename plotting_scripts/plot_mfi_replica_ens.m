function [] = plot_mfi_replica(represfile);
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

inorm1  = 1; %% Normalize profile marginals line by line
inorm2  = 1; %% Normalize profile marginals line by line
iar  = 0;
isave= 1;

NHYD   = 48;
NBANDS = 9;
bands = [400., 450., 500., 550.,...
         600., 650., 700., 750., 800.];
%NBANDS = 11;
%bands = [300., 350., 400., 450., 500., 550.,...
%         600., 650., 700., 750., 800.];

load(represfile);

filebase    = strrep(represfile,'repres.mat','');
plotfile1   = strcat(filebase,'datar.');
plotfile2   = strcat(filebase,'datai.');
plotfile3   = strcat(filebase,'datafr.');
plotfile4   = strcat(filebase,'datafi.');
plotfile5   = strcat(filebase,'totresr.');
plotfile6   = strcat(filebase,'totresi.');
plotfile7   = strcat(filebase,'rawresr.');
plotfile8   = strcat(filebase,'rawresi.');
plotfile9   = strcat(filebase,'arr.');
plotfile10  = strcat(filebase,'ari.');
plotfile11  = strcat(filebase,'axxr.');
plotfile12  = strcat(filebase,'axxi.');
plotfile13  = strcat(filebase,'axxr_raw.');
plotfile14  = strcat(filebase,'axxi_raw.');
plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

%NREP = length(rep)/NBANDS/2;
%
%zobs = complex(obs(1:NBANDS,:),obs(NBANDS+1:2*NBANDS,:));
%zobs = zobs.';
%
%for irep=1:NREP;
%
%   ireal = (irep-1)*NBANDS*2+1;
%   iimag = (irep-1)*NBANDS*2+1+NBANDS;
%   zrep(irep,:,:) = complex(rep(ireal:ireal+NBANDS-1,:),...
%                    rep(iimag:iimag+NBANDS-1,:));
%   zres(irep,:,:) = complex(res(ireal:ireal+NBANDS-1,:),...
%                    res(iimag:iimag+NBANDS-1,:));
%   if(iar == 1)
%      zar(irep,:,:) = complex(resar(ireal:ireal+NBANDS-1,:),...
%                      resar(iimag:iimag+NBANDS-1,:));
%   end;
%
%end;
%
%for irep=1:NREP;
%for ifr=1:NBANDS;
%   ztmp = squeeze(zrep(irep,ifr,:));
%   zrepscaled(irep,ifr,:) = ztmp*((ztmp'*zobs(:,ifr))/...
%                      (ztmp'*ztmp));
%%  zrepscaled=zrep*((zrep'*zobs)/(zrep'*zrep));
%
%end;
%end;

%%
%%  Data plots
%%
   nx = 3;
   ny = ceil(length(bands)/nx);
   xim = 0.01;
   yim = 0.05/ny;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

   fig1=figure(1);hold on;box on;
   title('Data misfit real part');
   hydplt = [0:NHYD+1];
   NV = 600;
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(real(zrepscaled(:,ifr,ihyd))/obsmxr(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for iz=1:NHYD
            Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([1:NHYD],real(zobs(:,ifr))/obsmxr(ifr),'-k');
      plot([1:NHYD],real(zobs(:,ifr))/obsmxr(ifr),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;
   end;

   fig2=figure(2);hold on;box on;
   title('Data misfit imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(imag(zrepscaled(:,ifr,ihyd))/obsmxi(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for iz=1:NHYD
            Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;

      plot([1:NHYD],imag(zobs(:,ifr))/obsmxi(ifr),'-k');
      plot([1:NHYD],imag(zobs(:,ifr))/obsmxi(ifr),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

%%
%%  Data plots (scaled in Fortran code)
%%
   fig3=figure(3);hold on;box on;
   title('Data misfit real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
%      vmin = -1.2;
%      vmax = 1.2;
      vmin = -obsmxr(ifr);
      vmax = obsmxr(ifr);
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
%         Nd(:,ihyd) = histc(real(zrepsc(:,ifr,ihyd))/obsmxr(ifr),vlim);
         Nd(:,ihyd) = histc(real(zrepsc(:,ifr,ihyd)),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for iz=1:NHYD
            Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;

%      plot([1:NHYD],real(zobs(:,ifr))/obsmxr(ifr),'-k');
%      plot([1:NHYD],real(zobs(:,ifr))/obsmxr(ifr),'--w');
      plot([1:NHYD],real(zobs(:,ifr)),'-k');
      plot([1:NHYD],real(zobs(:,ifr)),'--w');
%      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      set(gca,'XLim',[0 49])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

   fig4=figure(4);hold on;box on;
   title('Data misfit imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(imag(zrepsc(:,ifr,ihyd))/obsmxi(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for iz=1:NHYD
            Nd(:,iz) = Nd(:,iz)/max(Nd(:,iz));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;

      plot([1:NHYD],imag(zobs(:,ifr))/obsmxi(ifr),'-k');
      plot([1:NHYD],imag(zobs(:,ifr))/obsmxi(ifr),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

%%
%%  Residual plots
%%
   fig5=figure(5);hold on;box on;
   title('Data residuals real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(real(zres(:,ifr,ihyd))/obsmxr(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

   fig6=figure(6);hold on;box on;
   title('Data residuals imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(imag(zres(:,ifr,ihyd))/obsmxi(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

%     plot([1:NHYD],imag(zobs(:,ifr)),'-k');
%     plot([1:NHYD],imag(zobs(:,ifr)),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

%%
%%  Raw residual plots
%%
   fig7=figure(7);hold on;box on;
   title('Raw data residuals real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(real(zresraw(:,ifr,ihyd))/obsmxr(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

   fig8=figure(8);hold on;box on;
   title('Raw data residuals imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(imag(zresraw(:,ifr,ihyd))/obsmxi(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

%     plot([1:NHYD],imag(zobs(:,ifr)),'-k');
%     plot([1:NHYD],imag(zobs(:,ifr)),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;
%%
%%  AR plots
%%
  if(iar == 1)
   fig9=figure(9);hold on;box on;
   title('AR model series real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(real(zar(:,ifr,ihyd))/obsmxr(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

%     plot([1:NHYD],real(zobs(:,ifr)),'-k');
%     plot([1:NHYD],real(zobs(:,ifr)),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;

   fig10=figure(10);hold on;box on;
   title('AR model series imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
      NV = 400;
      vmin = -1.2;
      vmax = 1.2;
      vlim = vmin+cumsum((vmax-vmin)/NV*ones(1,NV));
      for ihyd=1:NHYD
         Nd(:,ihyd) = histc(imag(zar(:,ifr,ihyd))/obsmxi(ifr),vlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for ihyd=1:NHYD
            Nd(:,ihyd) = Nd(:,ihyd)/max(Nd(:,ihyd));
         end;
      elseif(inorm2 == 2)
         Nd = Nd/max(max(Nd));
      end;
      Nd = [zeros(length(vlim),1) Nd zeros(length(vlim),1)];
      pcolor(hydplt,vlim,Nd);shading flat;
      plot([0 49],[0 0],'-k');
      plot([0 49],[0 0],'--w');

%     plot([1:NHYD],imag(zobs(:,ifr)),'-k');
%     plot([1:NHYD],imag(zobs(:,ifr)),'--w');
      set(gca,'XLim',[0 49],'YLim',[-1.2 1.2])
      text(39,.95,[num2str(bands(ifr)) ' Hz'],'FontSize',14,'Color',[1,1,1])
      ytick = [-1.0 -0.5 0.0 0.5 1.0];
      if ((ifr == 1) | (ifr == nx+1) | (ifr == (2*nx)+1)| (ifr == (3*nx)+1))
         set(gca,'YTickLabel',ytick,'YTick',ytick);
         ylabel('Norm. pressure','FontSize',12);
      else
         set(gca,'YTickLabel',[],'YTick',ytick);
      end
      if (ifr > length(bands)-nx)
          set(gca,'XTickLabel',[10 20 30 40],'XTick',[10 20 30 40]);
          xlabel('Hydrophone No.');
      else
          set(gca,'XTickLabel',[],'XTick',[10 20 30 40]);
      end
      set(gca,'TickDir','out');
      clear Nd,vlim;

   end;
  end;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% PLOT Axx of total residuals (zres)
   %%
   fig11=figure(11);hold on;box on;
   title('AR model series imaginary part');
   hydplt = [0:NHYD+1];
   NMOD = size(zres,1);
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);
%      for imod=1:NMOD;
%        [axxr(imod,:),axxlagsr(imod,:)] = xcorr(real(zres(imod,ifr,:)),'coeff');
%        [axxi(imod,:),axxlagsi(imod,:)] = xcorr(imag(zres(imod,ifr,:)),'coeff');
%      end;
      NAX = 400;
      axmin = -.5;
      axmax = 1.2;
      axlim = axmin+cumsum((axmax-axmin)/NAX*ones(1,NAX));
      daxxlagsr = diff(axxlagsr(1,:));
%      axxlagsplt = [axxlagsr(1,1)-daxxlagsr(1), 
%                    axxlagsr(1,:), 
%                    axxlagsr(1,end)+daxxlagsr(end)]-([daxxlagsr(1),
%                    daxxlagsr(1),
%                    daxxlagsr,
%                    daxxlagsr(end)]./2.);
      axxlagsplt = axxlagsr(1,:);
      for iz=1:size(axxr,2)
        Naxxr(:,iz,ifr) = histc(axxr(:,iz,ifr),axlim);
        Naxxi(:,iz,ifr) = histc(axxi(:,iz,ifr),axlim);
        Naxxr_raw(:,iz,ifr) = histc(axxr_raw(:,iz,ifr),axlim);
        Naxxi_raw(:,iz,ifr) = histc(axxi_raw(:,iz,ifr),axlim);
      end;
      %
      % Normalize Histograms
      %
      if(inorm2 == 1)
         for iz=1:size(Naxxr,2)
            Naxxr(:,iz,ifr) = Naxxr(:,iz,ifr)/max(Naxxr(:,iz,ifr));
            Naxxi(:,iz,ifr) = Naxxi(:,iz,ifr)/max(Naxxi(:,iz,ifr));
            Naxxr_raw(:,iz,ifr) = Naxxr_raw(:,iz,ifr)/max(Naxxr_raw(:,iz,ifr));
            Naxxi_raw(:,iz,ifr) = Naxxi_raw(:,iz,ifr)/max(Naxxi_raw(:,iz,ifr));
         end;
      elseif(inorm2 == 2)
         Naxxr = Naxxr/max(max(max(Naxxr)));
         Naxxi = Naxxi/max(max(max(Naxxi)));
         Naxxr_raw = Naxxr_raw/max(max(max(Naxxr_raw)));
         Naxxi_raw = Naxxi_raw/max(max(max(Naxxi_raw)));
      end;

      imagesc(axxlagsplt',axlim,Naxxr(:,:,ifr));shading flat;

      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'-k','Linewidth',2);
      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'--w','Linewidth',2);
      set(gca,'XLim',[axxlagsr(1,1) axxlagsr(1,end)],'YLim',[axmin+.1 axmax-.1]);
      if(ifr==1 | ifr==4 | ifr==7 | ifr==10)
        ylabel('Auto corr.');
        set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
      else
        set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end;
      if(ifr>=9)
        set(gca,'XTickLabel',[-40 -20 0 20 40],'XTick',[-40 -20 0 20 40]);
        xlabel('Lag');
      else;
        set(gca,'XTickLabel',[],'XTick',[-40 -20 0 20 40]);
      end;
      set(gca,'TickDir','out');
   end;
   fig12=figure(12);hold on;box on;
   title('AR model series imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);

      imagesc(axxlagsplt',axlim,Naxxi(:,:,ifr));shading flat;

      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'-k','Linewidth',2);
      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'--w','Linewidth',2);
      set(gca,'XLim',[axxlagsr(1,1) axxlagsr(1,end)],'YLim',[axmin+.1 axmax-.1]);
      if(ifr==1 | ifr==4 | ifr==7 | ifr==10)
        ylabel('Auto corr.');
        set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
      else
        set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end;
      if(ifr>=9)
        set(gca,'XTickLabel',[-40 -20 0 20 40],'XTick',[-40 -20 0 20 40]);
        xlabel('Lag');
      else;
        set(gca,'XTickLabel',[],'XTick',[-40 -20 0 20 40]);
      end;
      set(gca,'TickDir','out');
   end;
   fig13=figure(13);hold on;box on;
   title('AR model series real part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);

      imagesc(axxlagsplt',axlim,Naxxr_raw(:,:,ifr));shading flat;

      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'-k','Linewidth',2);
      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'--w','Linewidth',2);
      set(gca,'XLim',[axxlagsr(1,1) axxlagsr(1,end)],'YLim',[axmin+.1 axmax-.1]);
      if(ifr==1 | ifr==4 | ifr==7 | ifr==10)
        ylabel('Auto corr.');
        set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
      else
        set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end;
      if(ifr>=9)
        set(gca,'XTickLabel',[-40 -20 0 20 40],'XTick',[-40 -20 0 20 40]);
        xlabel('Lag');
      else;
        set(gca,'XTickLabel',[],'XTick',[-40 -20 0 20 40]);
      end;
      set(gca,'TickDir','out');
   end;
   fig14=figure(14);hold on;box on;
   title('AR model series imaginary part');
   hydplt = [0:NHYD+1];
   for ifr=1:NBANDS;
      subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);hold on;box on;
      set(gca,'FontSize',14);

      imagesc(axxlagsplt',axlim,Naxxi_raw(:,:,ifr));shading flat;

      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'-k','Linewidth',2);
      plot([axxlagsr(1,1) axxlagsr(1,end)],[0 0],'--w','Linewidth',2);
      set(gca,'XLim',[axxlagsr(1,1) axxlagsr(1,end)],'YLim',[axmin+.1 axmax-.1]);
      if(ifr==1 | ifr==4 | ifr==7 | ifr==10)
        ylabel('Auto corr.');
        set(gca,'YTickLabel',[-0.4 0 0.4 0.8],'YTick',[-0.4 0 0.4 0.8]);
      else
        set(gca,'YTickLabel',[],'YTick',[-0.4 0 0.4 0.8]);
      end;
      if(ifr>=9)
        set(gca,'XTickLabel',[-40 -20 0 20 40],'XTick',[-40 -20 0 20 40]);
        xlabel('Lag');
      else;
        set(gca,'XTickLabel',[],'XTick',[-40 -20 0 20 40]);
      end;
      set(gca,'TickDir','out');
   end;

if(isave == 1)
%   print(fig1,'-painters','-r250',strcat(plotfile1,plotext2),'-dpng');
   print(fig1,'-painters','-r250',strcat(plotfile1,plotext3),'-depsc');
%   print(fig2,'-painters','-r250',strcat(plotfile2,plotext2),'-dpng');
   print(fig2,'-painters','-r250',strcat(plotfile2,plotext3),'-depsc');
%   print(fig3,'-painters','-r250',strcat(plotfile3,plotext2),'-dpng');
   print(fig3,'-painters','-r250',strcat(plotfile3,plotext3),'-depsc');
%   print(fig4,'-painters','-r250',strcat(plotfile4,plotext2),'-dpng');
   print(fig4,'-painters','-r250',strcat(plotfile4,plotext3),'-depsc');
%   print(fig5,'-painters','-r250',strcat(plotfile5,plotext2),'-dpng');
   print(fig5,'-painters','-r250',strcat(plotfile5,plotext3),'-depsc');
%   print(fig6,'-painters','-r250',strcat(plotfile6,plotext2),'-dpng');
   print(fig6,'-painters','-r250',strcat(plotfile6,plotext3),'-depsc');
%   print(fig7,'-painters','-r250',strcat(plotfile7,plotext2),'-dpng');
   print(fig7,'-painters','-r250',strcat(plotfile7,plotext3),'-depsc');
%   print(fig8,'-painters','-r250',strcat(plotfile8,plotext2),'-dpng');
   print(fig8,'-painters','-r250',strcat(plotfile8,plotext3),'-depsc');
  if(iar == 1)
%   print(fig9,'-painters','-r250',strcat(plotfile9,plotext2),'-dpng');
   print(fig9,'-painters','-r250',strcat(plotfile9,plotext3),'-depsc');
%   print(fig10,'-painters','-r250',strcat(plotfile10,plotext2),'-dpng');
   print(fig10,'-painters','-r250',strcat(plotfile10,plotext3),'-depsc');
  end;
%   print(fig11,'-painters','-r250',strcat(plotfile11,plotext2),'-dpng');
   print(fig11,'-painters','-r250',strcat(plotfile11,plotext3),'-depsc');
%   print(fig12,'-painters','-r250',strcat(plotfile12,plotext2),'-dpng');
   print(fig12,'-painters','-r250',strcat(plotfile12,plotext3),'-depsc');
%   print(fig13,'-painters','-r250',strcat(plotfile13,plotext2),'-dpng');
   print(fig13,'-painters','-r250',strcat(plotfile13,plotext3),'-depsc');
%   print(fig14,'-painters','-r250',strcat(plotfile14,plotext2),'-dpng');
   print(fig14,'-painters','-r250',strcat(plotfile14,plotext3),'-depsc');
end;

return;
