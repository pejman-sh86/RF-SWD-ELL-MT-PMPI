function []=auv_plot_track();
isyn = 1;
ipst = 1;
ipskip = 1;
set(0,'DefaultAxesFontSize',12);

load track.mat

if(isyn == 1);
   track_env = dlmread('track_environment.dat');
   ktru = track_env(:,1);
   htru = zeros(length(ktru),max(ktru));
   ztru = zeros(length(ktru),max(ktru));
   for i=1:length(ktru)
      idx=(ktru(i)+1)*4+1;
      sdtru(i) = track_env(i,idx);
      idx2 = [1:4:ktru(i)*4];
      htru(i,1:ktru(i)) = track_env(i,idx2+1);
      ztru(i,1:ktru(i)) = cumsum(htru(i,1:ktru(i)));
   end;
end;

NPING = size(kgl,1);
ipend = ipskip*NPING+ipst;
hmax = z(end);
hmax2 = z(end);
figw = 16;
figh = 14.;
%%
%%  PLOT WHOLE TRACK INFORMATION MAX
%%
   fig10=figure('visible','on');
   set(fig10,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
   nx = 1;
   ny = 4;
   if(isyn == 1);ny = 4;end;
   xim = 0.01;
   yim = 0.02;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   ML = .045;
   MR = .03;
   MB = .07;
   MT = .02;
   SP = .02;
   PAD = 0;
   FNT = 14;
   inter = 20;

   %%
   %% MAP ensemble velocity track profile
   %%
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,c_max);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Velocity (m/s)');
   set(gca,'ticklength',.5*get(gca,'ticklength'));


   %%
   %% MAP ensemble density track profile
   %%
   subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   hold on;box on;
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,r_max);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING]*10*.004,'TickDir','out')
   xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1.2 2.2]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Density (g/ccm)')
   set(gca,'ticklength',.5*get(gca,'ticklength'));

%%
%%  PLOT WHOLE TRACK INFORMATION MEAN
%%
   fig11=figure('visible','on');
   set(fig11,'PaperUnits','inches','PaperPosition',[0 0 figw figh])

   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,c_mean);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Velocity (m/s)');
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,r_mean);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING]*10*.004,'TickDir','out')
   xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1.2 2.2]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Density (g/ccm)')
   set(gca,'ticklength',.5*get(gca,'ticklength'));


%%
%%  PLOT WHOLE TRACK INFORMATION MEDIAN
%%
   fig12=figure('visible','on');
   set(gca,'FontSize',FNT);
   set(fig12,'PaperUnits','inches','PaperPosition',[0 0 figw figh])

   ny = 6;
   nx = 1;
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,c_mead);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Velocity (m/s)','FontSize',FNT);
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,ctr');shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1450 1750]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Velocity (m/s)','FontSize',FNT);
   set(gca,'ticklength',.5*get(gca,'ticklength'));
   
   subaxis(ny,nx,3,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,r_mead);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   %xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1.2 2.2]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Density (g/ccm)','FontSize',FNT)
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,4,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,rtr');shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   %xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[1.2 2.2]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Density (g/ccm)','FontSize',FNT)
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,5,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,a_mead);shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   %xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[0. 1.1]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Attenuation (dB/m/kHz)','FontSize',FNT)
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,6,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,atr');shading flat;
%   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING]*10*.004,'TickDir','out')
%   xlabel('Track distance (km)');
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING],'TickDir','out')
   xlabel('Track ping No.');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[0 1.1]);box on;
   colormap(jet);
   c2=colorbar;ylabel(c2,'Attenuation (dB/m/kHz)','FontSize',FNT)
   set(gca,'ticklength',.5*get(gca,'ticklength'));

%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  PLOT WHOLE TRACK INFORMATION k logL sigma
%%
   fig13=figure('visible','on');
   set(fig13,'PaperUnits','inches','PaperPosition',[0 0 figw figh])

   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],limgl,kgl');
   if(isyn == 1);
      stairs([1:NPING]',ktru(1:NPING),'k','LineWidth',2);
      stairs([1:NPING]',ktru(1:NPING),'--w','LineWidth',2);
   end;
   set(gca,'XTickLabel',[])
   ylabel('No. interfaces');
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[],'TickDir','out')
   xlabel('Track distance (km)');
   set(gca,'Layer','top');box on;
   set(gca,'YLim',[1 14],'CLim',[0 max(max(kgl))],'TickDir','out');box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'k probability')
   set(gca,'ticklength',.5*get(gca,'ticklength'));

   subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],logLlimgl,logLgl');shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING]*10*.004,'TickDir','out')
   xlabel('Track distance (km)');
   ylabel('log(L)');
   set(gca,'Layer','top');box on;
   set(gca,'YLim',[logLlimgl(1) logLlimgl(end)],'CLim',[0 max(max(logLgl))/4],'TickDir','out');box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'log(L) probability')
   set(gca,'ticklength',.5*get(gca,'ticklength'));
   set(gca,'Layer','top','YDir','reverse');box on;

%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  PLOT WHOLE TRACK INFORMATION INTERFACE PROBABILITY
%%
   figw = 14;
   figh = 6.6667;
   fig15=figure('visible','on');
   set(fig15,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   %%
   %% Interface probability
   %%
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   hold on;box on;
   imagesc([1:NPING],zi,hgl');shading flat;
   set(gca,'XTick',[1:inter:NPING],'XTickLabel',[0:inter:NPING]*10*.004,'TickDir','out')
   xlabel('Track distance (km)');
   set(gca,'YTick',[0:2:20])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[1 NPING],'YLim',[0 hmax2],'CLim',[0 0.1]);box on;
   set(gca,'ticklength',.5*get(gca,'ticklength'));
   colormap(jet);
   c1=colorbar;ylabel(c1,'Interface probability')

   if(isyn == 1);
      for i=1:length(ktru)
         for j=1:ktru(i)
            plot([i i+1],[ztru(i,j),ztru(i,j)]-.05,'-w','Linewidth',2);
%            plot([i i+1],[ztru(i,j),ztru(i,j)]-.05,':k','Linewidth',2);
         end;
      end;
   end;

%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  PLOT WHOLE TRACK INFORMATION UNCERTAINTY 95% HPDs
%%
if(1==2);
   nx = 1;
   ny = 3;
   xim = 0.01;
   yim = 0.02;
   xymarg = [0.07 0.04 0.04 0.14];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

   fig14=figure('visible','on');
   set(fig14,'PaperUnits','inches','PaperPosition',[0 0 figw figh])
   %%
   %% Velocity
   %%
   %subplot('Position',[loc(1,1) loc(2,1) spw sph]);hold on;box on;
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,nfc(:,:,2)-nfc(:,:,1));shading flat;
   set(gca,'XTickLabel',[])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[ipst ipend],'YLim',[0 hmax2],'CLim',[0 100]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Velocity uncertainty (m/s)')

   %%
   %% Density 
   %%
   %subplot('Position',[loc(1,2) loc(2,2) spw sph]);hold on;box on;
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,nfr(:,:,2)-nfr(:,:,1));shading flat;
   set(gca,'XTickLabel',[])
   ylabel('Depth (m)');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[ipst ipend],'YLim',[0 hmax2],'CLim',[0 0.5]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Density uncertainty (g/ccm)')


   %%
   %% Attenuation 
   %%
   %subplot('Position',[loc(1,3) loc(2,3) spw sph]);hold on;box on;
   subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
   set(gca,'FontSize',FNT);
   imagesc([1:NPING],z,nfa(:,:,2)-nfa(:,:,1));shading flat;
   ylabel('Depth (m)');
   xlabel('Ping no.');
   set(gca,'Layer','top','YDir','reverse');box on;
   set(gca,'XLim',[ipst ipend],'YLim',[0 hmax2],'CLim',[0 1]);box on;
   colormap(jet);
   c1=colorbar;ylabel(c1,'Attenuation uncertainty (dB/m/kHz)')
end;

%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  SAVE PLOTS 
%%

%   saveas(fig10,'track_results_v_r_max.fig','fig');
%   saveas(fig11,'track_results_v_r_mean.fig','fig');
%   saveas(fig12,'track_results_v_r_mead.fig','fig');
%   saveas(fig13,'track_results_k_logL.fig','fig');

   print(fig10,'-painters','-r300','track_results_v_r_max.png','-dpng');
   print(fig11,'-painters','-r300','track_results_v_r_mean.png','-dpng');
   print(fig12,'-painters','-r300','track_results_v_r_mead.png','-dpng');
   print(fig13,'-painters','-r300','track_results_k_logL.png','-dpng');
%   print(fig14,'-painters','-r300','track_results_uncertainty.png','-dpng');
   print(fig15,'-painters','-r300','track_results_int_prob.png','-dpng');

%   print(fig10,'-painters','-r300','track_results_v_r_max.eps','-depsc');
%   print(fig11,'-painters','-r300','track_results_v_r_mean.eps','-depsc');
   print(fig12,'-painters','-r300','track_results_v_r_mead.eps','-depsc');
%   print(fig13,'-painters','-r300','track_results_k_logL.eps','-depsc');
%   print(fig14,'-painters','-r300','track_results_uncertainty.eps','-depsc');
%   print(fig15,'-painters','-r300','track_results_int_prob.eps','-depsc');

return;
