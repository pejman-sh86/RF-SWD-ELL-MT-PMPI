function []=ffi_slab_mtw_slip_uncertainty(filename,sltot_hst2,bins,NR1,NR2,NZ1,NZ2,ieps,isyn);

filebase = strrep(filename,'_sample.mat','');
datfile  = strcat(filebase,'.hdf5');

%load ensemble_sl.mat;
NRAN   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP   = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);

plot_tw1 = strcat(filebase,'_sltot_marg_',num2str(NR2),'_',num2str(NZ2),'.');
disp([NR1,NR2,NZ1,NZ2]);
NR = (NR2-NR1)+1;
NZ = (NZ2-NZ1)+1;
fig=figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 24 14])
isub=1;
for iz = NZ1:NZ2;
  for ir = NR1:NR2;
    h(isub)=subaxis(NZ,NR,isub,'Spacing',0.005,'Padding',0,...
                   'ML', 0.04,'MR',0.005,'MB',.06,'MT',.005);
    hold on;box on;
    [xx,yy]=stairs(bins(4,3:end),squeeze(sltot_hst2(iz,ir,3:end)),'k');
    %patch(xx,yy,[0.8,0.8,0.8]);
    stairs(bins(4,:),squeeze(sltot_hst2(iz,ir,:)),'k');
    axis tight
    maxdisp(isub) = max(bins(4,:));
    mindisp(isub) = min(bins(4,:));
    ymx(isub)=max(sltot_hst2(iz,ir,:));
    set(gca,'XTick',[0:10:110]);
    if(iz<NZ2);
      set(gca,'XTickLabel',[]);
    else;
      xlabel('Displacement (m)');
    end;
    %if(ir>1);
      set(gca,'YTickLabel',[]);
    %end;
    %axis off
    isub = isub+1;
  end;
end;
isub=1;
for iz = NZ1:NZ2;
  for ir = NR1:NR2;
    %subplot(h(isub));
    axes(h(isub));
%    set(gca,'YLim',[0,max(ymx(1:end-1))],'TickDir','out');
    set(gca,'XLim',[0,90])
    %plot([0,0],[0,max(ymx(1:end-1))],'--k')
%    if(isyn == 1);
%      plot([sltru2(iz,ir,itw),sltru2(iz,ir,itw)],[0,max(ymx(1:end-1,itw))],'--r')
%    end;
    ylim=get(gca,'ylim');
    xlim=get(gca,'xlim');
    text(xlim(1)+(xlim(2)-xlim(1))/12,ylim(2)-(ylim(2)-ylim(1))/12,strcat('ix=',num2str(ir),', iy=',num2str(iz)));
    isub = isub + 1;
  end;
end;

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

%strcat(plot_tw1,plotext2)
%print(fig(1),'-painters','-r250',strcat(plot_tw1,plotext2),'-dpng');
%saveas(fig(1),strcat(plot_tw1,plotext1),'fig');
%if(ieps == 1);
%  print(fig(1),'-painters','-r250',strcat(plot_tw1,plotext3),'-depsc');
%end;

return;
