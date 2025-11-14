function []=ffi_tsunami_slip_uncertainty(filename,sl_ens,sltru2,NX1,NX2,NY1,NY2,ieps,isyn);

filebase = strrep(filename,'_sample.mat','');
datfile  = strcat(filebase,'.hdf5');

%load ensemble_sl.mat;
NTW  = h5readatt(datfile,'/Sensitivity_kernel','num_tw');
NRAN = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NTW  = cast(NTW,'like',1);
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);

plot_tw1 = strcat(filebase,'_sl1_marg_',num2str(NX2),'_',num2str(NY2),'.');
if(NTW>1);
  plot_tw2 = strcat(filebase,'_sl2_marg.');
end;
if(NTW>2);
  plot_tw3 = strcat(filebase,'_sl3_marg.');
end;
disp([NX1,NX2,NY1,NY2]);
size(sl_ens)
NX = (NX2-NX1)+1;
NY = (NY2-NY1)+1;
for itw=1:NTW;
    fig(itw)=figure();
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 24 14])
    isub=1;
    for iz = NY1:NY2;
      for ir = NX1:NX2;
        h(isub)=subaxis(NY,NX,isub,'Spacing',0.005,'Padding',0,...
                       'ML', 0.04,'MR',0.005,'MB',.06,'MT',.005);
        hold on;box on;
        %hist(squeeze(sl_ens(iz,ir,:,itw)),20);
        [n1sl,slout]=hist(squeeze(sl_ens(iz,ir,:,itw)),20);
        n1sl=n1sl/trapz(slout,n1sl);
        n1sl = [0, n1sl, 0];slout = [slout(1) slout slout(end)];
        [xx,yy]=stairs(slout,n1sl,'k');
        patch(xx,yy,[0.8,0.8,0.8]);
        stairs(slout,n1sl,'k');
        axis tight
        maxdisp(isub,itw) = max(slout);
        mindisp(isub,itw) = min(slout);
        ymx(isub,itw)=max(n1sl);
        set(gca,'XTick',[-10:2:10]);
        if(iz<NY2);
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
%end;
%for itw=1:NTW;
%    set(0, 'currentfigure', fig(itw));
    isub=1;
    for iz = NY1:NY2;
      for ir = NX1:NX2;
        %subplot(h(isub));
        axes(h(isub));
        set(gca,'YLim',[0,max(ymx(1:end-1,itw))],'TickDir','out');
        set(gca,'XLim',[-4,8.5])
%        set(gca,'XLim',[-1,2.5])
        plot([0,0],[0,max(ymx(1:end-1,itw))],'--k')
        if(isyn == 1);
          plot([sltru2(iz,ir,itw),sltru2(iz,ir,itw)],[0,max(ymx(1:end-1,itw))],'--r')
        end;
%        set(gca,'XLim',[min(min(mindisp(1:end-1,:))),max(max(maxdisp(1:end-1,:)))],'TickDir','out');
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(1)+(xlim(2)-xlim(1))/12,ylim(2)-(ylim(2)-ylim(1))/12,strcat('ix=',num2str(ir),', iy=',num2str(iz)));
        isub = isub + 1;
      end;
    end;
end;

plotext1    = 'fig';
plotext2    = 'png';
plotext3    = 'eps';

strcat(plot_tw1,plotext2)
print(fig(1),'-painters','-r250',strcat(plot_tw1,plotext2),'-dpng');
saveas(fig(1),strcat(plot_tw1,plotext1),'fig');
if(NTW > 1);
  strcat(plot_tw2,plotext2)
  print(fig(2),'-painters','-r250',strcat(plot_tw2,plotext2),'-dpng');
  saveas(fig(2),strcat(plot_tw2,plotext1),'fig');
end;
if(NTW > 2);
  print(fig(3),'-painters','-r250',strcat(plot_tw3,plotext2),'-dpng');
  saveas(fig(3),strcat(plot_tw3,plotext1),'fig');
end;
if(ieps == 1);
  print(fig(1),'-painters','-r250',strcat(plot_tw1,plotext3),'-depsc');
  if(NTW > 1);
    print(fig(2),'-painters','-r250',strcat(plot_tw2,plotext3),'-depsc');
  end;
  if(NTW > 2);
    print(fig(3),'-painters','-r250',strcat(plot_tw3,plotext3),'-depsc');
  end;
end;

return;
