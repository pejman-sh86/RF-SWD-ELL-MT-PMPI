function [fig]=ffi_fix_slip_uncertainty(filename,isyn,ilin,sl_ens,xx_ln,scale,lng,sltru2);

filebase = strrep(filename,'_sample.mat','');
datfile  = strcat(filebase,'.hdf5');

NTW = 1;
NRAN = h5readatt(datfile,'/Sensitivity_kernel','N_subf_x');
NDEP = h5readatt(datfile,'/Sensitivity_kernel','N_subf_y');
NRAN = cast(NRAN,'like',1);
NDEP = cast(NDEP,'like',1);

fig=figure();
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8])
isub=1;
ir1 = 6
ir2 = 4;
for iz = 1:NDEP-2;
  for ir = ir1:NRAN-ir2;

     h(isub)=subaxis(floor(NDEP-2),floor(NRAN-(ir1+ir2-1)),isub,'Spacing',0.005,'Padding',0,...
                   'ML', 0.04,'MR',0.005,'MB',.06,'MT',.005);
    hold on;box on;
    %hist(squeeze(sl_ens(iz,ir,:)),20);
    [n1sl,slout]=hist(squeeze(sl_ens(iz,ir,:)),20);
    n1sl=n1sl/trapz(slout,n1sl);
    n1sl = [0, n1sl, 0];slout = [slout(1) slout slout(end)];
    [xx,yy]=stairs(slout,n1sl,'k');
    patch(xx,yy,[0.8,0.8,0.8]);
    stairs(slout,n1sl,'k');
    axis tight
    maxdisp(isub) = max(slout);
    mindisp(isub) = min(slout);
    ymx(isub)=max(n1sl);
    %set(gca,'XTick',[0.05:.1:.4]);
    if(iz<NDEP-2);
      set(gca,'XTickLabel',[]);
    else;
%      xlabel('Moment (Nm)');
      xlabel('Slip (m)');
    end;
    if(ir>ir1);
      set(gca,'YTickLabel',[]);
    else;
      set(gca,'YTickLabel',[]);
      ylabel('Probability');
    end;
    %axis off
    isub = isub+1;
  end;
end;
isub=1;
for iz = 1:NDEP-2;
  for ir = ir1:NRAN-ir2;
    %subplot(h(isub));
    axes(h(isub));
    set(gca,'TickDir','out','FontSize',14);
    %set(gca,'YLim',[0,max(ymx(1:end-1))]);
    set(gca,'XLim',[0,5])
    if(isyn == 1);plot([sltru2(iz,ir),sltru2(iz,ir)],[0,(ymx(isub))],'--r');end;
    if(ilin == 1);plot(xx_ln/scale(iz,ir),lng(:,iz,ir),'-b');end;
%    set(gca,'XLim',[min(min(mindisp(1:end-1,:))),max(max(maxdisp(1:end-1,:)))],'TickDir','out');
    isub = isub + 1;
  end;
end;

return;
