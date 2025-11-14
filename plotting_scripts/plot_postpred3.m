
filebase = 'HON';
NDATASET = 3;
iadd_mappred = 1;
posterior_color = [.7 .7 .7];
posterior_color2 = [0. 1. 0. .004];
preds = importdata(strcat(filebase,'_postpred.dat'));
parfile = strcat(filebase,'_parameter.dat');
% predsRF_SM = importdata(strcat(filebase,'_SMpostpred_RF.dat'));

iwresRF = 0; %weighted residual
imake_CdRF = 0;

iwresSWD = 0;
imake_CdSWD = 0;

iwresELL = 0;
imake_CdELL = 0;

iwresMT = 0;
imake_CdMT = 0;

%load SWDerror_stds.mat   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
[IMAP,ICOV,iar,i_varpar,irv,itr,iswd, iell, imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL, NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVs,dVpVs, sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile(parfile);
NRF1 = ntr;
dt = sampling_dt;
%%

if irv == -1 
    NRF2 = NRF1; 
else
    NRF2 = 0;
end

if iswd == 1 
    NMODE2 = NMODE; 
else
    NMODE2 = 0;
end

if iell == 1 
    NMODE_ELL2 = NMODE_ELL; 
else
    NMODE_ELL2 = 0;
end

%%
istart = 1;
iend = NRF2*NTIME;
if(irv==-1); predrf = preds(:, istart:iend); end
% if(irv==-1); predrf_SM = predsRF_SM(:, istart:iend); end   %%%%%%%%%%%%%%%%%%%%%%
istart = iend + 1;
iend = iend + NMODE2*NDAT_SWD;
if(iswd==1); predswd = preds(:, istart:iend); end
istart = iend + 1;
iend = iend + NMODE_ELL2*NDAT_ELL;
if(iell==1); predell = preds(:, istart:iend); end
istart = iend + 1;
iend = iend + 2*NDAT_MT;
if(imt==1); predmt = preds(:, istart:iend); end

%%
%imt = 0; iswd=1; irv=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if irv == -1
    datarf = importdata(strcat(filebase, '_RF.txt')); % 2xNTIME array
%     datarf = importdata(strcat(filebase, '_RF_SMbest.txt')); % 2xNTIME array
%     datarf_SMbest = importdata(strcat(filebase, '_RF_SMbest.txt')); % 2xNTIME array
%     obsrf_SMbest = datarf_SMbest(1,:); % 1xNTIME array
%     datarf_mildHVZ_3_9 = importdata(strcat(filebase, '_RF_mildHVZ_3.9.txt')); % 2xNTIME array
%     obsrf_mildHVZ_3_9 = datarf_mildHVZ_3_9(1,:); % 1xNTIME array
%     datarf_sharpHVZ = importdata(strcat(filebase, '_RF_sharpHVZ.txt')); % 2xNTIME array
%     obsrf_sharpHVZ = datarf_sharpHVZ(1,:); % 1xNTIME array
%     datarf_LVZ = importdata(strcat(filebase, '_RF_LVZ.txt')); % 2xNTIME array
%     obsrf_LVZ = datarf_LVZ(1,:); % 1xNTIME array
    obsrf = datarf(1,:); % 1xNTIME array
    %resrf = obsrf - predrf;
    if iwresRF == 1
        if imake_CdRF == 1
            CdRF = importdata(strcat(filebase, '_Cd.dat'));
        else
            load covmatRF.mat;
            CdRF = F(1).Cd;
        end
        Lrf = chol(CdRF,'lower');
        resrf = transpose (inv(Lrf) * resrf' );
    end
end

if iswd == 1
    dataswd = importdata(strcat(filebase, '_SWD.dat')); % NDAT_SWDx2 array
%     dataswd = importdata(strcat(filebase, '_SWD_mildHVZ_3.9.dat')); % NDAT_SWDx2 array
%     dataswd_mildHVZ_3_9 = importdata(strcat(filebase, '_SWD_mildHVZ_3.9.dat')); % NDAT_SWDx2 array
%     obsswd_mildHVZ_3_9 = dataswd_mildHVZ_3_9(:,2); % NDAT_SWDx1 array
%     dataswd_sharpHVZ = importdata(strcat(filebase, '_SWD_sharpHVZ.dat')); % NDAT_SWDx2 array
%     obsswd_sharpHVZ = dataswd_sharpHVZ(:,2); % NDAT_SWDx1 array
%     dataswd_LVZ = importdata(strcat(filebase, '_SWD_LVZ.dat')); % NDAT_SWDx2 array
%     obsswd_LVZ = dataswd_LVZ(:,2); % NDAT_SWDx1 array
    obsswd = dataswd(:,2); % NDAT_SWDx1 array
    Tswd = dataswd(:,1); % NDAT_SWDx1 array
    %resswd = obsswd' - predswd;
     if iwresSWD == 1
         if imake_CdSWD == 1
             CdSWD = importdata(strcat(filebase, '_CdSWD.dat'));
         else
             load covmatSWD.mat;
             CdSWD = F(1).Cd;
         end
        Lswd = chol(CdSWD,'lower');
        resswd = transpose( inv(Lswd) * resswd' );
    end
end

if iell == 1
    dataell = importdata(strcat(filebase, '_ELL.dat')); % NDAT_ELLx2 array
    obsell = dataell(:,2); % NDAT_ELLx1 array
    Tell = dataell(:,1); % NDAT_ELLx1 array
    resell = obsell' - predell;
    if iwresELL == 1
        if imake_CdELL == 1
            CdELL = importdata(strcat(filebase, '_CdELL.dat'));
        else
            load covmatELL.mat;
            CdELL = F(1).Cd;
        end
        Lell = chol(CdELL,'lower');
        resell = transpose( inv(Lell) * resell' );
    end
end

if imt == 1
    datamt = importdata(strcat(filebase, '_MT.dat')); % NDAT_MTx3 array
    obsmt = reshape(datamt(:,2:3),[1,2*NDAT_MT]); % 1x2NDAT_MT array
    %obsmt = -obsmt;
    fmt = datamt(:,1); % NDAT_MTx1 array
    resmt = obsmt - predmt;
    complex_resmt = complex(resmt(:,1:NDAT_MT),resmt(:,NDAT_MT+1:2*NDAT_MT));
    if iwresMT == 1
        if imake_CdMT == 1
            AAA = importdata(strcat(filebase, '_CdMT.dat'));
            CdMT = complex( AAA(:,1:NDAT_MT), AAA(:,NDAT_MT+1:2*NDAT_MT) );
            CdMT(1:NDAT_MT, 1:NDAT_MT)=real(CdMT(1:NDAT_MT, 1:NDAT_MT));
        else
        load covmatMT.mat;
        CdMT = F(1).Cd;
        end
        Lmt = chol(CdMT,'lower');
        complex_resmt = transpose( inv(Lmt) * (complex_resmt.') );
    end
end

%%

fig100=figure;
   nx = NDATASET;
   ny = 2;
   xim = 0.05;
   yim = 0.08;
   xymarg = [0.1 0.04 0.04 0.1];
   [loct,spwt,spht] = get_loc(nx,ny,xim,yim,xymarg);
   loct(2,:) = loct(2,:) - 0.05; %real
   spwt2 = spwt;
   spwt = spwt*9./10.;
   nx = 1;
   [locb,spwb,sphb] = get_loc(nx,ny,xim,yim,xymarg);
   locb(2,:) = locb(2,:) - 0.05;  %real
   spwb = spwb - spwt2*1./10.;

   p21 = subplot('Position',[loct(1,1) loct(2,1) spwt spht]);
   hold on; box off;
   p22 = subplot('Position',[loct(1,2) loct(2,2) spwt spht]);
   hold on; box off;
   p23 = subplot('Position',[loct(1,3) loct(2,3) spwt spht]);
   hold on; box off;
%    p11 = subplot('Position',[locb(1,1+1) locb(2,1+1)+sphb/2-.04 spwb sphb/2]);
   p11 = subplot('Position',[locb(1,1+1) locb(2,1+1)+sphb/2-.06 spwb sphb/2]);
   hold on; box off;

   %% Plot posterior predictions

%    if irv == -1
%        for is = 1:size(predrf_SM,1)
%             subplot(p11)
%             hold on; box on
%             plot((0:NTIME-1)*dt-shift, predrf_SM(is,:),'Color', posterior_color2)
% %             scatter((0:NTIME-1)*dt-shift, predrf_SM(is,:),.1,...
% %                    'MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%        end
%    end

   for is = 1:size(preds,1)

       if imt == 1
        subplot(p21)
        hold on; box on
        plot(log10(1./fmt), predmt(is,1:NDAT_MT), 'Color', posterior_color)
    
        subplot(p22)
        hold on; box on
        plot(log10(1./fmt), predmt(is, NDAT_MT+1:2*NDAT_MT), 'Color', posterior_color)
        end

        if iswd == 1
            subplot(p23)
            hold on; box on
            p3=plot(Tswd, predswd(is,:), 'Color', posterior_color);
        end
%             if iswd == 1
%                         subplot(p21)
%                         hold on; box on
%                         p3=plot(Tswd, predswd(is,:), 'Color', posterior_color);
%             end
%             if iswd == 1
%                         subplot(p22)
%                         hold on; box on
%                         p3=plot(Tswd, predswd(is,:), 'Color', posterior_color);
%             end

        if irv == -1
            subplot(p11)
            hold on; box on
%             if(is<=size(predrf_SM,1)); plot((0:NTIME-1)*dt-shift, predrf_SM(is,:),'Color', posterior_color2); end
            p3=plot((0:NTIME-1)*dt-shift, predrf(is,:),'Color', posterior_color);
%             scatter((0:NTIME-1)*dt-shift, predrf(is,:),.1,...
%                    'MarkerFaceColor','r','MarkerEdgeColor','r',...
%                    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
        end

   end

%    if irv == -1       %%it was used for plot
%        for is = 1:size(predrf_SM,1)
%             subplot(p11)
%             hold on; box on
%             p4=plot((0:NTIME-1)*dt-shift, predrf_SM(is,:),'Color', posterior_color2);
% %                    plot((0:NTIME-1)*dt-shift, predrf_SM(is,:),.1,...
% %                    'MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
%        end
%    end

   %% Add observed values
    font_size = 30;
    font_size2 = 25;
    marker_size = 8;
    marker_size2 = 4;
    subplot(p21)
    hold on; box on
    if imt == 1
%         subplot(p21)
%         hold on; box on
        plot(log10(1./fmt), obsmt(1:NDAT_MT), 'ok', 'MarkerSize', marker_size)
%         semilogx(1./fmt, obsmt(1:NDAT_MT), 'ok', 'MarkerSize', marker_size)
%         plot(log10(1./fmt), obsmt(1:NDAT_MT), 'ow', 'MarkerSize', marker_size2)

        set(gca,'Fontsize',font_size);
        %set(gca, 'Xlim', [-0.2 3.5])
        %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
        %set(gca, 'XtickLabel', [])
%         set(gca, 'Ylim', [-.001 .03])   %sim7
%         set(gca, 'Ylim', [0. .009])   %sim8
%         set(gca, 'Ylim', [-.0015 .04])   %sim9
        set(gca, 'Ylim', [0. .009])   %sim10
    %     set(gca, 'YTick', [0.: 0.04: .12])   %sim
        %ylabel('Apparent resistivity (\Omega-m)')
        xlabel('log_{10} Period(s)')
        ylabel('REAL(Z)')
        %legend('fitted values', 'obsereved data')
    end

    subplot(p22)
    hold on; box on
    if imt == 1
%         subplot(p22)
%         hold on; box on
        plot(log10(1./fmt), obsmt(NDAT_MT+1:2*NDAT_MT), 'ok', 'MarkerSize', marker_size)
%         semilogx(1./fmt, obsmt(NDAT_MT+1:2*NDAT_MT), 'ok', 'MarkerSize', marker_size)
%         plot(log10(1./fmt), obsmt(NDAT_MT+1:2*NDAT_MT), 'ow', 'MarkerSize', marker_size2)

        set(gca,'Fontsize',font_size);
        %set(gca, 'Xlim', [-0.2 3.5])
        %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
        %set(gca, 'XtickLabel', [])
%         set(gca, 'Ylim', [-.001 .03])   %sim7
%         set(gca, 'Ylim', [0. 2.8e-3])   %sim8
%         set(gca, 'Ylim', [-.0015 .08])   %sim9
        set(gca, 'Ylim', [0. 2.8e-3])   %sim10
    %     set(gca, 'YTick', [0.: 0.04: .12])   %sim
        %ylabel('Phase(degree)')
        xlabel('log_{10} Period(s)')
        ylabel('IMAG(Z)')
        %legend('fitted values', 'obsereved data')
    end

    subplot(p23)
    hold on; box on
    if iswd == 1
%         subplot(p23)
%         hold on; box on
        p1=plot(Tswd, obsswd,'ok', 'MarkerSize', marker_size);
%         p2=plot(Tswd, obsswd_mildHVZ_3_9,'x-r', 'MarkerSize', marker_size2);
%         p2=plot(Tswd, obsswd_sharpHVZ,'x-r', 'MarkerSize', marker_size2);
%         p2=plot(Tswd, obsswd_LVZ,'x-r', 'MarkerSize', marker_size2);

        set(gca,'Fontsize',font_size);
        set(gca, 'Xlim', [0.0 260.])   %sim7,8 and real
        %set(gca, 'XtickLabel', [])
%         set(gca, 'Ylim', [2.8 5.])    %sim7
%         set(gca, 'Ylim', [2.68 5.25])    %sim8
%           set(gca, 'Ylim', [2.8 5.])    %sim9
        set(gca, 'Ylim', [2.68 5.5])    %sim10
%         set(gca, 'Ylim', [2.5 5.])    %sim7
        xlabel('Period(s)')
        ylabel('Phase Velocity(km/s)','FontSize',font_size2)
%         legend('fitted values', 'obsereved data')
%         legend([p1 p2 p3], {'observed data','SWD data by HVZ and LVZ','MT-SWD-RF posterior predictions'},'FontSize',18,'Location','best')
    end

    subplot(p11)
    hold on; box on
    if irv == -1
%          subplot(p11)
%        hold on; box on
%         shift = 2;
        p1= plot((0:NTIME-1)*dt-shift, obsrf, '*k');%, 'MarkerSize', marker_size)
%        p2=plot((0:NTIME-1)*dt-shift, obsrf_SMbest, 'g');%, 'MarkerSize', marker_size)
%         p2=plot((0:NTIME-1)*dt-shift, obsrf_mildHVZ_3_9, 'g');%, 'MarkerSize', marker_size)
%         p2=plot((0:NTIME-1)*dt-shift, obsrf_sharpHVZ, 'g');%, 'MarkerSize', marker_size)
%         p2=plot((0:NTIME-1)*dt-shift, obsrf_LVZ, 'g');%, 'MarkerSize', marker_size)
%         set(gca, 'Ylim', [-0.2 .58]) %sim8
%         set(gca, 'Ylim', [-0.23 .52])  %sim7
%         set(gca, 'Ylim', [-0.2 .48]) %sim9
        set(gca, 'Ylim', [-0.25 .58]) %sim10
%         set(gca, 'Ylim', [-0.15 .54])  %real
        ylim = get(gca,'Ylim');
        hold on
        plot([0. 0.],[ylim(1) ylim(2)],'--k')

        set(gca,'Fontsize',font_size);
%         set(gca, 'Xlim', [-2. 28.])  %sim7, sim9
        set(gca, 'Xlim', [-2. 28.2])  %sim8, sim10
%         set(gca, 'Xlim', [-1.5 20])   %real
        %set(gca, 'XtickLabel', [])
        xlabel('time(s)')
        ylabel('RF')
    %     legend([p1 p2 p3], {'observed RF','RF by MT-SWD Best fitting model','MT-SWD-RF posterior predictive'},'FontSize',17)
%         legend([p1 p4 p3], {'observed RF','RFs by MT-SWD posterior samples','RFs by RF-SWD-MT posterior samples'},'FontSize',18)
%         legend([p1 p2], {'observed RF','RF by a model with crustal HVZ and LVZ'},'FontSize',18)
% legend([p1 p2], {'observed RF','RF by a model with mild upper crustal HVZ'},'FontSize',18)
%           legend([p1 p2], {'observed RF','RF by a model with a sharp upper crustal HVZ'},'FontSize',18)
        %legend('fitted values', 'obsereved data')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iadd_mappred

    if imt == 1
        mappred_mt = importdata(strcat(filebase, '_mappredMT.dat'));
        subplot(p21)
        hold on; box on; 
        p3r = plot(log10(1./fmt), mappred_mt(1:NDAT_MT), '--r');
        subplot(p22)
        hold on; box on;
        p3i = plot(log10(1./fmt), mappred_mt(NDAT_MT+1:2*NDAT_MT), '--r');
    end   

    if iswd == 1
        mappred_swd = importdata(strcat(filebase, '_mappredSWD.dat'));
        subplot(p23)
        hold on; box on;
        p3 = plot(Tswd, mappred_swd,'--r');
    end    

    if irv == -1
        mappred_rf = importdata(strcat(filebase, '_mappred.dat'));
        subplot(p11)
        hold on; box on;
        p3 = plot((0:NTIME-1)*dt-shift, mappred_rf, '--r');
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
if iswd==1
    subplot(p22)
        hold on; box on
        p1=plot(Tswd, obsswd,'ok', 'MarkerSize', marker_size);
%         p2=plot(Tswd, obsswd_mildHVZ_3_9,'x-r', 'MarkerSize', marker_size2);
%         p2=plot(Tswd, obsswd_sharpHVZ,'x-r', 'MarkerSize', marker_size2);
        p2=plot(Tswd, obsswd_LVZ,'x-r', 'MarkerSize', marker_size2);

        set(gca,'Fontsize',font_size);
        set(gca, 'Xlim', [0.0 260.])   %sim and real
        %set(gca, 'XtickLabel', [])
        %set(gca, 'Ylim', [2.8 5.])    %sim
        xlabel('Period(s)')
        ylabel('Phase Velocity(km/s)','FontSize',font_size2)

subplot(p23)
        hold on; box on
        p1=plot(Tswd, SWDerror_stds,'ok', 'MarkerSize', marker_size);
%         p2=plot(Tswd, abs(obsswd_mildHVZ_3_9-obsswd),'x-r', 'MarkerSize', marker_size2);
%         p2=plot(Tswd, abs(obsswd_sharpHVZ-obsswd),'x-r', 'MarkerSize', marker_size2);
        p2=plot(Tswd, abs(obsswd_LVZ-obsswd),'x-r', 'MarkerSize', marker_size2);

        set(gca,'Fontsize',font_size);
        set(gca, 'Xlim', [0.0 260.])   %sim and real
        %set(gca, 'XtickLabel', [])
        %set(gca, 'Ylim', [2.8 5.])    %sim
        xlabel('Period(s)')
%         ylabel('Phase Velocity(km/s)')
        ylabel('Absolute residual(km/s)','FontSize',font_size2)
        legend([p1 p2], {'SWD estimated noise magnitude','model residual'},'FontSize',16)
end
end

