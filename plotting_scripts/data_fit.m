irv = -1;
dt = 0.1;
iwresRF = 0; %weighted residual
imake_CdRF = 1;

iswd = 1;
iwresSWD = 0;
imake_CdSWD = 1;

iell = 0;
iwresELL = 1;
imake_CdELL = 1;

imt = 1;
iwresMT = 0;
imake_CdMT = 1;

filebase = 'HON';

if irv == -1
%     NTIME = 349; 
    predrf = dlmread(strcat(filebase, '_mappred.dat')); % 1xNTIME array
    NTIME = length(predrf);
    datarf = dlmread(strcat(filebase, '_RF.txt')); % 2xNTIME array
    obsrf = datarf(1,:); % 1xNTIME array
    resrf = obsrf - predrf;
    if iwresRF == 1
        if imake_CdRF == 1
            CdRF = importdata(strcat(filebase, '_Cd.dat'));
        else
            load covmatRF.mat;
            CdRF = F(1).Cd;
        end
        Lrf = chol(CdRF,'lower');
        resrf = inv(Lrf) * resrf';
    end
end

if iswd == 1
    predswd = dlmread(strcat(filebase, '_mappredSWD.dat')); % 1xNDAT_SWD array
    dataswd = dlmread(strcat(filebase, '_SWD.dat')); % NDAT_SWDx2 array
    obsswd = dataswd(:,2); % NDAT_SWDx1 array
    Tswd = dataswd(:,1); % NDAT_SWDx1 array
    resswd = obsswd - predswd';
    NDAT_SWD = length(Tswd);
     if iwresSWD == 1
         if imake_CdSWD == 1
             CdSWD = importdata(strcat(filebase, '_CdSWD.dat'));
         else
             load covmatSWD.mat;
             CdSWD = F(1).Cd;
         end
        Lswd = chol(CdSWD,'lower');
        resswd = inv(Lswd) * resswd;
    end
end

if iell == 1
    predell = dlmread(strcat(filebase, '_mappredELL.dat')); % 1xNDAT_ELL array
    dataell = dlmread(strcat(filebase, '_ELL.dat')); % NDAT_ELLx2 array
    obsell = dataell(:,2); % NDAT_ELLx1 array
    Tell = dataell(:,1); % NDAT_ELLx1 array
    resell = obsell - predell';
    NDAT_ELL = length(Tell);
    if iwresELL == 1
        if imake_CdELL == 1
            CdELL = importdata(strcat(filebase, '_CdELL.dat'));
        else
            load covmatELL.mat;
            CdELL = F(1).Cd;
        end
        Lell = chol(CdELL,'lower');
        resell = inv(Lell) * resell;
    end
end

if imt == 1
    predmt = dlmread(strcat(filebase, '_mappredMT.dat')); % 1x2NDAT_MT array
    datamt = dlmread(strcat(filebase, '_MT.dat')); % NDAT_MTx3 array
    obsmt = reshape(datamt(:,2:3),size(predmt)); % 1x2NDAT_MT array
    %obsmt = -obsmt;
    fmt = datamt(:,1); % NDAT_MTx1 array
    resmt = obsmt - predmt;
    NDAT_MT = length(fmt);
    complex_resmt = complex(resmt(1:NDAT_MT),resmt(NDAT_MT+1:2*NDAT_MT));
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
        complex_resmt = inv(Lmt) * (complex_resmt.');
    end
end
%% 
    
left = 0.05; bottom = 0.10;
w1 = 0.275; h1 = 0.18; spw1 = 0.03; sph1 = 0.03;
w2 = 0.275; h2 = 0.3; spw2 = 0.03; sph2 = 0.06;
w3 = 0.91; h3 = 0.11; spw3 = 0.01; sph3 = 0.03;
w4 = 0.91; h4 = 0.17;

loc1 = [left,    left+w1+spw1, left+2*w1+2*spw1; ...
        bottom,      bottom,       bottom];

loc2 = [left,           left+w2+spw2,       left+2*w2+2*spw2; ...
        bottom+h1+sph1,    bottom+h1+sph1,   bottom+h1+sph1   ];
loc2(2,:) = loc2(2,:) - 0.018;

loc3 = [left; ...
          bottom+h1+sph1+h2+sph2];

loc4 = [left; ...
         bottom+h1+sph1+h2+sph2+h3+sph3];
loc4(2,:) = loc4(2,:) - 0.018;

fig100 = figure;hold on; box on;
if imt == 1
    p11 = subplot('Position',[loc1(1,1) loc1(2,1) w1 h1]);
    hold on; box off;
    p12 = subplot('Position',[loc1(1,2)+0.01 loc1(2,2) w1 h1]);
    hold on; box off;
end
if iswd ==1
    p13 = subplot('Position',[loc1(1,3)+0.02 loc1(2,3) w1 h1]);
    hold on; box off;
end
if iell ==1
    p13 = subplot('Position',[loc1(1,3)+0.02 loc1(2,3) w1 h1]);
    hold on; box off;
end
if imt == 1
    p21 = subplot('Position',[loc2(1,1) loc2(2,1)+0.06 w2 h2]);
    hold on; box off;
    p22 = subplot('Position',[loc2(1,2)+0.01 loc2(2,2)+0.06 w2 h2]);
    hold on; box off;
end
if iswd == 1
    p23 = subplot('Position',[loc2(1,3)+0.02 loc2(2,3)+0.06 w2 h2]);
    hold on; box off;
end
if iell == 1
    p23 = subplot('Position',[loc2(1,3)+0.02 loc2(2,3)+0.06 w2 h2]);
    hold on; box off;
end
if irv == -1
    p31 = subplot('Position',[loc3(1,1) loc3(2,1) w3 h3]);
    hold on; box off;

    p41 = subplot('Position',[loc4(1,1) loc4(2,1) w4 h4]);
    hold on; box off;
end
    %% 
if imt == 1
    subplot(p11)
    hold on; box on
    plot(log10(1./fmt), real(complex_resmt), 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-5.e-5 7.e-5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
    %set(gca, 'Ylim', [-7.e-5 9.e-5])
    xlabel('Period(s)')
    %ylabel('residual(\Omega-m)')
    ylabel('residual(REAL(Z))')

    subplot(p12)
    hold on; box on
    plot(log10(1./fmt), imag(complex_resmt), 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-0.2 3.5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
    %set(gca, 'Ylim', [-1.5e-4 1.5e-4])
    xlabel('Period(s)')
    %ylabel('residual(degree)')
    ylabel('residual(IMAG(Z))')
end
if iswd == 1
    subplot(p13)
    hold on; box on
    plot(Tswd, resswd, 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    %set(gca, 'Ylim', [-0.005  0.01])
    xlabel('Period(s)')
    ylabel('residual(Km/s)')
end
if iell == 1
    subplot(p13)
    hold on; box on
    %plot(Tell, resell, 'xk', 'MarkerSize', 12)
    semilogx(Tell, resell, 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    %set(gca, 'Ylim', [-0.005  0.01])
    xlabel('Period(s)')
    ylabel('residual')
end
if imt == 1
    subplot(p21)
    hold on; box on
    plot(log10(1./fmt), predmt(1:NDAT_MT), '-bo', log10(1./fmt), obsmt(1:NDAT_MT), 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-0.2 3.5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    %set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [-220 420])
    %ylabel('Apparent resistivity (\Omega-m)')
    ylabel('REAL(Z)')
    legend('fitted values', 'obsereved data')

    subplot(p22)
    hold on; box on
    plot(log10(1./fmt), predmt(NDAT_MT+1:2*NDAT_MT), '-bo', log10(1./fmt), obsmt(NDAT_MT+1:2*NDAT_MT), 'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-0.2 3.5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    %set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [-80.0 -55.0])
    %ylabel('Phase(degree)')
    ylabel('IMAG(Z)')
    legend('fitted values', 'obsereved data')
end
if iswd == 1
    subplot(p23)
    hold on; box on
    plot(Tswd, predswd, Tswd, obsswd,'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    set(gca, 'XtickLabel', [])
    set(gca, 'Ylim', [2.5 5.5])
    ylabel('Phase Velocty(Km/s)')
    legend('fitted values', 'obsereved data')
end
if iell == 1
    subplot(p23)
    hold on; box on
    %plot(Tell, predell, Tell, obsell,'xk', 'MarkerSize', 12)
    %plot(Tell, abs(predell), Tell, abs(obsell),'xk', 'MarkerSize', 12)
    %loglog(Tell, abs(predell), Tell, abs(obsell),'xk', 'MarkerSize', 12)
    semilogx(Tell, abs(predell), Tell, abs(obsell),'xk', 'MarkerSize', 12)
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [2.5 5.5])
    ylabel('Ellipticity')
    legend('fitted values', 'obsereved data')
    
end
if irv == -1
    subplot(p31)
    hold on; box on
    plot((0:NTIME-1)*dt, resrf, 'xk')
    set(gca,'Fontsize',8);
    %set(gca, 'Ylim', [-0.1 0.1])
    xlabel('Time(s)')
    ylabel('residual')

    subplot(p41)
    hold on; box on
    plot((0:NTIME-1)*dt, predrf, (0:NTIME-1)*dt, obsrf, 'xk')
    set(gca,'Fontsize',8);
    set(gca, 'XtickLabel', [])
    ylabel('RF')
    legend('fitted values', 'obsereved data')
end



%  if imt == 1
%      figure(200)
%      predZ2 = predmt(1:NDAT_MT).^2 + predmt(NDAT_MT+1:2*NDAT_MT).^2;
%      obsZ2 = obsmt(1:NDAT_MT).^2 + obsmt(NDAT_MT+1:2*NDAT_MT).^2;
%      resZ2 = obsZ2 - predZ2;
%     %subplot(p11)
%     subplot(2,1,1)
%     hold on; box on
%     plot(log10(1./fmt), predZ2, log10(1./fmt), obsZ2, 'xk')
%     set(gca,'Fontsize',8);
%     %set(gca, 'Xlim', [-0.2 3.5])
%     %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
%     %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
%     %set(gca, 'Ylim', [-240 440])
%     %xlabel('Period(s)')
%     %ylabel('residual(\Omega-m)')
%     ylabel('|Z|^2)')
% 
%     %subplot(p12)
%     subplot(2,1,2)
%     hold on; box on
%     plot(log10(1./fmt),resZ2, 'xk')
%     set(gca,'Fontsize',8);
%     %set(gca, 'Xlim', [-0.2 3.5])
%     %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
%     %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
%     %set(gca, 'Ylim', [-220 420])
%     xlabel('Period(s)')
%     %ylabel('residual(degree)')
%     ylabel('residual(|Z|^2)')
% 
%     figure(300)
%      predZ = sqrt(predZ2);
%      obsZ = sqrt(obsZ2);
%      resZ = obsZ - predZ;
%     %subplot(p11)
%     subplot(2,1,1)
%     hold on; box on
%     plot(log10(1./fmt), predZ, log10(1./fmt), obsZ, 'xk')
%     set(gca,'Fontsize',8);
%     %set(gca, 'Xlim', [-0.2 3.5])
%     %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
%     %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
%     %set(gca, 'Ylim', [-240 440])
%     %xlabel('Period(s)')
%     %ylabel('residual(\Omega-m)')
%     ylabel('|Z|)')
% 
%     %subplot(p12)
%     subplot(2,1,2)
%     hold on; box on
%     plot(log10(1./fmt),resZ, 'xk')
%     set(gca,'Fontsize',8);
%     %set(gca, 'Xlim', [-0.2 3.5])
%     %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
%     %set(gca, 'XtickLabel', {'1', '10', '100', '1000', '2000'})
%     %set(gca, 'Ylim', [-220 420])
%     xlabel('Period(s)')
%     %ylabel('residual(degree)')
%     ylabel('residual(|Z|)')
% 
%  end