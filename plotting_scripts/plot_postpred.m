
filebase = 'HON';
preds = importdata(strcat(filebase,'_postpred.dat'));
parfile = strcat(filebase,'_parameter.dat');

iwresRF = 0; %weighted residual
imake_CdRF = 0;

iwresSWD = 0;
imake_CdSWD = 0;

iwresELL = 0;
imake_CdELL = 0;

iwresMT = 0;
imake_CdMT = 0;

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

if irv == -1
    datarf = importdata(strcat(filebase, '_RF.txt')); % 2xNTIME array
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
        resrf = transpose (inv(Lrf) * resrf' );
    end
end

if iswd == 1
    dataswd = importdata(strcat(filebase, '_SWD.dat')); % NDAT_SWDx2 array
    obsswd = dataswd(:,2); % NDAT_SWDx1 array
    Tswd = dataswd(:,1); % NDAT_SWDx1 array
    resswd = obsswd' - predswd;
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

%% residuals and predicted values

for is = 1:size(preds,1)

    if imt == 1
        subplot(p11)
        hold on; box on
        plot(log10(1./fmt), real(complex_resmt(is,:)), 'Color', [.5 .5 .5])
    
        subplot(p12)
        hold on; box on
        plot(log10(1./fmt), imag(complex_resmt(is,:)), 'Color', [.5 .5 .5])
    end

    if iswd == 1
        subplot(p13)
        hold on; box on
        plot(Tswd, resswd(is,:), 'Color', [.5 .5 .5])
    end

    if iell == 1
        subplot(p13)
        hold on; box on
        %plot(Tell, resell, 'xk', 'Color', [.5 .5 .5])
        semilogx(Tell, resell(is,:), 'Color', [.5 .5 .5])
    end

    if imt == 1
        subplot(p21)
        hold on; box on
        plot(log10(1./fmt), predmt(is,1:NDAT_MT), 'Color', [.5 .5 .5])
    
        subplot(p22)
        hold on; box on
        plot(log10(1./fmt), predmt(is, NDAT_MT+1:2*NDAT_MT), 'Color', [.5 .5 .5])
    end

    if iswd == 1
        subplot(p23)
        hold on; box on
        plot(Tswd, predswd(is,:), 'Color', [.5 .5 .5])
    end

    if iell == 1
        subplot(p23)
        hold on; box onTell, 
        %plot(Tell, predell, 'Color', [.5 .5 .5])
        %plot(Tell, abs(predell), 'Color', [.5 .5 .5])
        %loglog(Tell, abs(predell), 'Color', [.5 .5 .5])
        semilogx(Tell, abs(predell(is,:)), 'Color', [.5 .5 .5])
    end

    if irv == -1
        subplot(p31)
        hold on; box on
        plot((0:NTIME-1)*dt, resrf(is,:), 'Color', [.5 .5 .5])
    
        subplot(p41)
        hold on; box on
        plot((0:NTIME-1)*dt, predrf(is,:),'Color', [.5 .5 .5])

    end

end   %is loop on samples

%% add observed values 

if imt == 1
    subplot(p21)
    hold on; box on
    plot(log10(1./fmt), obsmt(1:NDAT_MT), 'r-')
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-0.2 3.5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    %set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [-220 420])
    %ylabel('Apparent resistivity (\Omega-m)')
    ylabel('REAL(Z)')
    %legend('fitted values', 'obsereved data')

    subplot(p22)
    hold on; box on
    plot(log10(1./fmt), obsmt(NDAT_MT+1:2*NDAT_MT), 'r-')
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [-0.2 3.5])
    %set(gca, 'Xtick', [0, 1, 2, 3 , log10(2000)])
    set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [-80.0 -55.0])
    %ylabel('Phase(degree)')
    ylabel('IMAG(Z)')
    %legend('fitted values', 'obsereved data')
end

if iswd == 1
    subplot(p23)
    hold on; box on
    plot(Tswd, obsswd,'r-')
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [2.5 5.5])
    ylabel('Phase Velocty(Km/s)')
    %legend('fitted values', 'obsereved data')
end

if iell == 1
    subplot(p23)
    hold on; box on
    %plot(Tell, predell, Tell, obsell,'r-')
    %plot(Tell, abs(predell), Tell, abs(obsell),'r-')
    %loglog(Tell, abs(predell), Tell, abs(obsell),'r-')
    semilogx(Tell, abs(obsell),'r-')
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    set(gca, 'XtickLabel', [])
    %set(gca, 'Ylim', [2.5 5.5])
    ylabel('Ellipticity')
    %legend('fitted values', 'obsereved data')
    
end

if irv == -1
    subplot(p41)
    hold on; box on
    plot((0:NTIME-1)*dt, obsrf, 'r-')
    set(gca,'Fontsize',8);
    set(gca, 'XtickLabel', [])
    ylabel('RF')
    %legend('fitted values', 'obsereved data')
end

%% add line of residual zero to residual plots

if imt == 1
    subplot(p11)
    hold on; box on
    plot(log10(1./fmt), zeros(NDAT_MT),'r-')
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
    plot(log10(1./fmt), zeros(NDAT_MT),'r')
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
    plot(Tswd, zeros(NDAT_SWD), 'r-')
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
    semilogx(Tell, zeros(NDAT_ELL), 'r-')
    set(gca,'Fontsize',18);
    %set(gca, 'Xlim', [0.0 300])
    %set(gca, 'Ylim', [-0.005  0.01])
    xlabel('Period(s)')
    ylabel('residual')
end

if irv == -1
    subplot(p31)
    hold on; box on
    plot((0:NTIME-1)*dt, zeros(NTIME), 'r-')
    set(gca,'Fontsize',8);
    %set(gca, 'Ylim', [-0.1 0.1])
    xlabel('Time(s)')
    ylabel('residual')
end

