filebase = 'HON';
%%
NRF1 =1; NTIME = 300; sampling_dt = 0.1; sbRF = 3.e-2; scaleRF = 1.e-1; rRF = 0.9; fign = 0;
%% 
NMODE = 1; NDAT_SWD = 50; sbSWD = 2.0e-2; scaleSWD = 1.e-2; rSWD = 0.8; fign = 0;
%%
NMODE_ELL = 1; NDAT_ELL = 30; sbELL = 1.e-3; scaleELL = 1.e-2; rELL = 0.8; fign = 0;
%% 
NDAT_MT = 22; sbMT = 1.0e-4;  scaleMT = 1.e-2; rMT = 0.8; fign= 0;
%% 

DpredRF = importdata(strcat(filebase,'_mappred.dat')); %NRF1xNTIME array
indataRF = importdata(strcat(filebase,'_RF.txt')); %every second row:taper
taper_dpred = indataRF(2,:);
%% 
sdRF = sbRF + scaleRF * abs(DpredRF);
covRF = zeros(NTIME, NTIME);
for l = 1:NTIME
    for m = 1:NTIME
        covRF(l,m) = sdRF(l)*sdRF(m)*(rRF ^ abs(l-m));
    end
end
%%
fign = fign + 1;
figure(fign)
imagesc(covRF)
%%
LcovRF = chol(covRF, 'lower');
%%
noiseRF = LcovRF * randn(size(DpredRF))';
DobsRF = DpredRF + noiseRF';
indataRF(1,:) = DobsRF;
%% 

fign = fign + 1;
figure(fign)
plot([0:NTIME-1]*sampling_dt, DpredRF,[0:NTIME-1]*sampling_dt, DobsRF, '-*k') 
title('RF data')
xlabel('Time(s)')
ylabel('RF')
legend('true model', 'observed data')
%% 
inv_covRF = inv(covRF);
%%
% DpredRF = importdata(strcat(filebase,'_mappred.dat'));
% DobsRF = importdata(strcat(filebase,'_RF.txt'));
% inv_covRF = importdata(strcat(filebase, '_Cdi.dat'));

DresRF = DobsRF - DpredRF;
EtmpR = zeros(NRF1);

for iaz = 1:NRF1
    EtmpR(iaz) = -DresRF(iaz,:)*inv_covRF*DresRF(iaz,:)' / 2.;
end

logLR = sum(EtmpR);
fprintf('logLR = %f\n', logLR);
%%
save(strcat(filebase,'_RF.txt'), 'indataRF', '-ascii');
%%
save(strcat(filebase, '_Cd_true.dat'), 'covRF', '-ascii')
%%
save(strcat(filebase, '_Cdi_true.dat'), 'inv_covRF', '-ascii')
%% 

DpredSWD = importdata(strcat(filebase,'_mappredSWD.dat')); %NMODExNDAT_SWD array
indataSWD = importdata(strcat(filebase,'_SWD.dat')); % first column: periods
periods = indataSWD(:,1);
%% 
sdSWD = sbSWD + scaleSWD * abs(DpredSWD);
covSWD = zeros(NDAT_SWD, NDAT_SWD);
for l = 1:NDAT_SWD
    for m = 1:NDAT_SWD
        covSWD(l,m) = sdSWD(l)*sdSWD(m)*(rSWD ^ abs(l-m));
    end
end
%%
fign = fign + 1;
figure(fign)
imagesc(covSWD)
%%
LcovSWD = chol(covSWD, 'lower');
%%
noiseSWD = LcovSWD * randn(size(DpredSWD))';
DobsSWD = DpredSWD + noiseSWD';
indataSWD(:,2) = DobsSWD;

%% 

fign = fign + 1;
figure(fign)
plot(periods, DpredSWD, periods, DobsSWD, 'x-k')
title('SWD data');
xlabel('period (s)');
ylabel('Phase velocity (Km/s)')
legend('true model', 'observed data')
%% 
inv_covSWD = inv(covSWD);
%%
% DpredSWD = importdata(strcat(filebase,'_mappredSWD.dat'));
% indataSWD = importdata(strcat(filebase,'_SWD.dat'));
% DobsSWD = indataSWD(:,2)';
% inv_covSWD = importdata(strcat(filebase, '_CdiSWD.dat'));

DresSWD = DobsSWD - DpredSWD;
EtmpSWD = zeros(NMODE);

for imod = 1:NMODE
       EtmpSWD(imod) = -DresSWD(imod,:)*inv_covSWD*DresSWD(imod,:)' / 2.;
end

logLSWD = sum(EtmpSWD);
fprintf('logLSWD = %f\n', logLSWD);
%% 
save(strcat(filebase,'_SWD.dat'), 'indataSWD','-ascii');
%%
save(strcat(filebase, '_CdSWD_true.dat'), 'covSWD', '-ascii')
%%
save(strcat(filebase, '_CdiSWD_true.dat'), 'inv_covSWD', '-ascii')
%%
%temporary section for testing gpell outputs for which the first column is
%frequency
gpell_file = 'testmodel3gpell';
%gpell_file = 'testmodel_fail1gpell';
indataELL = importdata(gpell_file); % first column: frequencies
periods_ELL = 1 ./ indataELL(:,1);
indataELL(:,1) = periods_ELL;
indataELL = indataELL(100:5:224, :);
indataELL = indataELL(end:-1:1,:);
save(strcat(filebase,'_ELL.dat'), 'indataELL','-ascii');
%%
DpredELL = importdata(strcat(filebase,'_mappredELL.dat')); %NMODE_ELLxNDAT_ELL array
indataELL = importdata(strcat(filebase,'_ELL.dat')); % first column: periods
periods_ELL = indataELL(:,1);
%DobsELL = indataELL(:,2);
%%
noiseELL = sdELL * randn(size(DpredELL));
DobsELL = DpredELL + noiseELL;
indataELL(:,2) = DobsELL;
%%
save(strcat(filebase,'_ELL.dat'), 'indataELL','-ascii');

%%
fign = fign + 1;
figure(fign)
%plot(periods_ELL, DpredELL, periods_ELL, DobsELL, 'xk')
loglog(periods_ELL, abs(DpredELL), periods_ELL, abs(DobsELL), 'xk')
title('ELL data');
xlabel('period (s)');
ylabel('ell')
legend('true model', 'observed data')
%%

DresELL = log10(abs(DobsELL)) - log10(abs(DpredELL));
%DresELL = abs(DobsELL) - abs(DpredELL);
EtmpELL = zeros(NMODE_ELL);

for imod = 1:NMODE_ELL
       EtmpELL(imod) = log(1.0/(2.0*pi)^(NDAT_ELL/2.0)) ... 
                     -(sum(DresELL(imod,:).^2)./(2.0*sdELL(imod)^2) ...
                     +NDAT_ELL*log(sdELL(imod)));
end

logLELL = sum(EtmpELL);
fprintf('logLELL = %f\n', logLELL);
%% 
DpredMT = importdata(strcat(filebase,'_mappredMT.dat')); %1x2NDAT_MT array
DpredMT_cmplx = complex( DpredMT(1:NDAT_MT),DpredMT(NDAT_MT+1:2*NDAT_MT) );
indataMT = importdata(strcat(filebase,'_MT.dat')); % first column: frequencies
freqMT = indataMT(:,1);
%% 
sdMT = sbMT + scaleMT * abs(DpredMT_cmplx);
covMT = zeros(NDAT_MT, NDAT_MT, 'like', 1.j);
for l = 1:NDAT_MT
    for m = 1:NDAT_MT
        covMT(l,m) = sdMT(l)*sdMT(m)*(rMT ^ abs(l-m)) + 0.j;
    end
end
covMT = complex(covMT);
%%
fign = fign + 1;
figure(fign)
subplot(1,2,1)
imagesc(real(covMT))
subplot(1,2,2)
imagesc(imag(covMT))
%%
LcovMT = chol(covMT, 'lower');
%%
z_standardGauss = randn(size(DpredMT_cmplx))/sqrt(2.) + 1.j*( randn(size(DpredMT_cmplx))/sqrt(2) );
noiseMT_cmplx = complex(LcovMT) * z_standardGauss.';
DobsMT_cmplx = DpredMT_cmplx + noiseMT_cmplx.';

indataMT(:,2) = real(DobsMT_cmplx);
indataMT(:,3) = imag(DobsMT_cmplx);
%% 

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, DpredMT(1:NDAT_MT), 1./freqMT, real(DobsMT_cmplx), '-xk')
title('MT REAL(Z)');
xlabel('period (s)');
ylabel('REAL(Z)')

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, DpredMT(NDAT_MT+1:2*NDAT_MT), 1./freqMT, imag(DobsMT_cmplx), '-xk')
title('MT IMAG(Z)');
xlabel('period (s)');
ylabel('IMAG(Z)');
%%
inv_covMT = inv(covMT);
%% 

% DpredMT = importdata(strcat(filebase,'_mappredMT.dat')); %1x2NDAT_MT array
% DpredMT_cmplx = complex( DpredMT(1:NDAT_MT),DpredMT(NDAT_MT+1:2*NDAT_MT) );
% indataMT = importdata(strcat(filebase,'_MT.dat'));
% DobsMT_cmplx = transpose( complex(indataMT(:,2), indataMT(:,3)) );
% inv_covMT_aug = importdata(strcat(filebase, '_CdiZMT.dat'));
% inv_covMT = complex(inv_covMT_aug(:,1:NDAT_MT), inv_covMT_aug(:, NDAT_MT+1:2*NDAT_MT));


DresMT_cmplx = transpose(DobsMT_cmplx - DpredMT_cmplx);

logLMT = -DresMT_cmplx' * inv_covMT * DresMT_cmplx;

fprintf('logLMT = %f\n', logLMT);
%% 
logL = logLSWD + logLMT + logLR;

fprintf('logL = %f\n', logL);
%% 
save(strcat(filebase,'_MT.dat'), 'indataMT','-ascii');
%%
covMT_aug = zeros(NDAT_MT, 2*NDAT_MT);
covMT_aug(:, 1:NDAT_MT) = real(covMT);
covMT_aug(:, NDAT_MT+1:2*NDAT_MT) = imag(covMT);
save(strcat(filebase, '_CdZMT_true.dat'), 'covMT_aug', '-ascii')
%%
inv_covMT_aug = zeros(NDAT_MT, 2*NDAT_MT);
inv_covMT_aug(:, 1:NDAT_MT) = real(inv_covMT);
inv_covMT_aug(:, NDAT_MT+1:2*NDAT_MT) = imag(inv_covMT);
save(strcat(filebase, '_CdiZMT_true.dat'), 'inv_covMT_aug', '-ascii')
%% 

save simdata10.mat
