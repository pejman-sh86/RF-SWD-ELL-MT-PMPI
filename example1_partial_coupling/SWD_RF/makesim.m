filebase = 'HON';
%%
NRF1 =1; NTIME = 300; sampling_dt = 0.1; sdRF = 5.e-2; fign = 0;
%% 
NMODE = 1; NDAT_SWD = 50; sdSWD = 2.0e-2; fign = 0;
%%
NMODE_ELL = 1; NDAT_ELL = 30; sdELL = 1.e-3; fign = 0;
%% 
NDAT_MT = 22; sdMT = 1.0e-4; fign= 0;
%% 

DpredRF = importdata(strcat(filebase,'_mappred.dat')); %NRF1xNTIME array
indataRF = importdata(strcat(filebase,'_RF.txt')); %every second row:taper
taper_dpred = indataRF(2,:);
%% 

noiseRF = sdRF * randn(size(DpredRF));
DobsRF = DpredRF + noiseRF;
indataRF(1,:) = DobsRF;
%% 
save(strcat(filebase,'_RF.txt'), 'indataRF', '-ascii');
%% 

fign = fign + 1;
figure(fign)
plot([0:NTIME-1]*sampling_dt, DpredRF,[0:NTIME-1]*sampling_dt, DobsRF, 'xk') 
title('RF data')
xlabel('Time(s)')
ylabel('RF')
legend('true model', 'observed data')
%% 

DresRF = DobsRF - DpredRF;
EtmpR = zeros(NRF1);

for iaz = 1:NRF1
       EtmpR(iaz) = log(1.0/(2.0*pi)^(NTIME/2.0)) ... 
                     -(sum(DresRF(iaz,:).^2)./(2.0*sdRF(iaz)^2) ...
                     +NTIME*log(sdRF(iaz)));
end

logLR = sum(EtmpR);
fprintf('logLR = %f\n', logLR);
%% 

DpredSWD = importdata(strcat(filebase,'_mappredSWD.dat')); %NMODExNDAT_SWD array
indataSWD = importdata(strcat(filebase,'_SWD.dat')); % first column: periods
periods = indataSWD(:,1);
%% 

noiseSWD = sdSWD * randn(size(DpredSWD));
DobsSWD = DpredSWD + noiseSWD;
indataSWD(:,2) = DobsSWD;
%% 
save(strcat(filebase,'_SWD.dat'), 'indataSWD','-ascii');

%% 

fign = fign + 1;
figure(fign)
plot(periods, DpredSWD, periods, DobsSWD, 'xk')
title('SWD data');
xlabel('period (s)');
ylabel('Phase velocity (Km/s)')
legend('true model', 'observed data')
%% 
DresSWD = DobsSWD - DpredSWD;
EtmpSWD = zeros(NMODE);

for imod = 1:NMODE
       EtmpSWD(imod) = log(1.0/(2.0*pi)^(NDAT_SWD/2.0)) ... 
                     -(sum(DresSWD(imod,:).^2)./(2.0*sdSWD(imod)^2) ...
                     +NDAT_SWD*log(sdSWD(imod)));
end

logLSWD = sum(EtmpSWD);
fprintf('logLSWD = %f\n', logLSWD);
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
%cmplxDpredMT = complex( DpredMT(1:NDAT_MT),DpredMT(NDAT_MT+1:2*NDAT_MT) );
indataMT = importdata(strcat(filebase,'_MT.dat')); % first column: frequencies
freqMT = indataMT(:,1);
%% 

noiseMT = (sdMT/sqrt(2.)) .* randn(size(DpredMT));
%cmplxDobsMT = cmplxDpredMT + noiseMT;
DobsMT = DpredMT + noiseMT;
prealZ = DpredMT(1:NDAT_MT);
pimagZ = DpredMT(NDAT_MT+1:2*NDAT_MT);
nrealZ = DobsMT(1:NDAT_MT);
nimagZ = DobsMT(NDAT_MT+1:2*NDAT_MT);
indataMT(:,2) = nrealZ;
indataMT(:,3) = nimagZ;
%% 
save(strcat(filebase,'_MT.dat'), 'indataMT','-ascii');
%% 

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, prealZ, 1./freqMT, nrealZ, 'xk')
title('MT REAL(Z)');
xlabel('period (s)');
ylabel('REAL(Z)')

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, pimagZ, 1./freqMT, nimagZ, 'xk')
title('MT IMAG(Z)');
xlabel('period (s)');
ylabel('IMAG(Z)');
%% 

%noiseMT = sdMT*abs(DpredMT) .* randn(size(DpredMT));
noiseMT = sdMT .* randn(size(DpredMT));
DobsMT = DpredMT + noiseMT;
pappRes = DpredMT(1:NDAT_MT);
pphase = DpredMT(NDAT_MT+1:2*NDAT_MT);
nappRes = DobsMT(1:NDAT_MT);
nphase = DobsMT(NDAT_MT+1:2*NDAT_MT);
indataMT(:,2) = nappRes;
indataMT(:,3) = nphase;
%% 
save(strcat(filebase,'_MT.dat'), 'indataMT','-ascii');
%% 

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, pappRes, 1./freqMT, nappRes, 'xk')
title('MT apparent resistivity');
xlabel('period (s)');
ylabel('apparent resistivity (\sigma-m)')

fign = fign + 1;
figure(fign)
semilogx(1./freqMT, pphase, 1./freqMT, nphase, 'xk')
title('MT phase');
xlabel('period (s)');
ylabel('phase (degree)');
%% 
DresMT = DobsMT - DpredMT;

logLMT = -(2*NDAT_MT/2.)*log(2.*pi)- sum(DresMT.^2 ./(2.*(DobsMT*sdMT).^2.)) ...
         -2.*NDAT_MT*log(sdMT) - sum(log(abs(DobsMT)));

fprintf('logLMT = %f\n', logLMT);
%% 
logL = logLSWD + logLMT + logLR;

fprintf('logL = %f\n', logL);
%% 

save simdata8.mat
