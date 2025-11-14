function [] = refl_coeff_fft_analy();

load x_jan_1b.mat
NFFT = length(F(1).freq);

for i=1:length(F(1).freq)

    R(i,:)=F(i).dat;

end;
Df = F(1).freq(end)-F(1).freq(1);
Fs = mean(diff(F(1).freq))

%load fft_analy_syn.mat
NFFT = 8192;
%%NFFT = length(freq);
%Df = freq(end)-freq(1)
%Fs = mean(diff(freq))
%R = ref';
%size(R)

for i=80:10:80;

    FFTR = abs(fft(detrend(mean(R(:,i-2:i+2),2)),NFFT));
    [fe,P1]=efft(1/Fs,detrend(mean(R(:,i-2:i+2),2)),0,'r');
    f = Fs/2*linspace(0,1,NFFT/2);
%    f  = 1/Df * (0:nfft/2);
    f_ny = NFFT/(2*Df);

%    figure(i);hold on;
%    title(F(1).ang(i));
%    plot(1./f,FFTR(1:NFFT/2),'k');
    figure;hold on;
    size(P1)
    size(fe)
    plot(1./fe(1:length(P1)/2),P1,'k');
%    set(gca,'XLim',[0 1]);

end;
f(1)

return;
