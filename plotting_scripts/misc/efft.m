 function [f,P1]=efft(fs,ts,sens,colr)
%  [f,P1]=efft(fs,ts,sens,colr)
%[f,P1]=efft(x2.fout,x2.t_s(36,1:15000),-160.5,'r');
%
% plots the ENERGY SPECTRAL DENSITY for an arbitrary matrix
% fs=sample frequency (Hz) 
% ts=time series (column oriented)
% sensitivity "sens", color "colr" for plot

% find n=FFT size power of 2
sz=length(ts);
nf=14;
while sz>2^nf, nf=nf+1;end
n=2^nf


fft1=fft(ts,n);
P1=2*fft1.*conj(fft1)/(fs*fs);
f=fs*(0:(n/2)-1)/n;

figure; orient tall; box on;
  subplot(2,1,1); hold on
 plot((1:sz)/fs*1000,ts,colr);
 xlabel('Time (ms)')
 ylabel('Amplitude (V)')
 grid
 
%subplot(2,1,1); hold on
% plot(f,P1(1:n/2,:),colr);
% xlabel('Frequency (Hz)');ylabel('Energy Spectral Density (uPa^2 s/Hz)')
% grid; zoom

  subplot(2,1,2); hold on;
 plot(f/1000,10*log10(P1(:,1:n/2))-sens,colr);
 xlabel('Frequency (kHz)')
 ylabel('Energy Spectral Density (dB re 1 uPa^2 s/Hz)')
 grid
 zoom

%[xav] = sliding_win_dB(P1',0,25);
