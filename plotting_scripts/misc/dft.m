function [fr,fi] = dft(t,f)
% Basic unevenly spaced DFT routine.
% To use:  [fr,fi] = dft(t,f);
% where (t,f) represents a digital form
% of the times series f(t).
%
% by D. Holmgren holmgren@phobos.astro.uwo.ca

% prompt for input data...
nobs = length(f);
f0 = input('Starting frequency: ');
f1 = input('End frequency: ');
nfreq = input('Number of frequencies: ');
% remove mean from data...
f = f - mean(f); t = t - mean(t);
% set up arrays...
fr = zeros(nfreq,1);
fi = fr; wr = fr; wi = fi;
dft = fr;
freq = linspace(f0,f1,nfreq)';
twopi = 2 * pi;
% compute data and window amplitude spectra...
for j = 1:nfreq
    x = twopi * freq(j) .* t;
    c = cos(x);
    s = sin(x);
    fr(j) = sum( f .* c);
    fi(j) = sum( f .* s);
    wr(j) = sum(c); wi(j) = sum(s);
end
dft = (2/nobs) .* sqrt( fr .* fr + fi .* fi);
wft = (2/nobs) .* sqrt( wr .* wr + wi .* wi);
% find the maximum...
disp(' Maximum value of Fourier amplitude:')
[dftmax,i] = max(dft); disp(dftmax)
disp(' at the frequency:')
frmax = freq(i); disp(frmax)
% do plot of amplitude spectra...
subplot(211),plot(freq,dft)
title('Data DFT')
% xlabel('Frequency in cycles/unit time')
ylabel('DFT amplitude')
subplot(212),plot(freq,wft)
title('Window DFT')
xlabel('Frequency in cycles/unit time')
ylabel('DFT Amplitude')
% following line is for PC-Matlab...
% shg
