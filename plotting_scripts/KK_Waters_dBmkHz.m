function [alf1dB,vp]=KK_Waters_dBmkHz(n,fr,alfrdB,cr,f1)
% vp=KK_Waters(n,fr,alfr,cr,f1)
%
% %Kramers Kronig dispersion from Waters et al. (JASA 108, pg 556)
%
% INPUTS
% n = frequency exponent of attenuation assuming power law
% fr = reference frequency for attenuation
% alfr = attenuation in dB/m/kHz at reference frequency fr
% cr = sound speed at reference frequency, fr 
% f1 = frequency of interest (can be a vector)
%
% OUTPUTS
% alf1dB = attenuation in (dB/m/kHz) for frequency vector f1
% vp = compressional velocity in (m/s)


alfr = alfrdB*(fr/1000)^n/(20*log10(exp(1))) ; %reference attenuation in nepers/m
wo=2*pi*fr;w1=2*pi*f1;
%Q1 = quad(@w_int, w1, w2,[],[],n)

if n==1
 Q2 = -2/pi*log(f1/fr);
else
  p=n-1;  
 Q2 = (w1.^p-wo.^p).*tan(n*pi/2) ;
end

if n > 3; 'n must be less than 3!',  vp=[],alfp=[]
 else
    vp=1./(1./cr + Q2.*alfr./(2*pi*fr)^n );
    alf1dB= alfrdB*(f1/fr).^(n-1);
end

return;
