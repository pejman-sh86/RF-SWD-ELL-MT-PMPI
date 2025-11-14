function [cp,alfp,cs,alfs,rho]=VGSlambda_no_tau_s(Kg,Kf,rhos,rhof,b,Kp,Ks,n,tau,f)
% [cp,alfp,cs,alfs,rho]=VGSlambda_no_tau_s(Kg,Kf,rhos,rhof,b,Kp,Ks,n,tau,f)    
%[cp,alfp,cs,alfs,rho]=VGSlambda_no_tau_s(3.2e10,2.395e9,2690,1023,0.385,3.796e8,1.711E7,0.0617,1.2E-4,f); %test of Fig 1 
%
% Buckingham 2000 model (in JASA 127, 2010), where instead of requiring tau_s   
    % as an input, it is calculated directly.  Original March 20, 2013, corrected
    % April 4, 2014
%
% --------- Charles W. Holland April 4, 2014 ------------
%
% Inputs:
    % Kg bulk modulus of mineral grains (Pa)
    % Kf bullk modulus of interstitial fluid (Pa)
    % rhos density of grains (kg/m^3)
    % rhof density of fluid (kg/m^3)
    % b porosity
    % Kp - grain-grain compressional modulus (Pa)
    % Ks - grain-grain shear modulus (Pa)
    % n  strain hardening index  (0.05 -0.15)
    % tau - compressional wave viscoelastic time constant (s)
    % f frequency (Hz); can be vector
%
% Outputs:
    %cp - compressional velocity (m/s)
    %alfp - compressional attenuation (dB/m/kHz)
    %cs - shear velocity (m/s)
    %alfs - shear attenuation (dB/m/kHz)
    %rho bulk density (kg/m^3)


%constants
T=1;

w=2*pi*f; 
rho=b*rhof + (1-b)*rhos;
ko=1./(b/Kf + (1-b)/Kg );
co=sqrt(ko/rho);
cf=sqrt(Kf/rhof); %water sound speed

gw=(1 + 1./(1i*w*tau)).^(-1+n);
trm=(1 + gw.*(Kp+4/3*Ks)/(rho*co^2).*(1i*w*T).^n ).^(-1/2);
cp=co./real(trm);
alfp_Nepers=-w./co.*imag(trm);
alfp=alfp_Nepers*20*log10(exp(1))*1000./f;

%compute tau_s 
wT=1/tau; %threshold frequency for p-waves
gwT=(1 + 1./(1i*wT*tau)).^(-1+n);
trmT=(1 + gwT.*(Kp+4/3*Ks)/(rho*co^2).*(1i*wT*T).^n ).^(-1/2);
cpT=co./real(trmT); %sound speed at threshold frequency fTp=1/(2*pi*tau_p)

ts_num=cpT*cos((n+1)*pi/8)*tau;
ts_den=sqrt(Ks/rho)*2^((n-1)/4); 
tau_s=(ts_num./ts_den)^(2/(2-n));

gws=(1 + 1./(1i*w*tau_s)).^(-1+n);
trm_s=((1i*w*T).^n.*gws).^(-1/2);
cs=sqrt(Ks/rho)./real(trm_s); 
alfs_Nepers=-w*sqrt(rho./Ks).*imag(trm_s);
alfs=alfs_Nepers*20*log10(exp(1))*1000./f;

%figure
if 1==2;colr='c';
   subplot(221);semilogx(f/1e3,cp/cf,colr);hold on;ylabel 'sound speed ratio';axis([0.1 1e3 1.02 1.18]);grid
   subplot(223);loglog(f/1e3,alfp.*f/1000,colr);hold on;xlabel('Frequency (kHz)');ylabel('compr. attenuation (dB/m)');axis([0.1 1e3 1e-2 1e3]);grid
   subplot(222);semilogx(f/1e3,cs,colr);hold on;ylabel 'shear speed (m/s)';axis([0.1 1e3 0 180]);grid
   subplot(224);loglog(f/1e3,alfs.*f/1000,colr);hold on;xlabel('Frequency (kHz)');ylabel('shear attenuation (dB/m)');axis([0.1 1e3 1e0 1e4]);grid
end
