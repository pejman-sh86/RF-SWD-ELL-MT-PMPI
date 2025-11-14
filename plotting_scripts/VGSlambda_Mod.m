function [cp,alfp,cs,alfs,rho]=VGSlambda_Mod(Kg,Kf,rhos,rhof,b,Kp,n,tau,f)
%function [cp,alfp,cs,alfs,rho]=VGSlambda_Mod(Kg,Kf,rhos,rhof,b,Kp,n,tau,f)   
%
%[cp,alfp,cs,alfs,rho]=VGSlambda_Mod(3.2e10,2.395e9,2690,1023,0.385,3.796e8,0.0617,1.2E-4,f); test of Fig 1  
%
% Charles W. Holland October 12, 2011
%
% Buckingham 2010 model (in JASA 127, 2010) but *modified* so that 
    %1) the grain-to-grain shear modulus gso is not an input; rather a ratio with the 
       %compressional grain-grain modulus gpo; gs retains dependence on
       %porosity and mean grain size
%       
% Inputs:
    % Kg bulk modulus of mineral grains (Pa)
    % Kf bullk modulus of interstitial fluid (Pa)
    % rhos density of grains (kg/m^3)
    % rhof density of fluid (kg/m^3)
    % b porosity
    % Kp - grain-grain compressional modulus (Pa)
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

uo=1e-3;
bo=0.377;
%Buckingham mean grain size
b_min=0.37;
DELT=1e-6; %this is Richardson value (Hamilton would be 3e-6 which increases sound speed for a given gpo)
B=((1-b)./(1-b_min)).^(1/3);
u=2*DELT*(2*B-1)./(1-B);
%u=415e-6;%testing purposes
%moduli
g_trm=u*(1-b)/(uo*(1-bo) );
gpo=Kp./g_trm.^(1/3);
gso=gpo/10; %rough estimtae
Ks=gso*g_trm.^(2/3);   

w=2*pi*f;
rho=b*rhof + (1-b)*rhos;
ko=1./(b/Kf + (1-b)/Kg );
co=sqrt(ko/rho);
cf=sqrt(Kf/rhof); %water sound speed

gw=(1 + 1./(i*w*tau)).^(-1+n);
trm=(1 + gw.*(Kp+4/3*Ks)/(rho*co^2).*(i*w*T).^n ).^(-1/2);
cp=co./real(trm);
alfp_Nepers=-w./co.*imag(trm);
alfp=alfp_Nepers*20*log10(exp(1))*1000./f;

tau_s=cp(1)/sqrt(Ks/rho)*tau;  %this is a rough approximation (no freq dependence in tau_s); only affects shear atten
gws=(1 + 1./(i*w*tau_s)).^(-1+n);
trm_s=((i*w*T).^n.*gws).^(-1/2);
cs=sqrt(Ks/rho)./real(trm_s);
alfs_Nepers=-w*sqrt(rho./Ks).*imag(trm_s);
alfs=alfs_Nepers*20*log10(exp(1))*1000./f;

%if 1==2;colr='b';
%   subplot(221);semilogx(f/1e3,cp/cf,colr);hold on;ylabel 'sound speed ratio';axis([0.1 1e3 1.02 1.18]);grid
%   subplot(223);loglog(f/1e3,alfp.*f/1000,colr);hold on;xlabel('Frequency (kHz)');ylabel('compr. attenuation (dB/m)');axis([0.1 1e3 1e-2 1e3]);grid
%   subplot(222);semilogx(f/1e3,cs,colr);hold on;ylabel 'shear speed (m/s)';axis([0.1 1e3 0 180]);grid
%   subplot(224);loglog(f/1e3,alfs.*f/1000,colr);hold on;xlabel('Frequency (kHz)');ylabel('shear attenuation (dB/m)');axis([0.1 1e3 1e0 1e4]);grid
%end

