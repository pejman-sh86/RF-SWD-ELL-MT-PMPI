function [cp,alfp,cs,alfs,rho]=GS_Mod(Kg,Kf,rhog,rhof,b,Kp,n,f)
%function [cp,alfp,cs,alfs,rho]=GS_Mod(Kg,Kf,rhog,rhof,b,Kp,n,f)
%
%[cp,alfp,cs,alfs,rho]=GS_Mod(3.6e10,2.374e9,2730,1005,0.377,3.1543e+008,0.0851,f); test of Fig 2;
%[cp,alfp,cs,alfs,rho]=GS_Mod(3.6e10,2.374e9,2730,1005,0.63,gpo=3.3e8,0.0851,f); %test of silty clay
%
% Charles W. Holland September 26, 2011
%
% Buckingham Grain Shearing (GS) 2005 model (in JASA 117, 2005) but *modified* so that 
    %1) the grain-to-grain shear modulus gso is not an input; rather a ratio with the 
       %compressional grain-grain modulus gpo; gs retains dependence on porosity and mean grain size
    %2) mean grain size is not an input, rather it is taken from Eq 18 with DELT=1e-6 which agrees 
        %with Richardson data
    %3) depth dependence relations are removed
% It is the same as Buckingham2005_mod version, except that Kp is an input rather than gpo. This should 
% be more robust for inversion (Kp is independent of b, but using gpo, Kp is not ) 
%
% 
% Inputs:
    % Kg bulk modulus of mineral grains (Pa)
    % Kf bullk modulus of interstitial fluid (Pa)
    % rhog density of grains (kg/m^3)
    % rhof density of fluid (kg/m^3)
    % b porosity
    % Kp - grain-grain compressional modulus (Pa)
    % n  strain hardening index  (0.05 -0.15)
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
 % n=0.0851;  %this could be a parameter
uo=1e-3;
bo=0.377;
%gpo=3.888e8;% gpo - compressional coefficient (Pa)
%gso=4.588e7; % gso - shear coefficient (Pa)

%Bachman does not seem such a good fit to observations
%u_phi= -1.018e1 +3.811e-1*(b*100) -1.85e-3*(b*100).^2; u=2^-u_phi * 1e-3; %to get meters from mm

%Buckingham mean grain size
b_min=0.37;
DELT=1e-6; %this is Richardson value (Hamilton would be 3e-6 which increases sound speed for a given gpo)
B=((1-b)./(1-b_min)).^(1/3);
u=2*DELT*(2*B-1)./(1-B);
   %u=415e-6;%testing purposes

w=2*pi*f;
rho=b*rhof + (1-b)*rhog;
ko=1./(b/Kf + (1-b)/Kg );
co=sqrt(ko/rho);
cf=sqrt(Kf/rhof);           % water sound speed

%moduli 
g_trm=u*(1-b)/(uo*(1-bo));
gpo=Kp./g_trm.^(1/3);
gso=gpo/8.5;                % this is scaled based on Table 2 values
Ks=gso*g_trm.^(2/3);

trm=(1 + (Kp+4/3*Ks)/(rho*co^2)*(i*w*T).^n ).^(-1/2);
cp=co./real(trm);
alfp_Nepers=-w./co.*imag(trm);
alfp=alfp_Nepers*20*log10(exp(1))*1000./f;

cs=sqrt(Ks/rho).*(w*T).^(n/2)./cos(n*pi/4);
alfs_Nepers=w*sqrt(rho./Ks).*(w*T).^(-n/2).*sin(n*pi/4);
alfs=alfs_Nepers*20*log10(exp(1))*1000./f;

if 1==2;colr='g:';
   subplot(221);semilogx(f,cp,colr);hold on;ylabel 'compressional speed (m/s)'
   subplot(222);semilogx(f,alfp,colr);hold on;ylabel('compr. attenuation (dB/m/kHz)')
   subplot(223);semilogx(f,cs,colr);hold on;xlabel('Frequency (Hz)');ylabel 'shear speed (m/s)'
   subplot(224);semilogx(f,alfs,colr);hold on;xlabel('Frequency (Hz)');ylabel('shear attenuation (dB/m/kHz)')
end





    



