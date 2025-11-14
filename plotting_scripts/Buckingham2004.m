function [c_p,rho,alf_p]=Buckingham2004(Ks,Kf,rhos,rhof,B,n,fT,f)
%function [c_p,rho,alf_p]=Buckingham2004(Ks,Kf,rhos,rhof,B,n,fT,f)
%
% [c_p,rho,alf_p]=Buckingham2004(3.6e10,2.25e9,2650,1000,0.4,1.25,1e3,f);test case
%
% Buckingham 'fluid' model based on JASA, 116, 769 - 776, 2004 
%
% Charles W. Holland Sept 9, 2011
%
% INPUTS
    % Ks bulk modulus of mineral grains (Pa)
    % Kf bullk modulus of interstitial fluid (Pa)
    % rhos density of grains (kg/m^3)
    % rhof density of fluid (kg/m^3)
    % B porosity (0-1; for unconsolidated seds generally 0.35 - 0.9)
    % n pore tortuosity  (1-3)
    % fT transition frequency (Hz)
    % f frequency (Hz)
%
% OUPUTS
    %c_p compressional wave speed (m/s)
    %rho bulk density (kg/m^3)
    %alf_p compressional wave attenuation (dB/m/kHz)


%First compute co and rho from Woods Eqn
 rho=rhos*(1-B)+rhof*B;
 co=sqrt(Ks*Kf./( ((1-B)*Kf + B*Ks).*rho ));%low freq sound speed limit co
  
 u=1/(2*pi*fT);
 trm = 1 + (2*pi*f*u).^2;
 
 r1=n*(1-B)*rhos + B*(n-1)*rhof;
 r2=B*(1-B)*rhos + (n-2*B + B^2)*rhof;
 cI=co*sqrt(r2*rho/(r1*rhof)); %cI infinite frequency sound speed (m/s)
 
 c_p=1./( (1/co-1/cI)/sqrt(2)*sqrt( (sqrt(trm)+1)./trm) + 1/cI );
 alf_p_Nepers=2*pi*f/sqrt(2)*(1/co-1/cI).*sqrt( (sqrt(trm)-1)./trm) ;
 alf_p=alf_p_Nepers*20*log10(exp(1))*1000./f;

return; 
