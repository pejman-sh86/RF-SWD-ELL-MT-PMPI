function [c_b,alf_b,rho_eff]= Biot_EDFM(Ks,Kf,rhos,rhof,B,n,p,v,f)
%function [c_b,alf_b,rho_eff]= Biot_EDFM(Ks,Kf,rhos,rhof,B,n,p,v,f)
%
% [c_b,alf_b,rho_eff]=Biot_EDFM(3.6e10,2.25e9,2650,1000,0.4,1.25,1e-10,1e-3,f);%test case
%
% Effective Fluid Density Model from Kevin Williams JASA, 110, 2276, 2001
% code follows BUckingham JASA 116, pg 769, 2004

% Charles W. Holland Sept 9, 2011 
%
% INPUTS
 % Ks bulk modulus of mineral grains (Pa)
 % Kf bullk modulus of interstitial fluid (Pa)
 % rhos density of grains (kg/m^3)
 % rhof density of fluid (kg/m^3)
 % B porosity (0-1; for unconsolidated seds generally 0.35 - 0.9)
 % n pore tortuousity  (1-3)
 % p permeability
 % v fluid viscosity (~1e-3 Pa s)
 % f frequency (Hz)
%
% OUTPUTS
 %c_b compressional sound speed (m/s)
 %rho bulk density (g/cm^3)
 %alf_b compressional attenuation (dB/m/kHz)
 
 
 w=2*pi*f;
 rho=rhos*(1-B)+rhof*B; %Compute rho and H from Woods Eqn
 H=((1-B)/Ks + B/Kf)^-1;
 u=p*rhof/v;
   %crit_freq=eta*poros/(2*pi*perm*rhof);
   %porsize=sqrt(8*n*p/B); kap=porsize*sqrt(w.*rhof/v);
 kap=sqrt(8*n/B*w*u);
   barg=kap*sqrt(i);
   tkapd= -sqrt(i)*besselj(1,barg);
   %barg=kap.*exp(i*3*pi/4); %this was from my old code and worked there but maybe different e+iwt
   %tkapd= -exp(i*3*pi/4) * besselj(1,barg);
  %F=(kap./4)./(tkap + i*2./kap); ; %this was from my old code and worked there but maybe different e+iwt
   tkapn= besselj(0,barg);
  tkap=tkapn./tkapd;
  F=(kap./4)./(tkap - i*2./kap);
  %F=1;
 r1=n*(1-B)*rhos + B*(n-1)*rhof;
 r2=B*(1-B)*rhos + (n-2*B + B^2)*rhof; 
 rnum=rhof*(r1 + i*B*rho*F./(w*u));
 rden = r2 + i*B*rhof*F./(w*u);
 c=1./sqrt(rnum./rden/H); 
 
 c_b=real(c);
 alf_b_Nepers=imag(w./c);
 alf_b=alf_b_Nepers*20*log10(exp(1))*1000./f;
 rho_eff=H./c.^2;
