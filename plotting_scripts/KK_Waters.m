function vp=KK_Waters(n,fr,alfr,cr,f1)
% vp=KK_Waters(n,fr,alfr,cr,f1)
%
% %Kramers Kronig dispersion from Waters et al. (JASA 108, pg 556)
%
% n = freqeuncy exponent
% fr = reference frequency for attenuation
% alfr = attenuation in nepers/length at reference frequency fr
% cr = sound speed at reference frequency, fr 
% f1 = frequency of interest (can be a vector)

% Evans/Carey1 KK_Waters(1.5,50,3.03e-4,1560,100:100:6500)
% Evans/Carey2 KK_Waters(1.5,50,4e-4,1585,60,1400)
% O'Donnell KK_Waters(1.05,1e6,70,2000,1e6,(1:10)*1e6)
% HOrton KK_Waters(1.37,12e3,1/8.686,4566,12e3,200e3)
% Hamilton: [k,vp]=hamilton_k(36); vp_grav=KK_Waters(1,1000,k/8.686,vp,100,[otofc 40e3:20e3:1e7]);

wo=2*pi*fr;w1=2*pi*f1;
%Q1 = quad(@w_int, w1, w2,[],[],n)

if n==1
 Q2 = -2/pi*log(f1/fr);
else
  p=n-1;  
 Q2 = (w1.^p-wo.^p).*tan(n*pi/2) ;
end

if n > 3; 'n must be less than 3!',  vp=[]
 else
    vp=1./(1./cr + Q2.*alfr./(2*pi*fr)^n );
end
