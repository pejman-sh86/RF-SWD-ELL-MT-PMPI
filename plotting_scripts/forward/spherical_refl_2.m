function [Rs, Rp] = spherical_refl_2(geo,z_t,f,Th)
% function [Rs, Rp] = spherical_refl(geo,z_t,f,Th)
%          [Rs, Rp] = spherical_refl(geo,60,[500:500:2000],5:60)
%
% Ouputs
%----------------------------------------
%    Rs - Spherical Reflection Coefficient (complex)
%    Rp - Plane wave Reflection Coefficient (complex)
%
% Inputs
%--------------------------------------------------------------------------
%
% "geo" contains environmental data, e.g., 
%   thickness    speed    attenuation    density
%     (m)        (m/s)    (dB/m/kHz)      (g/cc)
%
% geo=[NaN        1511          0          1.03
%     1.7         1500        0.015        1.45
%     2.8         1555         0.1         1.7 
%     0.6         1600         0.1         2.0
%     NaN         1565         0.01        1.9];
%
% z_t is the sum of the source and receiver heights,  
% f is frequency (Hz): can be a vector
% Th is grazing angle (degrees): can be a vector
%
% for smoothing in frequency
% nsmooth=5;RU=sqrt(conv(abs(Rs).^2,ones(1,nsmooth))/nsmooth  );
%
% Notes: code written by Charles Holland and John Camin
% ---------------------------------------------------
%   -uses Simpsons rule for numerical integration
%   -analytic expression for N (integrand discretization)
%   -auto-scaling for umax; exponential term along imag axis falls to epsilon=1e-4
%     k*z_t*cos(pi/2-i*theta)=-ln(epsilon), sinh(theta)=-ln(1e-4)/k*z_t;
%   -arbitrary layering (fluid) & vectorization (~30% faster)
%   -used large arg to Bessel Function ~(100% faster)
%
% could speed up code perhaps another factor of 2 by integrating across all angles at once
% i.e., time is wasted sampling along the imaginary axis where generally the solution is not oscillatory
%
%Note 2: N & umax determine how finely the complex plane is sampled.  Since
%this is a 'brute force' calculation, this is directly related to the
%accuracy of the solution.  The path goes from 0 to +Pi and from 
%+Pi to +Pi - umax*j. 
%--------------------------------------------------------------------------

% plane wave reflection coefficient over sampled angles
Rp = zeros(length(f),length(Th));
Rp = ref_nlay3c(Th, geo,f).';

%Environment & transducer Constants...
c = geo(1,2);

%Spatial Variables
x = z_t./tan(Th*pi/180);   %offset
r = x;                     %xy plane, 2D range
R1 = sqrt(x.^2 + (z_t)^2); %Radial distance of Reflected Path

for kf=1:length(f)
    
   w = 2*pi*f(kf);
   k = w/c;
   RPs = exp(i*k*R1)./R1;
 
   % Set limit for integration along imag axis where exponential goes to 
   % 1e-4 of its value at pi/2;
   exp_arg=9.2/(k*z_t);
   % umax=1.5*asinh(exp_arg), %1.5*angle where exponential goes to 1e-4 
   % of its value at pi/2;
   umax=1.1*asinh(exp_arg); %changed to factor 1.1 on June 28, 2005

   %Fix N from analytic considerations, 5 samples per num oscillations 0->pi/2
   Nr_fr=floor(k*max(r+z_t)/(2*pi));
   Nr=max(Nr_fr*5,80);
   if mod(Nr,2)==0;Nr=Nr+1;end;

   % Theta Spacing in Complex Plane...
   % int_end= pi/2-.001i;
   % int_end= pi/2;rTh = linspace(0,int_end,N);
   % iTh = linspace(int_end,pi/2-i*umax,N);
   % drTh = mean(diff(rTh));diTh = mean(diff(iTh));

   drTh=pi/2/(Nr-1);
   diTh=-i*drTh;
   rTh=[0:drTh:pi/2  pi/2-i*(drTh:drTh:umax)]+1e-100; % the 1e-100 eliminates 
                                                      % warnings for divide by 0
   % force the total to also be odd then the int along imag will also 
   % be odd N-Nr+1
   if mod(length(rTh),2)==0; rTh=rTh(1:end-1);end ; 
   N=length(rTh);

   %Compute the reflection coef over the complex integration path
   Rr = ref_nlay3c(90-rTh*180/pi, geo,f(kf));

   %Numerical Integration masking scheme
   m = ones(1,Nr);      %Masking array for Simpsons Rule
   m(2:Nr-1) = 4;       %(1,4,2,4,2,4,1)
   m(3:2:Nr-1) = 2;

   mi = ones(1,N-Nr+1);      %Masking array for Simpsons Rule
   mi(2:N-Nr-1+1) = 4;       %(1,4,2,4,2,4,1)
   mi(3:2:N-Nr-1+1) = 2;
   exp_termR=exp(i.*k.*(z_t).*cos(rTh)).*sin(rTh);
   tmp = m.*exp_termR(1:Nr);
   tmpi = mi.*exp_termR(Nr:N);
   kr_sinRTh=k*r'*sin(rTh);

   % Spherical Reflection Coefficent RP
   btR=sqrt(2./(pi*kr_sinRTh)).*cos(kr_sinRTh-1/4*pi); %large arg approx >4;
   jj=find(kr_sinRTh<4);  btR(jj)=besselj(0,kr_sinRTh(jj) );

   % One step simpsons rule integration REAL
   rRP = drTh/3*btR(:,1:Nr)*(Rr(1:Nr).*tmp).';

%   rRP
%   drTh
%   tmp
%   Rr(1:Nr)
%   Nr

   % One step simpsons rule integration IMAG
   iRP = diTh/3*btR(:,Nr:N)*(Rr(Nr:N).*tmpi).';

   Rs(:,kf) = i*k*(rRP+iRP)./ RPs.';
%   save bla
end;
return;
