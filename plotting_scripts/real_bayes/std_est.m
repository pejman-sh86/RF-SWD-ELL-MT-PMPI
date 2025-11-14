% ------------------------------------------------------------------------
%                     estimate std dev
% ------------------------------------------------------------------------
%
%    Record:
%       Date            by            last  change
%     ========    ==============      ============
%     10/23/03     jand@uvic.ca         11/26/03
%
% ------------------------------------------------------------------------
%    Description:  
%                 ASSA by Stan Dosso. Forward model by Chrles Holland
%     input:
%    output: 
% ------------------------------------------------------------------------ 
function std_est;

format compact

%---------------------- GLOBAL VARIABLES ----------------------------------
global dat datm freq nfreq ang nang
global c1 rho1 znorm lay_thick sz sd i_nan
%--------------------------------------------------------------------------

npar = 7;
load blde4_2_ping_ida_inv.mat;
aincr   =  2;		%
a_start =  1;		%
nang    = 68;		% Nr of angles
fincr   =  1;		%
f_start =  6;		%
nfreq   =  9; 	% Nr of frequencies
%
% Forward model parameters:
%
c1 = 1511; rho1 = 1.029;
lay_thick = 0.1;    % discretization for layer thickness - controls 
                      % upper frequency limit
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
sz=length(znorm);
ndat = nfreq*nang;
%
% Model:
% h  rhot  rhob  nu  ct  cb  alpha
%
% Synthetic model parameters: mbest
%
mbest = [2.0930214 1.3508944 1.4841878 0.73839097 1472.8123 1467.5688 0.25152184];

%---------------------------------------------------------------------------
%  Load synthetic or real data.
%---------------------------------------------------------------------------

dat = zeros(ndat,1);
%
% REAL DATA
%
disp('I work on real data now!');
ang = zeros(1,nang);
freq = zeros(1,nfreq);

j = a_start;
for i = 1:nang
  ang(i) = xde4.ang(j);
  j = j + aincr;
end
j = f_start;
for i = 1:nfreq
  freq(i) = xde4.pref(j,1);
  j = j + fincr;
end
tmp = xde4.bl(f_start:fincr:f_start+(fincr*nfreq)-1,a_start:aincr:...
              a_start+(aincr*nang)-1);
dat = reshape(tmp',ndat,1);
ndat = zeros(nfreq,1);
for i = 1:nfreq
  ndat(i) = nang - length(find(isnan(dat(nang*(i-1)+1:nang*i))==1));
end
[i_nan] = find(isnan(dat)==1);
%dat = 10.^(-dat/20);
dat(i_nan) = 0;
ang
freq
[datm] = forward(mbest,freq,nfreq,ang,nang);

sd = zeros(nfreq,1);
for i = 1:nfreq
  sd(i) = sqrt(1/(ndat(i)-npar) * sum((datm(nang*(i-1)+1:nang*i)-...
          dat(nang*(i-1)+1:nang*i)).^2));
end

stdv_dB = [freq' sd];

save ml_std.mat stdv_dB

%===========================================================================
%   FORWARD.M
%===========================================================================

function [datm] = forward(m,freq,nfreq,ang,nang);

%--------------- GLOBAL VARIABLES ------------------------------------------
global c1 rho1 znorm lay_thick sz
global minlim maxlim sd i_nan
%----------------------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

% rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
rhos = m(2) + sin(znorm*pi/2).^m(4)*(m(3)-m(2));
% cs=ct+(cb-ct)*znorm;
cs=m(5)+(m(6)-m(5))*znorm;

geo_sin=[NaN c1 0 rho1; ...
lay_thick*m(1)*ones(sz,1) cs' m(7)*ones(sz,1) rhos'; ...
        NaN m(6) m(7) m(3)];
%
% Compute reflectivity:
%
%freq, ang, geo_sin

[ref] = ref_nlay3(ang, geo_sin,freq);% compute Reflection
%
% Careful with transpose! Matlab is doing conjugate transpose of 
% complex arrays! However, ref is real since we took the ABS...
%
tmp = ref';
datm = reshape(tmp,[nfreq*nang 1]);
datm(i_nan) = 0;
		  
return

%===========================================================================
%  REF_NLAY3.M
%===========================================================================

function [ref] = ref_nlay3(thd, geo,freq)
%    [ref] = ref_nlay2(thd, geo,freq)
%    [ref] = ref_nlay2(1:90, geo_sin,[200 1000 4000])
%
% compute the plane wave reflection coefficient for arbitrary number of isospeed layers
%
% thd is grazing angle in degrees
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
% "freq" is frequencies
% ref is magnitude of the plane wave pressure reflection coefficient (BL =-20 log10 (ref))
%
% Charles W. Holland - June 2003


c=geo(:,2);
alf=geo(:,3);
r=geo(:,4);
d=geo(2:end-1,1);

nfrq=length(freq);
ref=zeros(nfrq,length(thd));

clear i
dB2nep=2*pi*20*1000/log(10);
v=1./(1./c + alf*i/dB2nep); % force radiation condition to be satisfied
nlay=length(c)-1;
th1=thd*pi/180;

zz=zeros(length(c),length(thd));th=zz;
zz(1,:)=r(1)*c(1)./sin(th1); %since incident angles are real so must z1
for m=2:length(c)
  th(m,:)=acos( cos(th1)*v(m)/c(1) );
  zz(m,:)=r(m)*v(m)./sin(th(m,:));
end

zt=repmat(zz,1,nfrq);
z=reshape(zt,size(zz,1),size(zz,2),nfrq);clear zt;

k=(2*pi*repmat(freq,nlay+1,1)./repmat(v,1,nfrq)).';
  
  if nlay==1; % halfspace problem
    rfh=(z(2,:) - z(1,:) )./(z(2,:) + z(1,:));
    ref=repmat(rfh,nfrq,1);
  else
      
%COMPUTE input impedance
 tankd= tan( repmat(k(:,nlay),1,length(thd)).*d(nlay-1).*repmat(sin(th(nlay,:)),nfrq,1 ) );
 znlay=squeeze(z(nlay,:,:));  znlay1=squeeze(z(nlay+1,:,:));
 zin = znlay.*(znlay1-i*znlay.*tankd.')./(znlay - i*znlay1.*tankd.');

 for m=nlay:-1:3
    tankd= tan( repmat(k(:,m-1),1,length(thd)).*d(m-2).*repmat(sin(th(m-1,:)),nfrq,1 )  );
     zm1=squeeze(z(m-1,:,:));
     zin = zm1.*(zin -i*zm1.*tankd.')./(zm1 - i*zin.*tankd.');
 end
 z1=squeeze(z(1,:,:));
 ref=(zin-z1)./(zin+z1);
end 

ref=-20*log10(abs(ref.'));
%ref = abs(ref.');

% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
