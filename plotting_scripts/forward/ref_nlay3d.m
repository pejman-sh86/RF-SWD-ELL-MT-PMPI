function [ref] = ref_nlay3c(thd, geo,freq)
%    [ref] = ref_nlay3(thd, geo,freq)
%    [ref] = ref_nlay3(1:90, geo_sin,[200 1000 4000])
%
% compute the plane wave reflection coefficient for arbitrary number of isospeed layers
%  after Brekovskikh/Godin Waves in Layered Media I Eq 2.5.1-2.5.4  
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
% ref is the (complex) plane wave pressure reflection coefficient (BL =-20 log10 (ref))
%
% Charles W. Holland - June 2003
% eliminated do loop on frequency by adding another dimension to several matrices 
warning off
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
z=reshape(zt,size(zz,1),size(zz,2),nfrq);
clear zt;

squeeze(z(1,:,:))'

for kf=1:nfrq

k=(2*pi*repmat(freq,nlay+1,1)./repmat(v,1,nfrq)).';
%repmat(freq,nlay+1,1)
%repmat(v,1,nfrq)

if nlay==1; % halfspace problem
   rfh=(z(2,:,:) - z(1,:,:) )./(z(2,:,:) + z(1,:,:));
   if nfrq==1;ref=rfh.'; else ref=squeeze(rfh);end
%   rfh=(z(2,:) - z(1,:) )./(z(2,:) + z(1,:)); ref=repmat(rfh,nfrq,1);
else
      
%COMPUTE input impedance
 tankd= tan( repmat(k(:,nlay),1,length(thd)).*d(nlay-1).*repmat(sin(th(nlay,:)),nfrq,1 ) );
 znlay=reshape(z(nlay,:,:),length(thd),nfrq);  
 znlay1=reshape(z(nlay+1,:,:),length(thd),nfrq);
 zin = znlay.*(znlay1-i*znlay.*tankd.')./(znlay - i*znlay1.*tankd.');

size(znlay)
size(znlay1)
size(zin)

 for m=nlay:-1:3
    tankd= tan( repmat(k(:,m-1),1,length(thd)).*d(m-2).*...
           repmat(sin(th(m-1,:)),nfrq,1 )  );
     zm1=reshape(z(m-1,:,:),length(thd),nfrq);
     zin = zm1.*(zin -i*zm1.*tankd.')./(zm1 - i*zin.*tankd.');
 end
 z1=reshape(z(1,:,:),length(thd),nfrq);
 ref=(zin-z1)./(zin+z1);
end 

ref=ref.';
%at 0 degrees, the reflection coefficeint is -1
jj=find(thd==0);ref(jj)=-1;
return;
