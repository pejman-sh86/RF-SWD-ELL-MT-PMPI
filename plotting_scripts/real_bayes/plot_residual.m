function ref_vs_rep
%
%PLOT reflection model and data 
%
% 1) Uses MAP.mat from plot_hist!
% 2) Uses sd_mean.mat (created by ref_vs_data)
% 3) Saves residuals.dat
%

set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);
global c1 rho1 znorm lay_thick sz

sample1 = 'third_second_sample1.dat';
sample2 = 'third_second_sample2.dat';
plotfile1 = strrep(sample1,'sample1.dat','residual.eps');
plotfile2 = strrep(sample1,'sample1.dat','autocorr.eps');

load blde4_2_ping_ida_inv.mat;
load sd_mean.mat;
load map.mat;
aincr   =  1;		%
a_start =  1;		%
nang    =135;		% Nr of angles
fincr   =  1;		%
f_start =  6;		%
nfreq   =  9;		% Nr of frequencies
ang = zeros(1,nang);
fr = zeros(1,nfreq);
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
m = xml;

%xml(1)


j = a_start;
for i = 1:nang
  ang(i) = xde4.ang(j);
  j = j + aincr;
end
j = f_start;
for i = 1:nfreq
  fr(i) = xde4.pref(j,1);
  j = j + fincr;
end

tmp = xde4.bl(f_start:fincr:f_start+(fincr*nfreq)-1,a_start:aincr:...
              a_start+(aincr*nang)-1);
dat = tmp;

[rep] = forward(m,fr,nfreq,ang,nang);
for i = 1:nfreq
  diff(i,:) = (rep(i,:)-dat(i,:))/sd_mean(i);
end

fid1 = fopen('residuals.dat','w');
for i = 1:nfreq
  fprintf(fid1,'%10.4f',diff(i,:));
  fprintf(fid1,'\n');
end
fclose(fid1);

figure(1);
for k=1:length(fr)
   subplot(3,3,k);hold on;box on
   plot(ang, diff(k,:),'kx');
   hold on;
   
   title([num2str(fr(k)) ' Hz'] )
   axis([0 90 -5 5]); set(gca,'Xtick',[0:30:90]);
   if k>6;xlabel('Angle [deg]');end
   if k == 1;ylabel('Residual [\sigma]');end
   if k == 4;ylabel('Residual [\sigma]');end
   if k == 7;ylabel('Residual [\sigma]');end
   set(gca,'YGrid','on');
end

saveas(gca,plotfile1,'epsc2');

%
% Autocorrelation of residuals
%
figure(2);
for k=1:length(fr)
   j = 1;
   for i=1:nang
     if isnan(diff(k,i))~=1
       tmp2(j) = diff(k,i);
       j = j+1;
     end
   end
   mxlag = floor((j - 1)/2);
   [acor,lag]= xcorr(tmp2,mxlag,'coeff');
   subplot(3,3,k);hold on;box on
   plot(lag, acor,'-xk');
   hold on;
   
   title([num2str(fr(k)) ' Hz'] )
%   axis([0 90 -5 5]); set(gca,'Xtick',[0:30:90]);
   if k>6;xlabel('Lag');end
   set(gca,'YGrid','on');
%
% Perform run test:
%
%   z = runtest(tmp2')
end
saveas(gca,plotfile2,'epsc2');

%===========================================================================
%   FORWARD.M
%===========================================================================

function [datm] = forward(m,freq,nfreq,ang,nang);

%--------------- GLOBAL VARIABLES ------------------------------------------
global c1 rho1 znorm lay_thick sz
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
[ref] = ref_nlay2(ang, geo_sin,freq);% compute Reflection
%
% Careful with transpose! Matlab is doing conjugate transpose of 
% complex arrays! However, ref is real since we took the ABS...
%
datm = ref;
%datm = reshape(tmp,[nfreq*nang 1]);

return

%===========================================================================
%  REF_NLAY3.M
%===========================================================================
function [ref] = ref_nlay2(thd, geo,freq)
%    [ref] = ref_nlay2(thd, geo,freq)
%    [ref] = ref_nlay2(1:90, geo_sin,[200 1000 4000])
%
% compute the plane wave reflection coefficient for arbitrary number 
% of isospeed layers
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
% ref is magnitude of the plane wave pressure reflection coefficient 
% (BL =-20 log10 (ref))
%
% Charles W. Holland - June 2003


c = geo(:,2);		% Velocity profile
alf = geo(:,3);		% Attenuation profile
r = geo(:,4);		% Density profile
d = geo(2:end-1,1);	% Layer thicknesses
nrd = length(c); nra = length(thd); nfrq = length(freq);

clear i			% imaginary unit!
dB2nep=2*pi*20*1000/log(10);
v=1./(1./c + alf*i/dB2nep); % force radition condition to be satisfied
nlay=nrd-1;		% Nr layers is nr depths - 1
th1=thd*pi/180;

z=zeros(nrd,nra);th=z;
%
% First row is not eaqual to loop solution. However, the first
% is not used anyways?
%
th = acos( v*cos(th1)/c(1) );
z = ((repmat((v.*r),1,nra))./sin(th));
z(1,:)=r(1)*c(1)./sin(th1); %since incident angles are real so must z1

% treating k as matrix...
k=2*pi*repmat(freq,nrd,1)./repmat(v,1,nfrq);

if(nlay==1) % halfspace problem
  for n=1:nfrq
    ref(n,:)=(z(2,:) - z(1,:) )./(z(2,:) + z(1,:));
  end
else
  for n=1:nfrq
    %COMPUTE input impedance
    tankd= tan( k(nlay,n)*d(nlay-1)*sin(th(nlay,:)  ) );
    zin = z(nlay,:).*(z(nlay+1,:)-i*z(nlay,:).*tankd)./(z(nlay,:) -...
          i*z(nlay+1,:).*tankd);

    tankd = tan(repmat(k(2:nlay-1,n).* d(1:nlay-2),1,nra) .* sin(th(2:nlay-1,:)));
    for m=nlay:-1:3
      zin = z(m-1,:).*(zin -i*z(m-1,:).*tankd(m-2,:))./(z(m-1,:) -...
      i*zin.*tankd(m-2,:));
    end
    ref(n,:)=(zin-z(1,:))./(zin+z(1,:));
  end
end

ref=-20*log10(abs(ref));

%%plotting paramters
%amx = 90;  atk=15; ymn=0;  ymx=50; ylab=1:3:24; colr='bgrckmbgrckmbgrckm';
%
%if nfrq<15
% subplot(2,2,2)
% for n=1:nfrq  
%   plot(thd,abs(ref(n,:)),colr(n)); hold on 
% end
%   axis([0 amx ymn ymx]); 
%   legend([num2str(freq') repmat(' Hz',size(freq'))])  
%    xlabel('Grazing Angle (deg)');ylabel('|R| (dB)')
%   set(gca,'Position',[.5,0.2,0.4,0.4])
%end 


% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
