function [stat_fgs]=plot_residuals(nrbin);
%
% Plot residuals for ML model against 
%             
global c1 rho1 znorm lay_thick sz


load sd_mean.mat;
load blde4_2_ping_ida_inv.mat;
aincr   =  3;		%
a_start =  1;		%
nang    = 45;		% Nr of angles
fincr   =  1;		%
f_start =  6;	%
nfreq   =  9;		% Nr of frequencies
load map.mat;
m = xml;
%
% Forward model parameters:
%
c1 = 1511; rho1 = 1.029;
lay_thick = 0.1;    % discretization for layer thickness - controls 
                      % upper frequency limit
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
sz=length(znorm);
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
dat = tmp;
%dat = reshape(tmp',nfreq*nang,1);
[i_nan] = find(isnan(dat)==1);
dat(i_nan) = 0;

[rep] = forward(m,freq,nfreq,ang,nang);
rep(i_nan) = 0;
for i = 1:nfreq
  diff(i,:) = (rep(i,:)-dat(i,:))/sd_mean(i);
end
diff = reshape(diff,nfreq*nang,1);

j=1;
for i = 1:nfreq*nang
  if(diff(i) ~= 0)
    difff(j) = diff(i);
    j = j+1;
  end
end

ddiff = max(difff) - min(difff);
for i =1:nrbin
  edges1(i) = min(difff) + ddiff/nrbin * i;
end
t = -4:0.1:4;
size(t)
realgauss = 1/sqrt(2*pi)*exp(-t.^2/2);

figure(1);
col = ['r' 'b'];
n1 = histc(difff,edges1);
area = sum(n1 * ddiff/nrbin)
n1 = n1/area;
stairs(edges1,n1,col(1));
hold on;
plot(t,realgauss,'b');
xlabel('Residual [\sigma]');
legend('Data residuals','G(0,1)');

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
%datm = reshape(ref',[nfreq*nang 1]);

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
