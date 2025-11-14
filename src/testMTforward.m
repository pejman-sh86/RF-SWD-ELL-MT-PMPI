I_ZMT = 0;

f = 10.0 .^ [-8.0:0.1:2.0];
%f = fmt';
%f = 10.0 .^ [-4.0:1:4.0];
%sigma = [1.0e-1, 1.0e-3, 1.0e-2];
%sigma = [1.0e-2, 1.0e-2, 1.0e-2, 1.e-2];
%sigma = [1.0e-1, 1.0e-3, 1.0e0];
sigma = [1.e-1, 1.e-2, 1.e-3, 1.e-4];
%sigma = [1.0e-2, 1.0e-2, 1.0e-2, 1.e-2];
%h = [1000, 49000.0];
%h = [5000.0, 195000.0];
h = [40000., 60000., 100000.];
w = 2*pi*f;
u0 = 4*pi*(1.0e-7);
wu0 = w*u0;
%k = zeros(len(sigma), len(f));
%k = sqrt( 1i* (sigma'*wu0) ); %each column is at each frequency for different sigmas
k_real = sqrt( sigma'*wu0  /2.0 ); %each column is at each frequency for different sigmas
k = complex(k_real, k_real);

a = sqrt(sigma(2)/sigma(3)) * tanh( 1*k(3,:)*h(3) + atanh( sqrt(sigma(3)/sigma(4)) ) );
b = sqrt(sigma(1)/sigma(2)) * tanh( 1*k(2,:)*h(2) + atanh(a) );
c = tanh( 1*k(1,:)*h(1) + atanh(b) );
Z = 1i*wu0 ./ k(1,:) .* c;
%Z = conj(Z);

appRes = abs(Z) .^2 ./ wu0;
phase = atan2(imag(Z), real(Z)) * 180.0/pi;
%% 

figure(1)
subplot(2,1,1)
loglog(sqrt(1./f), appRes, '--*r')
%loglog(1./f, appRes, '--*r')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Amplitude (\Omega-m)')
hold on
subplot(2,1,2)
semilogx(sqrt(1./f), phase, '--*r')
%semilogx(1./f, phase, '--*r')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Phase (deg)')
hold on
%% 

pred = load('Athabasca_mappredMT.dat');
data = load("Athabasca_MT.dat");
fre = data(:,1);
NDAT_MT = length(fre);
amp = pred(1:NDAT_MT);
phs = pred(NDAT_MT+1:2*NDAT_MT);

%figure(2)
subplot(2,1,1)
loglog(sqrt(1./fre), amp,'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Amplitude (\Omega-m)')
legend('Analytic', 'my code')
subplot(2,1,2)
semilogx(sqrt(1./fre), phs, 'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Phase (deg)')


%% 

figure(1)
subplot(2,1,1)
%semilogx(sqrt(1./f), real(Z), '--r*')
semilogx(1./f, real(Z), '--r*')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Real(Z)')
hold on
subplot(2,1,2)
%semilogx(sqrt(1./f), imag(Z), '--r*')
semilogx(1./f, imag(Z), '--r*')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Imag(Z)')
hold on

if I_ZMT == 1
pred = load('Athabasca_mappredMT.dat');
data = load("Athabasca_MT.dat");
fre = data(:,1);
NDAT_MT = length(fre);
rZ = pred(1:NDAT_MT);
iZ = pred(NDAT_MT+1:2*NDAT_MT);
%figure(2)
subplot(2,1,1)
%semilogx(sqrt(1./fre), rZ,'b-')
semilogx(1./fre, rZ,'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Real(Z)')
legend('Analytic', 'my code')
subplot(2,1,2)
%semilogx(sqrt(1./fre), iZ, 'b-')
semilogx(1./fre, iZ, 'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Imag(Z)')

end
%%
fre = f;
amp = zeros(1,length(fre));
phs = zeros(1,length(fre));
for n = 1:length(fre)
[amp(n), phs(n)] = modelMT(1./sigma, h, f(n));
end

subplot(2,1,1)
loglog(sqrt(1./fre), amp,'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Amplitude (\Omega-m)')
legend('Analytic', 'my code')
subplot(2,1,2)
semilogx(sqrt(1./fre), phs, 'b-')
xlabel('sqrt(T) (s^{1/2})')
ylabel('Phase (deg)')