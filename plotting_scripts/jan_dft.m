function test;

equal = 0;
N = 101;
fac = 1;

if(equal ~= 1)
  dt = rand(1,N-1);
  dt = sort(dt);
  dt = dt/sum(dt)*fac;
  t(1) = 0;
  for i = 2:N
    t(i) = t(i-1)+dt(i-1);
  end
else
  d = 1/(N-1)*fac;
  t = 0:d:(N-1)*d;
  for i = 1:N-1
    dt(i) = d;
  end
end
figure(10)
plot(t)

x = sin(2*pi*20*t)+sin(2*pi*40*t);

T = t(end) + dt(1);
%df = 1/(N*dt(1));
df = 1/T;
f = [0:df:(N-1)*df];



[XM] = fft(x,N);
[x2] = ifft(XM,N);
x2 = real(x2);
PXM = XM .* conj(XM)/N;
PHM = atan2(imag(XM),real(XM));

[XJ] = jan_dft(x,t,dt,f,df);
[xj] = jan_idft(XJ,t,dt,f,df);
%xj
xj = real(xj);
PXJ = XJ .* conj(XJ);
PHJ = atan2(imag(XJ),real(XJ));

%save bla XM XJ;

figure(1);
subplot(3,1,1);
plot(x);
subplot(3,1,2);
plot(x-x2);
subplot(3,1,3);
plot(xj);

figure(2); 
subplot(2,1,1);
plot(f(1:floor(length(t)/2)+1),PXM(1:floor(length(t)/2)+1));
subplot(2,1,2);
plot(f(1:floor(length(t)/2)+1),PXJ(1:floor(length(t)/2)+1));

figure(3); 
subplot(2,1,1);
plot(f,PHM);
subplot(2,1,2);
plot(f,PHJ);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X] = jan_dft(x,t,dt,f,df)

N = length(t);

clear i;
for k = 1:N 
  X(k) = 0;
  for j = 1:N-1
    X(k) = X(k) + x(j) * exp(-i*2*pi*f(k)*t(j)) * dt(j);
%    X(k) = X(k) + x(j) * exp(-i*2*pi*(k-1)*(j-1)/N) * dt(j);
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = jan_idft(X,t,dt,f,df)

N = length(t);

clear i;
for j = 1:N 
  x(j) = 0;
  for k = 1:N
    x(j) = x(j) + X(k) * exp(i*2*pi*f(k)*t(j));
%    x(j) = x(j) + X(k) * exp(i*2*pi*(k-1)*(j-1)/N);
  end
  x(j) = df * x(j);
%  x(j) = 1/N * x(j);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=ddft(t,x,f)
% function X=dft(t,x,f)
% Compute DFT (Discrete Fourier Transform) at frequencies given
%   in f, given samples x taken at times t:
%     X(f) = sum { x(k) * e**(2*pi*j*t(k)*f) }
%             k

shape = size(f);
t = t(:); % Format 't' into a column vector
x = x(:); % Format 'x' into a column vector
f = f(:); % Format 'f' into a column vector

% It's just this simple:
W = exp(-2*pi*j * f*t');
X = W * x;
X = reshape(X,shape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=iddft(t,x,f)
% function X=dft(t,x,f)
% Compute DFT (Discrete Fourier Transform) at frequencies given
%   in f, given samples x taken at times t:
%     X(f) = sum { x(k) * e**(2*pi*j*t(k)*f) }
%             k

shape = size(f);
t = t(:); % Format 't' into a column vector
x = x(:); % Format 'x' into a column vector
f = f(:); % Format 'f' into a column vector

% It's just this simple:
W = exp(2*pi*j * f*t');
X = W * x;
X = reshape(X,shape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%EOF
