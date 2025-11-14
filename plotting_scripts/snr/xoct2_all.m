function [x1] = xoct2_all(x, tsq_lim1, n, l, block_l,duss_sen)
%
% processes and filters time series data
% x is data structure
% tsq_lim are the time limits of integration
% "n" is center frequencies, it is assumed that the bandwidth is simply the mean of the frequency intervals
% if "n" is negative 1/3 octave filters are used
% "l" is parameter of how far back from the inital integration time, the noise is estimated
% "block_l" is the FFT block size
% duss_sen is the system sensitivity
%
% in this version the noise is computed looking over a longer time window
% from teh beginning of the file until 0.01+"l" (el) before the direct arrival


if n<0; n=abs(n); x1.pref = [n', 10.^([n', n'-.5, n'+.5]/10)];
  else 
      bw=[diff(n) n(end)-n(end-1) ];
      x1.pref=[n',n',n'-bw'/2,n'+bw'/2];
end

x1.t = x.t;
x1.r = x.r;
x1.ang = x.ang;
res = x.fout/block_l;
lim = [ceil(x1.pref(:, 3)/res), floor(x1.pref(:, 4)/res)]+1;
w = 10*log10(lim(:, 2)-lim(:, 1)+1);

for i1 = 1:size(tsq_lim1, 2)
     t_lim = sqrt(tsq_lim1(:, i1)+x.r(i1)^2/x.c^2);
     i2 = find(t_lim(1) <= x.t_ax & x.t_ax <= t_lim(2));
     tmp = fft(x.t_s(i1, i2), block_l);
     tmp = tmp.*conj(tmp);
     for i = 1:size(x1.pref, 1)
        if duss_sen==0;
           sen(i) = vsystsen(0, 0, x.eq(i1), (x1.pref(i, 3)+x1.pref(i, 4))/2); %gains have already been accounted for
         else
           sen(i)=duss_sen; 
        end
       x1.oct(i, i1) = 10*log10(2/x.fout^2*sum(tmp(lim(i, 1):lim(i, 2))))...
               -w(i)-sen(i);          
     end
     t = x.t(i1)-l;
     %i3 = find(t-t_lim(2)+t_lim(1) <= x.t_ax & x.t_ax <= t);
     i3 = find(0 <= x.t_ax & x.t_ax <= (t-.01)  );
     tmp = fft(x.t_s(i1, i3), block_l);
     tmp = tmp.*conj(tmp);
     correct = 10*log10(length(i2)/length(i3));
     for i = 1:size(x1.pref, 1)
          x1.SNR(i, i1) = x1.oct(i, i1)-...
               10*log10(2/x.fout^2*sum(tmp(lim(i, 1):lim(i, 2))))...
               +w(i)+sen(i)-correct;
     end
end

