function r=acf(x,normflg)
%
% r=acf(x,normflg)
%
%  This function computes the acf r(k) k=0:length(x)
%  CAUTION the first value r(1) is hat{gamma}(0) IE the first lag is 0
% x = time series vector (column)
% normflg 0 to divide by nr=length(x)
%         1 to divide by nr*sample variance
%         2 to divide by nr-k
%         3 to divide by (nr-k)*sample variance

[nr,nc]=size(x);
if nc > nr
    error('x must be a column vector');
end
r=conv(flipud(x),x);
r=r(nr:end);
if normflg==0
    r=r/nr;
elseif normflg==1
    r=r/nr;
    r=r/r(1);
elseif normflg==2
    den=[nr:-1:1]';
    r=r./den;
elseif normflg==3
    den=[nr:-1:1]';
    r=r./den;
    r=r/r(1);
end
