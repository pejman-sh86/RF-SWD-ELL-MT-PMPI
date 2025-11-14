function [A,window] = mattaper(X,tf,w,s)
% MATTAPER  Tapering of 2D matrices.
% [A,WINDOW] = MATTAPER(X,TF,W,S) performs gradual tapering by applying 
% either a student's t cumulative distribution function or a cosine 
% function to the four borders of the 2D matrix X. 
% For cosine: W is width of cosine in percent of width of X. S is unused.
% For t-distribution: S defines the (integer) rate of decreasing border 
% values, and W the (integer) position of the half height of the applied 
% border window.
% 
%    tf =
% 1: t-distribution
% 2: cosine
% 
% Example: 
% X = ones(60); % Uniform distribution 
% figure,surf(mattaper(X))
% 
% Written by Christian Lundmand Jensen
% Date: 15th of February, 2013.
% 
% Update: 4th of April 2013
%   -  Added cosine taper

A = zeros(size(X));
if nargin == 1
    tf = 2;
    s = 3;
    w = .2;
end
if nargin == 2
    s = 3;
    w = .2;
end
if nargin == 3
    s = 3;
end

[n,m,p] = size(X);
x = (-w:1:w)'; 
range = length(x); % Problems will arise if X is smaller than 2*range
if n < 2*range+1 || m < 2*range+1
    fprintf('Warning! To avoid asymmetrical borders the first two dimensions of X should not be smaller than 4*w+3. \nw = %g -> 4*w+3 = %g. Size(X) = [%g,%g,%g].',w,4*w+3,n,m,p)
%     return
end

%% t-distribution
if tf == 1
    t = tcdf(x,s); % t-distribution
    t(1) = 0; % Forcing zero at the end

    vert1D = ones(n,1);
    vert1D(1:range) = t;
    vert1D(n-range+1:n) = flipud(t); % 1D
    vert1Drep = repmat(vert1D,1,m); % Expanding into 2D
    horz1D = ones(m,1);
    horz1D(1:range) = t;
    horz1D(m-range+1:m) = flipud(t); % 1D
    horz1Drep = repmat(horz1D,1,n); % Expanding into 2D
    window = vert1Drep.*rot90(horz1Drep);
    
%% cosine taper
elseif tf == 2
    window = repmat(tukeywin(n,w),1,m).*rot90(repmat(tukeywin(m,w),1,n));
end

for i = 1:p
    A(:,:,i) = X(:,:,i).*window;
end

end