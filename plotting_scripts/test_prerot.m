function [Cr] = test_prerot();

%A = [ 2  3 -1 -1  2  2  3; ...
%     -1 -2 -1 -1  7 -1  5; ...
%      8  4 -1 -1 -5  8  2; ...
%      4  6 -1 -1  7  4  1; ...
%      9  3 -1 -1 -5  9  7; ...
%      7  4 -1 -1  3  7  4; ...
%      6  9 -1 -1 -1  6  9;];
A = [ 2  3 -1 -1  2  2  3; ...
     -1 -2 -1 -1  7 -1  5; ...
      8  4 +1 -1 -5  8 -2; ...
      4  6 -1 -1  7  4  1; ...
     -9  3 -1 -1 -5  9  7; ...
      7  4 -1 -1  3  7  4; ...
      6  9 -1 -1 -1  6  9;];

x = [4 3 2 5 7 3 6]'


b = A*x;
b = b + 1*randn(size(b));

b'

x = pinv(A)*b;

x'

C = .1*.1* pinv((A'*A));
sing = svd(A)


nmod = 7;
for im=1:nmod,
   for jm=1:nmod,
      norm = sqrt(C(im,im)*C(jm,jm));
      if (norm > 0)
         Cr(im,jm) = C(im,jm)/norm;
      else
         Cr(im,jm) = 0;
      end
   end
end

Cr

load fort.77
for i=1:nmod

    fort(:,i) = fort(:,i)-mean(fort(:,i));

end

COV = zeros(nmod,nmod);

for i = 1:nmod
    for j = 1:nmod

        COV(i,j) = sum(fort(:,i).*fort(:,j))/(size(fort,1));

    end
end

for im=1:nmod,
   for jm=1:nmod,
      norm = sqrt(COV(im,im)*COV(jm,jm));
      if (norm > 0)
         Cr2(im,jm) = COV(im,jm)/norm;
      else
         Cr2(im,jm) = 0;
      end
   end
end

Cr2

return;
