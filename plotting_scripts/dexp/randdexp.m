function [x] = randdexp(m,n,sd);

  x = random('Exponential',sd,m,n);
  y = sd * randn(m,n);
  
  for i = 1:m
    for j = 1:n
      if(y(i,j) < 0);
        x(i,j) = -1*x(i,j);
      end
    end
  end

  clear y;
return;
