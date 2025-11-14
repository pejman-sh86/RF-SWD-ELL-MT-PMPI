function [alpha] = atten_fit(f);

f = f/1000;
alpha = 3.3e-3 + 0.11*f.^2/(1+f.^2) + 44*f.^2/(4100+f.^2)+...
        3.0e-4*f.^2;
return;
