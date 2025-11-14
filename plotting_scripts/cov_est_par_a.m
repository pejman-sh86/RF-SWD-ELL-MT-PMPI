function [C] = cov_est_par_a(ang,sigma,lambda);

nang = length(ang);

for i = 1:nang
    for j = 1:nang
        C(i,j) = sigma^2 * exp(-lambda * abs(ang(i)-ang(j)));
    end;
end;

return;
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
