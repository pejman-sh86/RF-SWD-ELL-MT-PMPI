function [C] = cov_est_par_d(ang,sigma,lambda,beta);

nang = length(ang);

for i = 1:nang
    for j = 1:nang
        C(i,j) = sigma^2 * (1 + (ang(i)-ang(j))^2/beta^2)^(-lambda^2);
    end;
end;

return;
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
