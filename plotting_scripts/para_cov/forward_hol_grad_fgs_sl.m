% =====================================================================
%   FORWARD.M
%   Forward model for Holland data to invert for data correlations
%   Correlated errors are treated with parametrized approach
%   sigma and lambda are inversion parameters
% =====================================================================
function [E] = forward_hol_grad_fgs_sl(m,F);

minlim = F(1).minlim;
maxlim = F(1).maxlim;
lay_thick = F(1).lay_thick;
znorm = F(1).znorm;
sz = F(1).sz;
nfreq = length(F(1).freq);

% ---------------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

% rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
rhos = F(1).msim(2) + sin(znorm*pi/2).^F(1).msim(4)*(F(1).msim(3)-F(1).msim(2));
% cs=ct+(cb-ct)*znorm;
cs=F(1).msim(5)+(F(1).msim(6)-F(1).msim(5))*znorm;

geo_sin=[NaN F(1).cw 0 F(1).rw; ...
lay_thick*F(1).msim(1)*ones(sz,1) cs' F(1).msim(7)*ones(sz,1) rhos'; ...
        NaN F(1).msim(6) F(1).msim(7) F(1).msim(3)];
%
% Compute reflectivity:
%
  E = 0;
  for ifreq=1:nfreq

    F(ifreq).csave = cov_est_par_a(F(ifreq).ang,m(ifreq),...
                      m(nfreq+ifreq));
    F(ifreq).c_inv = inv(F(ifreq).csave);
  
    [ref] = ref_nlay3(F(ifreq).ang,geo_sin,F(1).freq(ifreq));% compute BL
    ref         = -20*log10(abs(ref));
    Fm(ifreq).dat = ref';
  
    E = E + (F(ifreq).dat-Fm(ifreq).dat)'*F(ifreq).c_inv*...
        (F(ifreq).dat-Fm(ifreq).dat)/2 + 1/2 * log(det(F(ifreq).csave));

  end

return;
% --------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
