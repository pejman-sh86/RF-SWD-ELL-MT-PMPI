% =====================================================================
%   FORWARD.M
% =====================================================================

function [E] = forward_hol_grad_assa_ssl(m,F);

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
rhos = m(2) + sin(znorm*pi/2).^m(4)*(m(3)-m(2));
% cs=ct+(cb-ct)*znorm;
cs=m(5)+(m(6)-m(5))*znorm;

geo_sin=[NaN F(1).cw 0 F(1).rw; ...
lay_thick*m(1)*ones(sz,1) cs' m(7)*ones(sz,1) rhos'; ...
        NaN m(6) m(7) m(3)];
%
% Compute reflectivity:
%
if((m > minlim) & (m < maxlim))
  E = 0;
  for ifreq=1:nfreq

    Fm(ifreq).csave = cov_est_par_a(F(ifreq).ang,m(7+ifreq),...
                      m(7+nfreq+ifreq));
    Fm(ifreq).c_inv = inv(Fm(ifreq).csave);
  
    [ref] = ref_nlay3(F(ifreq).ang,geo_sin,F(1).freq(ifreq));% compute BL
    ref         = -20*log10(abs(ref));
    Fm(ifreq).dat = ref';
  
    E = E + (Fm(ifreq).dat-F(ifreq).dat)'*Fm(ifreq).c_inv*...
        (Fm(ifreq).dat-F(ifreq).dat)/2 + 1/2 * log(det(Fm(ifreq).csave));

  end
%  E = E/10^15;
else
  disp('Out of bounds')
%  for i=1:nfreq
%    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
%    ref         = -20*log10(abs(ref));
%    Fm(i).dat = zeros(length);
%  end
  E = 10^50.;
end

return;
% --------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
