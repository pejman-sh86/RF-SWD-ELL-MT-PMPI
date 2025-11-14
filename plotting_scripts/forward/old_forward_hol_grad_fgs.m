% =====================================================================
%   FORWARD.M
% =====================================================================
function [E,Fm] = old_forward_hol_grad_fgs(m,F);

minlim = F(1).minlim;
maxlim = F(1).maxlim;
lay_thick = F(1).lay_thick;
znorm = F(1).znorm;
sz = F(1).sz;
nfreq = length(F(1).ifreq);

%----------------------------------------------------------------
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
E = 0;
for ifreq=1:nfreq
  [ref] = ref_nlay3(F(ifreq).ang,geo_sin,F(1).freq(ifreq));% compute Reflection
  Fm(ifreq).dat = ref';

  E = E + (Fm(ifreq).dat-F(ifreq).dat)'*F(ifreq).c_inv*(Fm(ifreq).dat-F(ifreq).dat)/2;
end

return
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF

