% =====================================================================
%   FORWARD.M
% =====================================================================

function [E,Fm] = forward_hol_lingrad_assa(m,F);

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
if(isfield(F,'h'))
  rhos = m(1) + (m(2)-m(1))*znorm;
  cs   = m(3) + (m(4)-m(3))*znorm;

  geo_sin=[NaN F(1).cw 0 F(1).rw; ...
  lay_thick*F(1).h*ones(sz,1) cs' m(5)*ones(sz,1) rhos'; ...
          NaN m(4) m(5) m(2)];
else
  rhos = m(2) + (m(3)-m(2))*znorm;
  cs   = m(4) + (m(5)-m(4))*znorm;

  geo_sin=[NaN F(1).cw 0 F(1).rw; ...
  lay_thick*m(1)*ones(sz,1) cs' m(6)*ones(sz,1) rhos'; ...
          NaN m(5) m(6) m(3)];
end
%
% Compute reflectivity:
%
  E = 1;
  for i=1:nfreq
    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
    ref         = -20*log10(abs(ref));
    Fm(i).dat = ref';
  end

return;
% --------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
