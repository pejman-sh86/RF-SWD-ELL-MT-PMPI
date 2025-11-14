% =====================================================================
%   FORWARD.M
% =====================================================================

function [E,Fm] = forward_hol_grad_assa(m,F);

global icorr

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
  rhos = m(1) + sin(znorm*pi/2).^m(3)*(m(2)-m(1));
  cs=m(4)+(m(5)-m(4))*znorm;

  geo_sin=[NaN F(1).cw 0 F(1).rw; ...
  lay_thick*F(1).h*ones(sz,1) cs' m(6)*ones(sz,1) rhos'; ...
          NaN m(5) m(6) m(2)];
else
  % rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
  rhos = m(2) + sin(znorm*pi/2).^m(4)*(m(3)-m(2));
  % cs=ct+(cb-ct)*znorm;
  cs=m(5)+(m(6)-m(5))*znorm;

  geo_sin=[NaN F(1).cw 0 F(1).rw; ...
  lay_thick*m(1)*ones(sz,1) cs' m(7)*ones(sz,1) rhos'; ...
          NaN m(6) m(7) m(3)];
end

%
% Compute reflectivity:
%
if((m > minlim) & (m < maxlim))
  E = 0;
  for i=1:nfreq
    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
    ref         = -20*log10(abs(ref));
    Fm(i).dat = ref';
%    E = E * sum((Fm(i).dat-F(i).dat).^2);
    if icorr == 1
        E = E + (Fm(i).dat-F(i).dat)'*F(i).c_inv*(Fm(i).dat-F(i).dat)/2;
    else
        E = E + length(F(i).dat)/2*log(sum((Fm(i).dat-F(i).dat).^2));
    end
  end
%  E = E/10^28;
else
  disp('Out of bounds')
  for i=1:nfreq
    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
    ref         = -20*log10(abs(ref));
    Fm(i).dat = ref';
  end
  E = 10^50.;
end

return;
% --------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
