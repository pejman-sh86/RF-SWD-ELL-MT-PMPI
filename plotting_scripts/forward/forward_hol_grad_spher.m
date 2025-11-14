% =====================================================================
%   FORWARD.M
% =====================================================================

function [E,Fm] = forward_hol_grad_spher(m,F);

z_h = 140;
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
  E = 1;
  for i=1:nfreq

  thd = F(i).ang;
  d3 = .01;
  thdm = thd-d3;
  thdp = thd+d3;

  ref=ref_nlay3(thd,geo_sin,F(1).freq(i));
  refm=ref_nlay3(thdm,geo_sin,F(1).freq(i));
  refp=ref_nlay3(thdp,geo_sin,F(1).freq(i));
  [thd_deg,spher_r_num]=spherical_ref_Brek_num_b(thd,...
                     thdm,thdp,F(1).cw,ref,refm,refp,F(1).freq(i),z_h);


  spher_r_num = -20*log10(abs(spher_r_num));

  Fm(i).dat = spher_r_num';
  Fm(i).th_deg = thd_deg;
  Fm(i).sph_r = spher_r_num;

  end

return;
% --------------------------------------------------------------------
% ...this is the end my fiend.
% EOF
