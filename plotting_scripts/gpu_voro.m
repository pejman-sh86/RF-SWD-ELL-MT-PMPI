function []=gpu_voro();


dx = x_evn(ir)-xn;
dxsq = dx.*dx;
dz = z_evn(iz)-zn;
dzsq = dz.*dz;
if(idx_bott(ir) < iz);
  d  = sqrt(dxsq + dzsq);
  [dmin,iv] = min(d);
  c_ev(iz,ir) = c(iv);
  r_ev(iz,ir) = r(iv);
  a_ev(iz,ir) = a(iv);
end;




return;
