function [th_deg,spher_r_num]=spherical_ref_Brek_num_b(thd,thdm,thdp,...
                              c, ref,refm,refp,f,z_h);
% [spher_num]=spherical_ref_Brek_num(1:90,1511, ref2,[100 200 400 800 1600 3200], 120);
%
% Computes spherical reflection coefficient
%  from Brekovskikh, Waves in Layered Media,  Eq 28.13 
%  given pressure reflection coeffiecnt "ref" and grazing angles "thd" (degrees)
%  c is the sound speed in the upper halfspace
%  at frequency "f" (Hz ) and source+receiver height "z_h" (m)
%
% ref is computed from
%
%  ref=ref=ref_nlay3(1:90,geo,[100 200 400 800 1600 3200]);
%
% NB:
%   - "ref" can be a single row and f a vector (if a halfspace problem)
%      but generally the rows of "ref" correspond to the frequencies in "f"
%   - the derivatives are computed numerically so the sampling in "thd" is assumed adequate
%
% coded by Charles W. Holland June 2003
% Modified by Jan Dettmer for uneven spacing Jun 2005

th1=thd*pi/180;
th1m=thdm*pi/180;
th1p=thdp*pi/180;


del_rad_b=(thdp-thdm)*pi/180;

k=2*pi*f/c;
R= z_h./sin(th1);

if length(f)>1
  del_rad_b = repmat(del_rad_b,length(f),1);
  disp('WARNING: same ang. accross freq!\n')
end

  derv_b=(refp-refm)./del_rad_b;
  derv2_b=(refp+refm-2.*ref)./(del_rad_b./2).^2; 
  th=repmat(th1,size(ref,1),1); 
  N_b = (derv2_b + derv_b.*tan(th) )./2;
  th_deg=thd;

for jk=1:length(f); 

   if min(size(ref))==1
     spher_r_num(jk,:)=ref - i * N_b./(k(jk)*R);  
   else 
     spher_r_num(jk,:)=ref(jk,:) - i * N_b(jk,:)./(k(jk)*R);
   end
  
end

return;
