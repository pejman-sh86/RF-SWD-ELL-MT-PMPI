function [th_deg,spher_r_num]=spherical_ref_Brek_num(thd,c, ref,f,z_h);
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

th1=thd*pi/180;
%del_rad=diff(thd(1:2))*pi/180;
del_rad=diff(thd)*pi/180;
k=2*pi*f/c;
R= z_h./sin(th1(3:end));

if length(f)>1
  del_rad = repmat(del_rad,length(f),1);
  disp('WARNING: same ang. accross freq!\n')
end

  derv=diff(ref,1,2)./del_rad;
  del_rad = del_rad(:,2:end);
  derv2=diff(derv,1,2)./del_rad; 
  th=repmat(th1(3:end),size(ref,1),1); 
  N = (derv2 + derv(:,2:end).*tan(th) )./2;
  th_deg=thd(3:end);
  
for jk=1:length(f); 
  if min(size(ref))==1
     spher_r_num(jk,:)=ref(3:end) - i * N./(k(jk)*R);  
   else 
     spher_r_num(jk,:)=ref(jk,3:end) - i * N(jk,:)./(k(jk)*R);
end
  
%spher_r_num = -20*log10(abs(spher_r_num(jk,:)));

%%if length(f)<2; figure 
% if length(f)<15; figure 
%  plot(th_deg,-20*log10(abs(spher_r_num(jk,:))),'r:')
%  xlabel('Grazing Angle (deg)'); ylabel('Reflection Loss (dB)')
%  title([num2str(f(jk)) ' Hz'])
%  set(gca,'Xlim',[0 90],'XTick',[0:30:90])
% end

end

return;

%if 1==2
%figure
% for jk=1:length(f); subplot(2,3,jk)
%  plot(thd,-20*log10(abs(ref)),'k');hold on
%  plot(th_deg,-20*log10(abs(spher_r_num(jk,:))),'r:')
%  xlabel('Grazing Angle (deg)'); ylabel('Reflection Coefficient')
% end
% legend('pw','spher')
%end 
