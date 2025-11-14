%===========================================================================
%   FORWARD.M
%===========================================================================

function [E,Fm] = forward_har_assa(m,F,Bft);

minlim = (F(1).minlim)';
maxlim = (F(1).maxlim)';

if(F(1).nmod == 7)
  % 2-Layer model
  h1  = [   0 m(1)     ];
  c1  = [1500 m(2) m(5)];
  r1  = [   1 m(3) m(6)];
  a1  = [   0 m(4) m(7)];
elseif(F(1).nmod == 11)
  %  3-Layer model
  h1  = [   0 m(1) m(5)     ];
  c1  = [1500 m(2) m(6) m(9)];
  r1  = [   1 m(3) m(7) m(10)];
  a1  = [   0 m(4) m(8) m(11)];
end

freq = F(1).freq;
ang = F(1).ang;
if F(1).ave == 0
  spacing = [0];
else
  spacing = [-23.43 -11.72 0 11.72 23.43];
end

for ifreq = 1:length(freq)

  fr = spacing + freq(ifreq);

  [t1,V] = arb_layer_new(c1,h1,r1,a1,fr,ang);

  V2 = fliplr(-20*log10(abs(V)));
  t2 = fliplr(round(t1));
%
% RMS average
%
  V2 = V2.^2;
  V3(ifreq,:) = sqrt(mean(V2(1:length(spacing),:),1));
end 

%
% Simulate smearing of beam forming
%
ship = 0.;
[RLcalc,angle,VB,newang] = Smudge(freq,t2,V3,1500,1500,1500,ship,Bft,'yes',t2);

if((m > minlim) & (m < maxlim))
%  E = sum((datm-dat).^2./(2*sd*sd));
  E = 0;
  for ifreq=1:length(freq)
    Fm(ifreq).dat = (10.^(-VB(ifreq,:)/20))';
  
    E = E + sum((Fm(ifreq).dat(F(1).astart:F(1).aend)-...
        F(ifreq).dat(F(1).astart:F(1).aend)).^2./(F(1).sd.*F(1).sd*2));
%    E = E + (Fm(ifreq).dat(F(1).astart:F(1).aend)-F(ifreq).dat(F(1).astart:...
%        F(1).aend))' * F(ifreq).c_inv * (Fm(ifreq).dat(F(1).astart:F(1).aend)-...
%        F(ifreq).dat(F(1).astart:F(1).aend))/2;
%    E = E * sum((Fm(i).dat-F(i).dat).^2);
  end
else
  disp('Out of bounds')
  for ifreq=1:length(freq)
    Fm(ifreq).dat = 10.^(-VB(ifreq,:)/20);
  end
  E = 10^50.;
end

return;

