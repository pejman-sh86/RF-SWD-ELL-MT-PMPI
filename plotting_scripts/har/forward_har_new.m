%===========================================================================
%   FORWARD.M
%===========================================================================

function [E,Fm] = forward_har_new(m,F,Bft);

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
%spacing = [-50 -33.33 -16.67 0 16.67 33.33 50];
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
  ship = 0.;

%
% RMS average
%
  V2 = V2.^2;
  V3(ifreq,:) = sqrt(mean(V2(1:length(fr),:),1));

end 

%
% Simulate smearing of beam forming
%
[RLcalc,angle,VB,newang] = Smudge(freq,t2,V3,1500,1500,1500,ship,Bft,'yes',t2);

E = 0;
for ifreq=1:length(freq)
  Fm(ifreq).dat = 10.^(-VB(ifreq,:)/20);
  E = E + sum((Fm(ifreq).dat(F(1).astart:F(1).aend)-...
      F(ifreq).dat(F(1).astart:F(1).aend)).^2./(F(1).sd.*F(1).sd*2));
%  E = E + (Fm(i).dat-F(i).dat)'*F(i).c_inv*(Fm(i).dat-F(i).dat)/2;
end

return;

