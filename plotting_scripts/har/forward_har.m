%===========================================================================
%   FORWARD.M
%===========================================================================

function [E,VB,t1] = forward_har(m,freq,ang,Bft);

%global fstart fend astart aend sd

% 2-Layer model
h1  = [   0 m(1)     ];
c1  = [1500 m(2) m(5)];
r1  = [   1 m(3) m(6)];
a1  = [   0 m(4) m(7)];

%  3-Layer model
%h1  = [   0 m(1) m(5)     ];
%c1  = [1500 m(2) m(6) m(9)];
%r1  = [   1 m(3) m(7) m(10)];
%a1  = [   0 m(4) m(8) m(11)];

[t1,V] = arb_layer(c1,h1,r1,a1,freq);

V2 = fliplr(-20*log10(abs(V)));
t2 = fliplr(round(t1));
ship = 0.;

%
% Simulate smearing of beam forming
%
[RLcalc,angle,VB,newang] = ...
    Smudge(freq,t2,V2,1500,1500,1500,ship,Bft,'yes',t2);
VB = 10.^(-VB/20);
%VB = V2;

%ndat = length(dat);
%datm = zeros(ndat,1);
%id = 1 ;
%for ifreq=fstart:fend
%   for iang=astart:aend
%      if (Weight(ifreq,iang) == 1);
%         datm(id) = VB(ifreq,iang);
%         id = id+1;
%      end
%   end
%end


E = 1;
%E = sum((datm-dat).^2./(sd.*sd*2));


%E = (datm-dat)'*Cd_inv*(datm-dat)/2;
%E = 0;
%for i  = 1:nfreq
%  E = E + ndat(i)/2 * log(sum((datm(nang*(i-1)+1:nang*i)-...
%  dat(nang*(i-1)+1:nang*i)).^2));
%end
%E

return;

